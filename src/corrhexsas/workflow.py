"""Configuration, shared-structure fitting, and auditable local run outputs."""
from __future__ import annotations
from dataclasses import asdict, replace
from pathlib import Path
from datetime import datetime, timezone
import json
import platform
import subprocess
import warnings
import numpy as np
import scipy
from scipy.optimize import least_squares, lsq_linear
from threadpoolctl import threadpool_limits
from .model import Structure, Numerics, Evaluator
from .profiles import Disk, Shells, TabulatedProfile, SmoothLayers, conserved_solvents
from .data import load_data, sha256


STRUCTURE_NAMES = {"a", "sigma_over_a", "ell_b_over_a", "ell_z_over_a",
                   "grain_radius_over_a", "grain_length_over_a"}
CHEMICAL_NAMES = ("water_position", "core_position", "oil_position", "outer_position",
                  "inner_width_position", "outer_width_position")
PARAMETER_NAMES = STRUCTURE_NAMES | set(CHEMICAL_NAMES) | {"radius_over_a"}


def _resolve(value, parameters):
    if isinstance(value, str):
        if value not in parameters:
            raise ValueError(f"Unknown parameter reference {value!r}")
        return parameters[value]
    if isinstance(value, list):
        return [_resolve(v, parameters) for v in value]
    return value


def structure_from_parameters(p, config):
    a = p["a"]
    axial = p["ell_z_over_a"] if config.get("axially_constant", False) is False else np.inf
    return Structure(a, p["sigma_over_a"]*a, p["ell_b_over_a"]*a,
                     axial*a, p["grain_radius_over_a"]*a,
                     p["grain_length_over_a"]*a, config.get("grain_radius_spread", 0.0))


def profile_from_config(spec, p, base):
    kind = spec["type"]
    if kind == "disk":
        return Disk(p["radius_over_a"]*p["a"], _resolve(spec.get("sld", 1.0), p), _resolve(spec.get("medium", 0.0), p)), {}
    if kind == "shells":
        return Shells(_resolve(spec["radii"], p), _resolve(spec["sld"], p), _resolve(spec.get("medium", 0.0), p)), {}
    if kind == "smooth_layers":
        return SmoothLayers(_resolve(spec["starts"], p), _resolve(spec["widths"], p), _resolve(spec["levels"], p)), {}
    if kind == "tabulated":
        path = base/spec["file"]
        lines = [s.strip() for s in path.read_text(encoding="utf-8-sig").splitlines() if s.strip() and not s.lstrip().startswith("#")]
        if not lines:
            raise ValueError("Empty radial profile file")
        delimiter = "," if "," in lines[0] else None
        first = lines[0].replace(",", " ").split()[0]
        try:
            float(first)
            header = 0
        except ValueError:
            header = 1
        from io import StringIO
        table = np.loadtxt(StringIO("\n".join(lines)), delimiter=delimiter, skiprows=header, ndmin=2)
        if table.shape[1] < 2:
            raise ValueError("Radial file needs r and at least one SLD column")
        return TabulatedProfile(table[:, 0], table[:, 1:].T, _resolve(spec.get("medium", table[-1, 1:].tolist()), p)), {}
    if kind == "conserved_solvents":
        options = {k: _resolve(v, p) for k, v in spec.items() if k != "type"}
        return conserved_solvents(p["a"], [p[k] for k in CHEMICAL_NAMES], **options)
    raise ValueError(f"Unknown profile type {kind!r}")


def _bounds(values):
    if len(values) != 2:
        raise ValueError("A bound must contain a lower and upper value")
    lo, hi = (-np.inf if values[0] is None else float(values[0])), (np.inf if values[1] is None else float(values[1]))
    if not lo < hi:
        raise ValueError("Bounds must increase strictly")
    return lo, hi


class Problem:
    def __init__(self, config_path, *, numerical_override=None):
        self.path = Path(config_path).resolve()
        self.base = self.path.parent
        self.config = cfg = json.loads(self.path.read_text(encoding="utf-8"))
        self.parameters = cfg["parameters"]
        def strings(obj):
            if isinstance(obj, str):
                return {obj}
            if isinstance(obj, list):
                return set().union(*(strings(v) for v in obj)) if obj else set()
            if isinstance(obj, dict):
                return strings(list(obj.values()))
            return set()
        profile_references = strings(cfg.get("profiles", cfg.get("profile", {})))
        unknown = set(self.parameters)-PARAMETER_NAMES-profile_references
        if unknown:
            raise ValueError(f"Unknown parameter names {sorted(unknown)}")
        if not STRUCTURE_NAMES <= set(self.parameters):
            raise ValueError(f"Missing structure parameters {sorted(STRUCTURE_NAMES-set(self.parameters))}")
        self.initial = {k: float(v["value"]) for k, v in self.parameters.items()}
        if not np.all(np.isfinite(list(self.initial.values()))):
            raise ValueError("Initial parameter values must be finite; use axially_constant for infinite ell_z")
        self.free = [k for k, v in self.parameters.items() if v.get("vary", False)]
        self.limits = {}
        for key in self.free:
            lo, hi = _bounds(self.parameters[key]["bounds"])
            if not np.isfinite(lo+hi) or not lo <= self.initial[key] <= hi:
                raise ValueError(f"{key} needs finite bounds containing its initial value")
            if self.parameters[key].get("log", False) and lo <= 0:
                raise ValueError(f"Log bounds for {key} must be positive")
            self.limits[key] = (lo, hi)
        if cfg.get("axially_constant", False) and "ell_z_over_a" in self.free:
            raise ValueError("ell_z cannot vary in an axially_constant model")
        if self.initial["sigma_over_a"] == 0 and "sigma_over_a" not in self.free and set(self.free) & {"ell_b_over_a", "ell_z_over_a"}:
            raise ValueError("Correlation lengths are inactive when sigma is fixed at zero")
        self.data, self.inputs = [], {self.path.name: sha256(self.path)}
        for i, ds in enumerate(cfg["datasets"]):
            resolution = ds.get("resolution", cfg.get("resolution", {}))
            path = self.base/ds["file"]
            d = load_data(path, q_unit=ds.get("q_unit", cfg.get("q_unit", "1/A")),
                          relative_error=ds.get("relative_error"), label=ds.get("label", str(i)),
                          resolution_relative=resolution.get("relative", 0.0),
                          resolution_absolute=resolution.get("absolute", 0.0))
            self.data.append(d)
            self.inputs[f"dataset_{i}:{path.name}"] = d.source_hash
        if not self.data:
            raise ValueError("At least one dataset is required")
        self.specs = cfg.get("profiles", [cfg.get("profile")]*len(self.data))
        if len(self.specs) != len(self.data) or any(s is None for s in self.specs):
            raise ValueError("Supply one profile or one profiles entry per dataset")
        for i, spec in enumerate(self.specs):
            if spec["type"] == "tabulated":
                path = self.base/spec["file"]
                self.inputs[f"profile_{i}:{path.name}"] = sha256(path)
        self.numerics = numerical_override or Numerics(**cfg.get("numerics", {}))
        cache, self.engines = {}, []
        for d in self.data:
            key = (d.q.tobytes(), d.sigma_q.tobytes())
            if key not in cache:
                cache[key] = Evaluator(d.q, d.sigma_q, self.numerics)
            self.engines.append(cache[key])
        self.scales = [np.ones_like(d.q) for d in self.data]
        self.regions = []
        self._set_objective()

    def decode(self, x):
        out = dict(self.initial)
        for key, t in zip(self.free, x):
            lo, hi = self.limits[key]
            out[key] = float(np.exp(np.log(lo)+t*np.log(hi/lo)) if self.parameters[key].get("log", False) else lo+t*(hi-lo))
        return out

    def encode(self):
        out = []
        for key in self.free:
            lo, hi = self.limits[key]
            value = self.initial[key]
            out.append(np.log(value/lo)/np.log(hi/lo) if self.parameters[key].get("log", False) else (value-lo)/(hi-lo))
        return np.asarray(out)

    def raw(self, p, *, part="total"):
        structure = structure_from_parameters(p, self.config)
        curves, info = [], []
        for i, (spec, engine) in enumerate(zip(self.specs, self.engines)):
            profile, decoded = profile_from_config(spec, p, self.base)
            result = engine.evaluate(structure, profile, part=part)
            if len(result) not in (1, len(self.data)):
                raise ValueError("A multi-contrast profile must provide one SLD per dataset")
            curves.append(result[0 if len(result) == 1 else i])
            info.append(decoded)
        return curves, info

    def project(self, curves):
        counts = [len(d.q) for d in self.data]
        size, n = sum(counts), len(counts)
        matrix = np.zeros((size, n+1))
        matrix[:, 0] = np.concatenate(curves)
        start = 0
        for i, count in enumerate(counts):
            matrix[start:start+count, i+1] = 1
            start += count
        y = np.concatenate([d.intensity for d in self.data])
        weight = np.concatenate([s/d.uncertainty for s, d in zip(self.scales, self.data)])
        linear = self.config.get("linear", {})
        slo, shi = _bounds(linear.get("scale_bounds", [0, None]))
        blo, bhi = _bounds(linear.get("background_bounds", [0, None]))
        result = lsq_linear(matrix*weight[:, None], y*weight, bounds=(np.r_[slo, np.full(n, blo)], np.r_[shi, np.full(n, bhi)]), tol=1e-11)
        if not result.success:
            raise RuntimeError(f"Linear scale/background fit failed: {result.message}")
        prediction = result.x[0]*np.concatenate(curves)+np.repeat(result.x[1:], counts)
        return np.split(prediction, np.cumsum(counts)[:-1]), result.x

    def _set_objective(self):
        fit = self.config.get("fit", {})
        objective = fit.get("objective", "chi2")
        if objective == "chi2":
            return
        if objective != "balanced":
            raise ValueError("fit.objective must be chi2 or balanced")
        edges = np.asarray(fit.get("region_edges", []), float)
        if edges.ndim != 1 or np.any(np.diff(edges) <= 0) or not np.all(np.isfinite(edges)):
            raise ValueError("region_edges must be finite and strictly increasing")
        baseline, _ = self.project(self.raw(self.initial)[0])
        for i, (d, ref) in enumerate(zip(self.data, baseline)):
            groups = np.digitize(d.q, edges)
            for j in range(len(edges)+1):
                mask = groups == j
                count = int(mask.sum())
                if not count:
                    raise ValueError("Every balanced region must contain data in every dataset")
                normal = max(1.0, float(np.mean(((ref[mask]-d.intensity[mask])/d.uncertainty[mask])**2))) if fit.get("normalize_to_initial", False) else 1.0
                self.scales[i][mask] = 1/np.sqrt(len(self.data)*(len(edges)+1)*count*normal)
                self.regions.append(dict(dataset=i, region=j, count=count, normalization=normal))

    def residual(self, x):
        prediction, _ = self.project(self.raw(self.decode(x))[0])
        return np.concatenate([(pred-d.intensity)/d.uncertainty*s for pred, d, s in zip(prediction, self.data, self.scales)])


def fit_problem(problem):
    settings = problem.config.get("fit", {})
    starts = int(settings.get("starts", 1))
    if starts < 1:
        raise ValueError("fit.starts must be positive")
    x0 = problem.encode()
    rng = np.random.default_rng(settings.get("seed", 20261110))
    trials, best = [], None
    for i in range(starts if len(x0) else 1):
        start = np.clip(x0 if i == 0 else x0+rng.normal(0, 0.1, len(x0)), 1e-9, 1-1e-9)
        if len(x0):
            res = least_squares(problem.residual, start, bounds=(np.zeros(len(x0)), np.ones(len(x0))),
                                x_scale="jac", diff_step=settings.get("diff_step", 1e-4),
                                max_nfev=settings.get("max_nfev", 200), ftol=1e-8, xtol=1e-8, gtol=1e-8)
            trial = dict(x=res.x.tolist(), objective=float(res.fun@res.fun), success=bool(res.success),
                         nfev=res.nfev, message=res.message, singular_values=np.linalg.svd(res.jac, compute_uv=False).tolist())
        else:
            residual = problem.residual([])
            trial = dict(x=[], objective=float(residual@residual), success=True, nfev=1,
                         message="No nonlinear coordinates selected", singular_values=[])
        trials.append(trial)
        if best is None or trial["objective"] < best["objective"]:
            best = trial
    p = problem.decode(best["x"])
    raw, info = problem.raw(p)
    prediction, linear = problem.project(raw)
    chi2 = float(sum(np.sum(((v-d.intensity)/d.uncertainty)**2) for v, d in zip(prediction, problem.data)))
    dof = sum(len(d.q) for d in problem.data)-len(problem.free)-len(linear)
    return dict(parameters=p, scale=float(linear[0]), backgrounds=linear[1:].tolist(),
                objective_name=settings.get("objective", "chi2"), objective=best["objective"],
                chi2=chi2, nominal_dof=dof, nominal_reduced_chi2=chi2/dof if dof > 0 else None,
                success=best["success"], trials=trials, decoded_profiles=info,
                free_parameters=problem.free, frozen_region_normalizations=problem.regions,
                inference_note="Local optimization only. Balanced loss is not a likelihood. No confidence intervals or identifiability are inferred from convergence."), prediction


def write_run(problem, result, prediction, output, *, kind="fit"):
    from . import __version__
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    (out/"config.json").write_text(json.dumps(problem.config, indent=2)+"\n", encoding="utf-8", newline="\n")
    result["numerics"] = asdict(problem.numerics)
    result["kind"] = kind
    (out/"result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    for i, (d, fit) in enumerate(zip(problem.data, prediction)):
        np.savetxt(out/f"curve_{i:02d}.csv", np.c_[d.q, d.intensity, d.uncertainty, d.sigma_q, fit, (fit-d.intensity)/d.uncertainty],
                   delimiter=",", header="q_Ainv,I_cm_inv,dI_cm_inv,sigma_q_Ainv,calculated_cm_inv,standardized_residual", comments="")
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), version=__version__,
                      environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                      inputs=problem.inputs, outputs={p.name: sha256(p) for p in out.iterdir() if p.is_file()})
    try:
        repo = Path(__file__).resolve().parents[2]
        provenance["git_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True, stderr=subprocess.DEVNULL).strip()
        provenance["git_dirty"] = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repo, text=True, stderr=subprocess.DEVNULL).strip())
    except (OSError, subprocess.CalledProcessError):
        provenance["git_commit"] = None
    (out/"provenance.json").write_text(json.dumps(provenance, indent=2)+"\n", encoding="utf-8", newline="\n")
    return out


def plot_run(problem, prediction, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ValueError("Plotting needs pip install 'corrhexsas[plot]'") from exc
    fig, axes = plt.subplots(2, len(problem.data), figsize=(5*len(problem.data), 6), squeeze=False, layout="constrained", height_ratios=[3, 1])
    for i, (d, fit) in enumerate(zip(problem.data, prediction)):
        axes[0, i].errorbar(d.q, d.intensity, yerr=d.uncertainty, fmt=".", ms=3, label=d.label)
        axes[0, i].plot(d.q, fit, label="calculated")
        if np.all(d.intensity > 0) and np.all(fit > 0):
            axes[0, i].set_yscale("log")
        axes[0, i].set(ylabel="I(q) [cm$^{-1}$]")
        axes[0, i].legend()
        axes[1, i].axhline(0, color="0.5", lw=0.7)
        axes[1, i].plot(d.q, (fit-d.intensity)/d.uncertainty)
        axes[1, i].set(xlabel="q [A$^{-1}$]", ylabel="residual / dI")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    stamp = Path(path).parent/"provenance.json"
    if stamp.exists():
        provenance = json.loads(stamp.read_text(encoding="utf-8"))
        provenance["outputs"][Path(path).name] = sha256(path)
        stamp.write_text(json.dumps(provenance, indent=2)+"\n", encoding="utf-8", newline="\n")


def run_fit(config, output, *, plot=False):
    if Path(output).exists():
        raise FileExistsError(f"Run directory already exists: {output}")
    with threadpool_limits(limits=1):
        problem = Problem(config)
        result, prediction = fit_problem(problem)
        write_run(problem, result, prediction, output)
        if plot:
            plot_run(problem, prediction, Path(output)/"fit.png")
    return result


def check_convergence(config, fit_result, output, factor=1.5):
    if Path(output).exists():
        raise FileExistsError(f"Run directory already exists: {output}")
    if factor <= 1:
        raise ValueError("Refinement factor must exceed one")
    with threadpool_limits(limits=1):
        coarse = Problem(config)
        saved = json.loads(Path(fit_result).read_text(encoding="utf-8"))
        p, scale, bg = saved["parameters"], saved["scale"], saved["backgrounds"]
        a = [scale*y+b for y, b in zip(coarse.raw(p)[0], bg)]
        n = coarse.numerics
        # Increase cutoffs and nodes together; shell cutoff separately doubles.
        refined = replace(n, ns=int(np.ceil(n.ns*factor)), nz=int(np.ceil(n.nz*factor)),
                          s_max=n.s_max*factor, z_max=n.z_max*factor, max_shell=2*n.max_shell)
        fine = Problem(config, numerical_override=refined)
        b = [scale*y+v for y, v in zip(fine.raw(p)[0], bg)]
        metrics = [dict(dataset=i, max_change_in_dI=float(np.max(np.abs(y-x)/d.uncertainty)),
                        rms_change_in_dI=float(np.sqrt(np.mean(((y-x)/d.uncertainty)**2))))
                   for i, (x, y, d) in enumerate(zip(a, b, coarse.data))]
        result = dict(parameters=p, scale=scale, backgrounds=bg, metrics=metrics,
                      coarse_numerics=asdict(n), note="Fixed structural and linear parameters. No refit. A selected-case refinement is not uniform convergence over the parameter domain.")
        fine.inputs["fit_result"] = sha256(fit_result)
        write_run(fine, result, b, output, kind="convergence")
    return result
