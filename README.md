# corrhexsas

Small-angle scattering from a hexagonal SLD texture deformed by correlated
transverse displacements, with finite grains, powder averaging and radial
instrumental resolution. This is the reusable implementation associated with
the manuscript *Correlated displacements and finite grains in small-angle
scattering from hexagonal phases — what powder data can determine*.

The model connects the radial SLD profile, lateral/axial displacement
correlations and grain dimensions through one density covariance. It returns
Bragg and structured diffuse contributions within the same calculation.
The fitting interface reads ordinary intensity files and supports a shared
structure across multiple contrasts. No private workspace, measured dataset,
browser service or network connection is needed for analysis.

## Install

Python 3.10 or newer is required. We recommend a separate virtual environment.

```bash
git clone https://github.com/cw-do/corrhexsas.git
cd corrhexsas
python -m venv .venv
# Activate .venv using the command appropriate to your shell.
python -m pip install -e ".[plot]"
corrhexsas --version
```

On PowerShell, activate with `.\.venv\Scripts\Activate.ps1`. On macOS/Linux,
use `source .venv/bin/activate`. If activation is unavailable, invoke the
environment's Python directly as `python -m corrhexsas`.

For a pinned release without a working clone,

```bash
python -m pip install "corrhexsas[plot] @ git+https://github.com/cw-do/corrhexsas.git@v0.1.0"
```

The required numerical dependencies are NumPy, SciPy and threadpoolctl.
Matplotlib is optional and used only by `--plot`. This release is installed
from GitHub; it is not a PyPI publication.

## First fit in three commands

```bash
corrhexsas example --output demo
corrhexsas fit demo/disk.json --output runs/disk_fit --plot
corrhexsas check demo/disk.json --fit runs/disk_fit/result.json --output runs/disk_check
```

The example generates seeded synthetic observations with known parameters.
It fits lattice spacing, disk radius and displacement amplitude, while
keeping covariance lengths and grain dimensions fixed. It is a working
example, not a claim that all parameters are identifiable in an arbitrary
powder dataset. Every output directory must be new, so a second attempt
cannot overwrite the first.

## Fit your own I(q) file

Copy `demo/disk.json`, change `datasets[0].file` to your file, and set
physically appropriate starting values, bounds and fixed assumptions.
Data paths are relative to the configuration file. A CSV or whitespace
file may have one header line and `#` comments.

```text
q,I,dI,sigma_q
0.020,1.23,0.04,0.0017
0.021,1.30,0.04,0.0017
...
```

| Column | Meaning |
|---|---|
| `q` | Scattering-vector magnitude |
| `I` | Measured intensity, normally cm^-1 |
| `dI` | Positive one-standard-deviation intensity uncertainty |
| `sigma_q` | Optional nonnegative radial one-standard-deviation q resolution |

Use `"q_unit": "1/A"` or `"q_unit": "1/nm"`. The same conversion applies
to sigma_q. Internal lengths are always Angstrom, and numerical SLDs
always use `1e-6 A^-2`. Intensity is not converted by the loader.
If input intensity is in arbitrary units, the fitted scale absorbs that
convention and cannot be interpreted as an absolute volume fraction.

For three-column data, sigma_q defaults to zero unless `resolution` is
specified. For example, with q in `1/A`,

```json
"resolution": {"relative": 0.04, "absolute": 0.0015}
```

uses `sigma_q = hypot(0.04*q, 0.0015)`. A fourth data column takes
precedence. The absolute term has the dataset's q unit. This is a radial
Gaussian response, not a vector instrument resolution model. Convert
FWHM to a standard deviation before supplying it.

Two-column q/I input requires an explicit `relative_error` in its dataset
entry. This assumes `dI = relative_error*abs(I)` and is not a substitute
for measured uncertainties. Zero-intensity points therefore need an
explicit positive dI column. Malformed or nonfinite data raise an error;
rows are not silently dropped or reordered.

## Parameters and fixed assumptions

Every configuration parameter has a `value`. It varies only when
`"vary": true` is present. Varying coordinates require finite
`"bounds": [lower, upper]`. Add `"log": true` for logarithmic fitting
coordinates with positive bounds.

| Parameter | Definition |
|---|---|
| `a` | Nearest-neighbor hexagonal lattice spacing in A |
| `radius_over_a` | Disk radius divided by a, used only by the disk profile |
| `sigma_over_a` | Standard deviation of **one transverse displacement component**, divided by a |
| `ell_b_over_a` | Lateral displacement covariance length divided by a |
| `ell_z_over_a` | Axial displacement covariance length divided by a |
| `grain_radius_over_a` | Circular lateral grain radius divided by a |
| `grain_length_over_a` | Total axial grain length divided by a |

Ratios stay fixed when a changes unless those ratios are selected to vary.
The Python API instead accepts all dimensions directly through `Structure`.
For displacements constant along the column direction, set
`"axially_constant": true` and keep `ell_z_over_a` fixed. This retains
lateral disorder and finite grain length. It is different from fixing
`sigma_over_a` at zero. Correlation lengths cannot be fitted when sigma
is fixed at zero because they have no effect in that limit.

The common nonnegative scale and one constant background per dataset
are solved by constrained linear least squares at each structural trial.
Their default bounds are `[0, infinity]`. Override them through
`linear.scale_bounds` and `linear.background_bounds`, using JSON `null`
for an infinite endpoint. Report these bounds with your results.

## SLD profiles and contrast series

The package supports a uniform disk, sharp concentric shells, compact
smooth transitions, a tabulated radial SLD, and the composition-conserving
two-solvent profile used in the manuscript. See the runnable Python
examples and definitions in [docs/profiles.md](docs/profiles.md).

Add datasets to share the structural parameters and one intensity scale.
Supply `profiles`, with one specification per dataset, or a single
multi-contrast `profile`. Different q grids are supported. Each dataset
uses its own dI, sigma_q and fitted constant background.

`examples/chemical_template.json` is a six-contrast template that uses
the manuscript's component and mean SLDs. It requires your data files.
Those SLDs and the 0.34 polar area fraction are system-specific assumptions,
not defaults to transfer blindly to another material. The template
does not contain experimental data or reproduce the paper's archived
fit simply by being run on a different file.

## Objectives and outputs

The default objective is ordinary weighted chi-square. To give disjoint
q regions equal total importance after accounting for dI, use

```json
"fit": {
  "objective": "balanced",
  "region_edges": [0.040, 0.075, 0.098, 0.120],
  "normalize_to_initial": false,
  "starts": 3,
  "max_nfev": 250,
  "seed": 20261110
}
```

Edges are always internal `1/A` units. Every region must contain points
in every dataset. With `normalize_to_initial: false`, the objective is
the mean of the per-region mean standardized squared residuals across
datasets. With `true`, each region is additionally divided by the larger
of one and its starting-model mean squared standardized residual. These
normalizations are frozen **before optimization** and saved in the result.
They are not updated for each trial. The manuscript's archived reference
used a separately frozen earlier model; this optional starting-model
normalization does not silently recreate those weights.

A balanced objective is a declared scientific weighting choice, not a
likelihood or a confidence interval. Curve fitting is a local numerical
optimization. Multiple starts improve exploration but do not prove a
global optimum. Varying overall SLD contrast and an unconstrained intensity
scale together can also create a direct amplitude degeneracy.

Each fit writes

| File | Contents |
|---|---|
| `result.json` | Fitted/fixed parameter values, scale, backgrounds, objective, ordinary chi2, solver status, trials and decoded physical profile dimensions |
| `curve_00.csv`, etc. | q, observations, dI, sigma_q, calculated intensity and standardized residual |
| `config.json` | Snapshot of the requested analysis configuration |
| `provenance.json` | Input/output hashes, package/environment versions and source git commit when available |
| `fit.png` | Optional fit and residual plot |

Keep the original input files and their configuration together. Data paths
in the copied config remain relative to the original configuration's
directory; the run directory is a result record, not a complete data deposit.
Solver success does not establish parameter uniqueness. Nominal degrees
of freedom do not account for constraints, correlation or model discrepancy.

## Numerical reliability

The finite shell sum and quadrature need convergence checks for **your**
parameter and q ranges. Default settings are not a universal accuracy claim.
`Numerics` declares `ns`, `nz`, `s_max`, `z_max`, `max_shell`, the
`direct`/`series` route, q blocking and response-cache memory.

The physical support is `s <= 2*grain_radius` and `z <= grain_length`.
If the selected cutoffs omit support, evaluation stops unless
`allow_truncated_support: true` was explicitly selected. That option emits
a warning and must be justified by range refinements. It does not change
the grain dimensions. Negative calculated intensities produce an error
instead of being clipped. Very narrow peaks and small sigma_q require
more nodes. Sharp grain windows also create genuine finite-size fringes.

`corrhexsas check` keeps both structural and linear fit parameters fixed,
increases integration cutoffs and node counts by a factor of 1.5, doubles
the shell cutoff, and reports curve changes in supplied dI units. It does
not refit or automatically certify a tolerance. For rigorous work, also
increase nodes, range and shell cutoff separately and test representative
locations across the fitting bounds. Finer settings can require substantial
RAM and compute time. See [docs/model.md](docs/model.md) for the absolute
series-tail bound and cache behavior.

## Python API

```python
import numpy as np
from corrhexsas import Structure, Numerics, Evaluator, Disk

q = np.linspace(0.02, 0.16, 100)
sigma_q = np.hypot(0.04*q, 0.0015)
structure = Structure(a=140, sigma_u=5.6, ell_b=98, ell_z=210,
                      grain_radius=560, grain_length=1120)
evaluator = Evaluator(q, sigma_q, Numerics())
intensity = evaluator.evaluate(structure, Disk(radius=28))[0]
coherent = evaluator.evaluate(structure, Disk(radius=28), part="coherent")[0]
diffuse = evaluator.evaluate(structure, Disk(radius=28), part="diffuse")[0]
```

For your own radial model, supply an object with
`coefficients(G, a)` returning an array shaped `(contrasts, shells)`.
The evaluator and fitting workflow contain no hard-coded research paths.

## Validation and scientific limits

Run the numerical and fitting tests with

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m build
```

Tests include complete shell enumeration, independent resolution/radial
integrals, direct/series agreement, zero-displacement limits, data intake,
shared-contrast synthetic recovery, and regression against synthetic curves
frozen before this package existed. The additional private-data comparison
checks only stored calculated curves and does not redistribute observations.
See [validation/README.md](validation/README.md) for scope and provenance.

The Gaussian displacement covariance is a prescribed kinematic model,
not a derived elastic theory or a guaranteed material-tube construction.
Powder data can leave axial correlation and grain dimensions weakly
identified. Low fitting residuals alone do not make these measured
material constants. The code does not perform Bayesian inference or
provide a GISANS/GISAXS distorted-wave calculation. Prior scattering
models are credited in [docs/model.md](docs/model.md).

## Citation and license

Use the GitHub release tag and commit for reproducibility. `CITATION.cff`
provides software metadata. The associated manuscript is a working draft
by Changwoo Do, Chi-Huan Tung and Wei-Ren Chen; no journal publication or
DOI is claimed here. A persistent archival software/data DOI is a separate
step from this GitHub release. The code is distributed under the MIT license.
