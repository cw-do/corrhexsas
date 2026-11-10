"""Generate a reproducible synthetic disk dataset. No measured data included."""
from pathlib import Path
from dataclasses import asdict
import json
import numpy as np
from threadpoolctl import threadpool_limits
from .model import Structure, Numerics, Evaluator
from .profiles import Disk
from .data import sha256


def create_example(output):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    q = np.linspace(0.02, 0.16, 90)
    dq = np.hypot(0.04*q, 0.0015)
    p = Structure(140.0, 5.6, 98.0, 210.0, 560.0, 1120.0)
    n = Numerics(ns=512, nz=256, s_max=1400, z_max=1400, max_shell=64)
    with threadpool_limits(limits=1):
        raw = Evaluator(q, dq, n).evaluate(p, Disk(28.0))[0]
    truth = 0.8*raw+0.02
    dy = 0.02*truth+0.001
    y = truth+np.random.default_rng(20261110).normal(size=len(q))*dy
    np.savetxt(out/"synthetic.csv", np.c_[q, y, dy, dq], delimiter=",", header="q,I,dI,sigma_q", comments="")
    config = dict(datasets=[dict(file="synthetic.csv", label="synthetic disk")], q_unit="1/A",
                  profile=dict(type="disk", sld=1.0, medium=0.0),
                  parameters=dict(a=dict(value=138.0, vary=True, bounds=[132, 148]),
                                  radius_over_a=dict(value=0.19, vary=True, bounds=[0.16, 0.23]),
                                  sigma_over_a=dict(value=0.035, vary=True, bounds=[0.025, 0.055]),
                                  ell_b_over_a=dict(value=0.7), ell_z_over_a=dict(value=1.5),
                                  grain_radius_over_a=dict(value=4.0), grain_length_over_a=dict(value=8.0)),
                  numerics=asdict(n), linear=dict(scale_bounds=[0, None], background_bounds=[0, None]),
                  fit=dict(objective="chi2", max_nfev=120, starts=1, seed=20261110))
    (out/"disk.json").write_text(json.dumps(config, indent=2)+"\n", encoding="utf-8", newline="\n")
    (out/"truth.json").write_text(json.dumps(dict(kind="synthetic, not measured data", structure=asdict(p),
        radius=28.0, scale=0.8, background=0.02, seed=20261110,
        data_sha256=sha256(out/"synthetic.csv")), indent=2)+"\n", encoding="utf-8", newline="\n")
    return out/"disk.json"
