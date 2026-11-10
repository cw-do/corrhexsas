import json
from dataclasses import asdict
import numpy as np
import pytest
from threadpoolctl import threadpool_limits
from corrhexsas import Evaluator, Numerics, Structure, Disk, load_data
from corrhexsas.workflow import Problem, fit_problem


def test_data_units_header_and_missing_errors(tmp_path):
    f = tmp_path/"data.csv"
    f.write_text("q,I,dI,sigma_q\n0.2,2,0.1,0.01\n0.3,3,0.2,0.02\n", encoding="utf-8")
    data = load_data(f, q_unit="1/nm")
    np.testing.assert_allclose(data.q, [.02, .03])
    np.testing.assert_allclose(data.sigma_q, [.001, .002])
    f.write_text("0.02 2\n0.03 3\n", encoding="utf-8")
    with pytest.raises(ValueError, match="relative_error"):
        load_data(f)
    np.testing.assert_allclose(load_data(f, relative_error=.05).uncertainty, [.1, .15])
    f.write_text("q I dI\n0.02 2 .1\nBROKEN\n", encoding="utf-8")
    with pytest.raises(ValueError, match="line 3"):
        load_data(f)


def test_shared_fit_recovers_synthetic_geometry_and_two_backgrounds(tmp_path):
    q = np.linspace(.025, .15, 45)
    dq = np.hypot(.04*q, .0015)
    n = Numerics(ns=384, nz=192, s_max=1100, z_max=1100, max_shell=28)
    p = Structure(140, 5.6, 98, 210, 420, 840)
    with threadpool_limits(limits=1):
        raw = Evaluator(q, dq, n).evaluate(p, Disk(28))[0]
    datasets = []
    for i, (contrast, background) in enumerate([(1, .02), (.7, .03)]):
        f = tmp_path/f"data_{i}.csv"
        np.savetxt(f, np.c_[q, .8*contrast**2*raw+background, np.full_like(q, .01), dq], delimiter=",")
        datasets.append(dict(file=f.name))
    cfg = dict(datasets=datasets, profiles=[dict(type="disk", sld=1), dict(type="disk", sld=.7)],
               parameters=dict(a=dict(value=138, vary=True, bounds=[134, 146]),
                               radius_over_a=dict(value=.19, vary=True, bounds=[.17, .23]),
                               sigma_over_a=dict(value=.04), ell_b_over_a=dict(value=.7),
                               ell_z_over_a=dict(value=1.5), grain_radius_over_a=dict(value=3),
                               grain_length_over_a=dict(value=6)), numerics=asdict(n))
    path = tmp_path/"fit.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    with threadpool_limits(limits=1):
        result, curves = fit_problem(Problem(path))
    assert result["success"]
    assert result["chi2"] < 1e-6
    np.testing.assert_allclose(result["parameters"]["a"], 140, atol=1e-4)
    np.testing.assert_allclose(result["parameters"]["radius_over_a"], .2, atol=1e-5)
    np.testing.assert_allclose(result["scale"], .8, atol=1e-5)
    np.testing.assert_allclose(result["backgrounds"], [.02, .03], atol=1e-5)
    cfg["fit"] = dict(objective="balanced", region_edges=[.06, .10], normalize_to_initial=True)
    path.write_text(json.dumps(cfg), encoding="utf-8")
    with threadpool_limits(limits=1):
        problem = Problem(path)
        before = [s.copy() for s in problem.scales]
        problem.residual(problem.encode())
    for a, b in zip(before, problem.scales):
        np.testing.assert_array_equal(a, b)
    assert len(problem.regions) == 6
