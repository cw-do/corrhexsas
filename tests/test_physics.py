from dataclasses import replace
import numpy as np
import pytest
from scipy.special import roots_hermitenorm, j0, roots_legendre
from threadpoolctl import threadpool_limits
from corrhexsas import Evaluator, Numerics, Structure, Disk, Shells, SmoothLayers, conserved_solvents
from corrhexsas.model import reciprocal_shells, resolution_kernel, NumericalError


def test_complete_shell_enumeration():
    shell, multiplicity = reciprocal_shells(7)
    np.testing.assert_array_equal(shell, [1, 3, 4, 7])
    np.testing.assert_array_equal(multiplicity, [6, 6, 6, 12])
    shell, multiplicity = reciprocal_shells(256)
    assert multiplicity[-1] > 0
    # Independent wider box catches the historical fixed +/-10 enumeration issue.
    expected = {}
    for h in range(-40, 41):
        for k in range(-40, 41):
            v = h*h-h*k+k*k
            if 0 < v <= 256:
                expected[v] = expected.get(v, 0)+1
    assert dict(zip(shell, multiplicity)) == expected


def test_resolution_kernel_against_independent_gaussian_integral():
    q, dq = np.array([0, 0.05, 0.14]), np.array([0.002, 0.003, 0.006])
    r = np.array([0, 0.01, 1, 50, 200, 700])
    x, w = roots_hermitenorm(256)
    expected = np.array([w@np.sinc((qi+di*x[:, None])*r[None, :]/np.pi)/np.sqrt(2*np.pi) for qi, di in zip(q, dq)])
    np.testing.assert_allclose(resolution_kernel(q, dq, r), expected, atol=2e-11, rtol=1e-9)
    np.testing.assert_allclose(resolution_kernel(q, np.zeros(3), r), np.sinc(q[:, None]*r/np.pi), atol=0)


def test_smooth_profile_coefficients_against_independent_radial_integral():
    profile = SmoothLayers([18, 30, 38], [8, 0, 10], [[4, 1, 2, 3], [3, 1, 1.5, 2]])
    a, g = 140.0, np.array([0.05, 0.11, 0.24, 0.45])
    x, w = roots_legendre(128)
    numerical = np.zeros((2, len(g)))
    for lo, hi in zip([0, 18, 26, 30, 38], [18, 26, 30, 38, 48]):
        r = (lo+hi)/2+(hi-lo)/2*x
        rho = profile.radial(r)-np.array([3, 2])[:, None]
        numerical += (hi-lo)/2*(rho*w*r)@j0(g[:, None]*r).T
    numerical *= 2*np.pi/(np.sqrt(3)*a*a/2)
    np.testing.assert_allclose(profile.coefficients(g, a), numerical, atol=1e-12)
    np.testing.assert_allclose(SmoothLayers([28], [0], [1, 0]).coefficients(g, a), Disk(28).coefficients(g, a), atol=1e-14)


def test_solvent_balance_and_bounded_fixed_outer_width():
    oil, apolar = np.array([5.84, 5.02]), np.array([3.61, 3.14])
    profile, info = conserved_solvents(140, [.5, .6, .5, .5, .8, .4], oil_sld=oil, apolar_mean=apolar, outer_width=10)
    assert info["outer_width"] == 10
    rp = info["polar_radius"]
    ci, co = info["core_equivalent_radius"], info["outer_equivalent_radius"]
    fw, wc = info["water_in_EO"], info["water_in_core"]
    np.testing.assert_allclose(fw+(wc-fw)*ci**2/rp**2, (3.34-.57)/(6.33-.57))
    area = np.sqrt(3)*140**2/2
    aa, annulus = .66*area, np.pi*(co*co-rp*rp)
    expected = (apolar-.35)/(oil-.35)
    np.testing.assert_allclose((info["oil_in_PO"]*annulus+np.array(info["matrix_oil_fraction"])*(aa-annulus))/aa, expected)


def test_routes_parts_zero_displacement_and_streaming():
    q = np.linspace(.025, .16, 25)
    dq = np.hypot(.04*q, .0015)
    p = Structure(140, 5.6, 98, 210, 420, 840)
    n = Numerics(ns=384, nz=192, s_max=1000, z_max=1000, max_shell=28)
    with threadpool_limits(limits=1):
        e = Evaluator(q, dq, n)
        direct = Evaluator(q, dq, replace(n, route="direct", cache_mb=0))
        motif = Shells([20, 28], [2, 1])
        total = e.evaluate(p, motif)
        np.testing.assert_allclose(total, direct.evaluate(p, motif), rtol=1e-9, atol=1e-10)
        np.testing.assert_allclose(total, e.evaluate(p, motif, part="coherent")+e.evaluate(p, motif, part="diffuse"), rtol=1e-9, atol=1e-10)
        assert e.last_tail_bound < 1e-10
        p0 = replace(p, sigma_u=0)
        np.testing.assert_array_equal(e.evaluate(p0, motif, part="diffuse"), np.zeros_like(total))
        np.testing.assert_allclose(e.evaluate(p0, motif), e.evaluate(replace(p0, ell_b=1, ell_z=np.inf), motif), rtol=1e-14)
        with pytest.raises(NumericalError, match="support"):
            e.evaluate(replace(p, grain_length=2000), motif)
