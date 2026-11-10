"""Finite-shell covariance quadrature. No morphology simulation is required.

Lengths are Angstrom, q and sigma_q are inverse Angstrom, and numerical
SLDs are in 1e-6 Angstrom^-2. Intensities are returned in cm^-1.
"""
from __future__ import annotations
from dataclasses import dataclass
import warnings
import numpy as np
from scipy.special import erfcx, j0, roots_legendre, gammaln, gammainc


class NumericalError(ValueError):
    """The requested numerical representation is inadequate or too large."""


@dataclass(frozen=True)
class Structure:
    a: float
    sigma_u: float
    ell_b: float
    ell_z: float
    grain_radius: float
    grain_length: float
    grain_radius_spread: float = 0.0

    def __post_init__(self):
        if not np.isfinite(self.a) or self.a <= 0:
            raise ValueError("a must be finite and positive")
        if not np.isfinite(self.sigma_u) or self.sigma_u < 0:
            raise ValueError("sigma_u must be finite and nonnegative")
        for name in ("ell_b", "ell_z", "grain_radius", "grain_length"):
            v = getattr(self, name)
            if v <= 0 or np.isnan(v) or (np.isinf(v) and name != "ell_z"):
                raise ValueError(f"{name} must be positive (only ell_z may be infinity)")
        if not 0 <= self.grain_radius_spread < 1:
            raise ValueError("grain_radius_spread must lie in [0,1)")


@dataclass(frozen=True)
class Numerics:
    ns: int = 512
    nz: int = 256
    s_max: float = 2000.0
    z_max: float = 2000.0
    max_shell: int = 64
    route: str = "series"
    cache_mb: float = 256.0
    q_block: int = 8
    allow_truncated_support: bool = False

    def __post_init__(self):
        for name in ("ns", "nz", "max_shell", "q_block"):
            v = getattr(self, name)
            if not isinstance(v, int) or isinstance(v, bool) or v < 1:
                raise ValueError(f"{name} must be a positive integer")
        if min(self.s_max, self.z_max) <= 0 or not np.all(np.isfinite([self.s_max, self.z_max, self.cache_mb])) or self.cache_mb < 0:
            raise ValueError("Finite positive separation cutoffs and nonnegative cache_mb are required")
        if self.route not in ("direct", "series"):
            raise ValueError("route must be direct or series")


def reciprocal_shells(cutoff):
    """Complete integer shell box, including both signs, without a second x6."""
    bound = int(np.ceil(np.sqrt(4 * cutoff / 3)))
    h, k = np.meshgrid(np.arange(-bound, bound + 1), np.arange(-bound, bound + 1))
    nu = h*h - h*k + k*k
    return np.unique(nu[(nu > 0) & (nu <= cutoff)], return_counts=True)


def resolution_kernel(q, sigma_q, r):
    """Gaussian radial average of sinc(q' r), with the even-q convention.

    This is a radial scalar response, not a vector resolution ellipsoid.
    Zero sigma and zero separation are handled analytically.
    """
    q, dq, r = np.asarray(q, float), np.asarray(sigma_q, float), np.asarray(r, float)
    result = np.sinc(q[:, None] * r[None, :] / np.pi)
    rows = np.flatnonzero(dq > 0)
    cols = np.flatnonzero(r > 0)
    if len(rows) and len(cols):
        x = dq[rows, None] * r[None, cols] / np.sqrt(2)
        y = q[rows, None] / (np.sqrt(2) * dq[rows, None])
        z = x + 1j*y
        result[np.ix_(rows, cols)] = np.sqrt(np.pi)/(2*x) * (
            np.exp(-y*y) - np.exp(-x*x)*np.real(np.exp(-2j*x*y)*erfcx(z)))
    return result


def disk_overlap(s, radius):
    t = np.clip(np.asarray(s)/(2*radius), 0, 1)
    return 2/np.pi * (np.arccos(t) - t*np.sqrt(1-t*t))


class Evaluator:
    """Cache the response for one q/sigma_q grid, reuse it during fitting.

    A profile supplies coefficients(G, a) with shape (contrasts, shells).
    See profiles.py for supported motifs or implement that small interface.
    """
    def __init__(self, q, sigma_q=None, numerics=None):
        self.q = np.asarray(q, float)
        self.dq = np.zeros_like(self.q) if sigma_q is None else np.broadcast_to(np.asarray(sigma_q, float), self.q.shape).copy()
        if self.q.ndim != 1 or not len(self.q) or not np.all(np.isfinite(self.q)) or np.any(self.q < 0):
            raise ValueError("q must be a finite nonnegative one-dimensional array")
        if np.any(self.dq < 0) or not np.all(np.isfinite(self.dq)):
            raise ValueError("sigma_q must be finite and nonnegative")
        self.numerics = n = numerics or Numerics()
        x, w = roots_legendre(n.ns)
        self.s, ws = (x+1)*n.s_max/2, w*n.s_max/2
        x, w = roots_legendre(n.nz)
        self.z, wz = (x+1)*n.z_max/2, w*n.z_max/2
        self.weights = 4*np.pi*(self.s*ws)[:, None]*wz[None, :]
        self.r = np.hypot(self.s[:, None], self.z[None, :]).ravel()
        self.shells, self.multiplicity = reciprocal_shells(n.max_shell)
        self.kernel = None
        size_mb = len(self.q)*self.r.size*8/1024**2
        if size_mb <= n.cache_mb:
            self.kernel = np.empty((len(self.q), self.r.size))
            for sl in self._blocks():
                self.kernel[sl] = resolution_kernel(self.q[sl], self.dq[sl], self.r)
        self.last_tail_bound = None
        self._warned_cutoff = False

    def _blocks(self):
        for start in range(0, len(self.q), self.numerics.q_block):
            yield slice(start, start+self.numerics.q_block)

    def _window(self, p):
        n = self.numerics
        if n.s_max < 2*p.grain_radius*(1+p.grain_radius_spread) or n.z_max < p.grain_length:
            if not n.allow_truncated_support:
                raise NumericalError("Separation cutoff is smaller than the physical grain support. Increase s_max/z_max or explicitly set allow_truncated_support after a range check.")
            if not self._warned_cutoff:
                warnings.warn("Using explicitly truncated grain support. Check larger integration ranges for this parameter domain.", RuntimeWarning)
                self._warned_cutoff = True
        spread = p.grain_radius_spread
        if spread == 0:
            wp = disk_overlap(self.s, p.grain_radius)
        else:
            # Historical optional uniform VOLUME-weighted radius average.
            x, w = roots_legendre(48)
            lo, hi = (1-spread)*p.grain_radius, (1+spread)*p.grain_radius
            lower = np.maximum(lo, self.s/2)
            span = np.maximum(hi-lower, 0)
            radii = lower[:, None] + span[:, None]*(x+1)/2
            wp = span/(2*(hi-lo))*np.sum(w*disk_overlap(self.s[:, None], radii), axis=1)
        return wp[:, None]*np.maximum(1-self.z[None, :]/p.grain_length, 0)

    def evaluate(self, structure, profile, *, scale=1.0, background=0.0, part="total"):
        """Return (contrasts, q) intensities. part may be total/coherent/diffuse.

        The tail bound is absolute in cm^-1 before scale/background. It does
        not bound quadrature or reciprocal-shell error.
        """
        if part not in ("total", "coherent", "diffuse"):
            raise ValueError("part must be total, coherent or diffuse")
        p, n = structure, self.numerics
        window = self._window(p)
        g = 4*np.pi/(np.sqrt(3)*p.a)*np.sqrt(self.shells)
        coefficients = np.atleast_2d(profile.coefficients(g, p.a))
        if coefficients.shape[1] != len(g) or not np.all(np.isfinite(coefficients)):
            raise ValueError("Profile coefficients must be finite with one column per shell")
        powers = np.abs(coefficients)**2*self.multiplicity
        bessel = j0(g[:, None]*self.s[None, :])
        v = (g*p.sigma_u)**2
        coherent = part == "coherent"
        if n.route == "direct":
            h = np.exp(-0.5*((self.s[:, None]/p.ell_b)**2+(self.z[None, :]/p.ell_z)**2))
            cov = np.zeros((len(powers), len(self.s), len(self.z)))
            for j in range(len(g)):
                phase = np.full_like(h, np.exp(-v[j])) if coherent else np.exp(-v[j]*(1-h))
                if part == "diffuse":
                    phase -= np.exp(-v[j])
                cov += powers[:, j, None, None]*bessel[j, :, None]*phase
            self.last_tail_bound = 0.0
        else:
            vmax = float(v.max())
            last = 0 if coherent or vmax == 0 else int(np.ceil(vmax+12*np.sqrt(vmax+1)))
            if last > 4096:
                raise NumericalError("Series has more than 4096 orders. Use route=direct or reduce the trial domain.")
            order = np.arange(last+1)
            prob = np.ones((len(g), 1)) if vmax == 0 else np.exp(-v[:, None] + order*np.log(v[:, None]) - gammaln(order+1))
            if part == "diffuse":
                prob[:, 0] = 0
            lateral = np.exp(-0.5*(self.s[:, None]/p.ell_b)**2*order)
            axial = np.exp(-0.5*(self.z[:, None]/p.ell_z)**2*order)
            cov = np.stack([((bessel.T@(power[:, None]*prob))*lateral)@axial.T for power in powers])
            volume = np.pi*p.grain_radius**2*p.grain_length*(1+p.grain_radius_spread**2/3)
            self.last_tail_bound = 0.0 if coherent or vmax == 0 else float(1e-4*powers.sum(axis=1).max()*gammainc(last+1, vmax)*volume)
        weighted = (cov*window*self.weights).reshape(len(powers), -1)
        if self.kernel is not None:
            intensity = 1e-4*(self.kernel@weighted.T).T
        else:
            intensity = np.empty((len(powers), len(self.q)))
            for sl in self._blocks():
                intensity[:, sl] = 1e-4*(resolution_kernel(self.q[sl], self.dq[sl], self.r)@weighted.T).T
        if not np.all(np.isfinite(intensity)):
            raise NumericalError("Nonfinite intensity. Check numerical settings and trial parameters.")
        if np.min(intensity) < -1e-9*max(float(np.max(np.abs(intensity))), 1e-100):
            raise NumericalError("Negative calculated intensity. Refine quadrature, integration range and shell cutoff; intensities are not clipped.")
        bg = np.asarray(background, float)
        if bg.ndim == 1:
            bg = bg[:, None]
        return scale*intensity+bg
