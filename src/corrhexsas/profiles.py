"""Radial SLD motifs. SLDs use 1e-6 A^-2 and radii use A."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.special import j0, j1, roots_legendre


def _radii(radii, a):
    r = np.asarray(radii, float)
    if r.ndim != 1 or not len(r) or not np.all(np.isfinite(r)) or np.any(r <= 0) or np.any(np.diff(r) <= 0) or r[-1] >= a/2:
        raise ValueError("Radii must be positive, increasing, and strictly below a/2")
    return r


@dataclass(frozen=True)
class Disk:
    radius: float
    sld: float = 1.0
    medium: float = 0.0

    def coefficients(self, g, a):
        r = _radii([self.radius], a)[0]
        return np.atleast_2d((self.sld-self.medium)*2*np.pi*r*j1(g*r)/(g*(np.sqrt(3)*a*a/2)))


@dataclass(frozen=True)
class Shells:
    """Each radius ends a constant-SLD shell. Final outside SLD is medium.

    sld may be (shells,) or (contrasts, shells); medium is scalar or per contrast.
    """
    radii: object
    sld: object
    medium: object = 0.0

    def coefficients(self, g, a):
        r = _radii(self.radii, a)
        rho = np.atleast_2d(np.asarray(self.sld, float))
        if rho.shape[1] != len(r):
            raise ValueError("One SLD per shell is required")
        medium = np.broadcast_to(self.medium, (len(rho),))
        jumps = -np.diff(np.c_[rho, medium], axis=1)
        amplitude = 2*np.pi*r[:, None]*j1(r[:, None]*g)/(g*(np.sqrt(3)*a*a/2))
        return jumps@amplitude


@dataclass(frozen=True)
class TabulatedProfile:
    """Piecewise-linear profile, constant outside its last radius.

    radius begins at zero. sld has shape (radii,) or (contrasts, radii).
    Duplicate radii are not allowed. Use Shells/SmoothLayers for exact jumps.
    """
    radius: object
    sld: object
    medium: object
    nodes: int = 32

    def coefficients(self, g, a):
        r = np.asarray(self.radius, float)
        rho = np.atleast_2d(np.asarray(self.sld, float))
        if r.ndim != 1 or len(r) < 2 or r[0] != 0 or np.any(np.diff(r) <= 0) or r[-1] >= a/2 or not np.all(np.isfinite(r)):
            raise ValueError("Tabulated radii must begin at zero, increase, and end below a/2")
        medium = np.broadcast_to(self.medium, (len(rho),))
        if rho.shape[1] != len(r) or not np.all(np.isfinite(rho)) or not np.allclose(rho[:, -1], medium, rtol=0, atol=1e-10):
            raise ValueError("Tabulated SLD must have one entry per radius and end at medium SLD")
        x, w = roots_legendre(self.nodes)
        out = np.zeros((len(rho), len(g)))
        for i in range(len(r)-1):
            rr = (r[i]+r[i+1])/2+(r[i+1]-r[i])/2*x
            level = rho[:, i, None]+(rho[:, i+1]-rho[:, i])[:, None]*(x+1)/2
            out += (r[i+1]-r[i])/2*((level-medium[:, None])*(w*rr))@j0(g[:, None]*rr).T
        return 2*np.pi/(np.sqrt(3)*a*a/2)*out


@dataclass(frozen=True)
class SmoothLayers:
    """Compact smoothstep transitions and exact sharp jumps.

    starts/widths locate each SLD change; levels has one more entry than
    starts. Width zero is a sharp jump. All transitions must end below a/2.
    """
    starts: object
    widths: object
    levels: object
    nodes: int = 48

    def coefficients(self, g, a):
        starts, widths = np.asarray(self.starts, float), np.asarray(self.widths, float)
        rho = np.atleast_2d(np.asarray(self.levels, float))
        if starts.ndim != 1 or widths.shape != starts.shape or rho.shape[1] != len(starts)+1:
            raise ValueError("Transition arrays and SLD levels have inconsistent shapes")
        if not np.all(np.isfinite(np.r_[starts, widths, rho.ravel()])) or np.any(starts < 0) or np.any(widths < 0) or np.any(starts[1:] < (starts+widths)[:-1]-1e-9) or np.max(starts+widths) >= a/2:
            raise ValueError("Smooth transition supports overlap or leave the reference disk")
        x, w = roots_legendre(self.nodes)
        t, dw = (x+1)/2, w/2*6*((x+1)/2)*(1-(x+1)/2)
        amplitude = []
        for start, width in zip(starts, widths):
            rr = start+width*t
            amplitude.append(np.sum(dw*rr*j1(g[:, None]*rr), axis=1))
        return (rho[:, :-1]-rho[:, 1:])@np.asarray(amplitude)*2*np.pi/(g*(np.sqrt(3)*a*a/2))

    def radial(self, r):
        r = np.asarray(r, float)
        levels = np.atleast_2d(self.levels)
        out = np.broadcast_to(levels[:, 0, None], (len(levels), len(r))).copy()
        for i, (start, width) in enumerate(zip(self.starts, self.widths)):
            t = np.clip((r-start)/width, 0, 1) if width > 0 else (r >= start).astype(float)
            out += (levels[:, i+1]-levels[:, i])[:, None]*(t*t*(3-2*t))
        return out


def maximum_width(c, lower, upper):
    start = (-0.5*lower+np.sqrt(0.25*lower**2+0.3*(c*c-lower**2)))/0.3
    disc = upper**2-1.2*(upper**2-c*c)
    end = (upper-np.sqrt(disc))/0.6 if disc >= 0 else np.inf
    return max(0.0, min(start, end))


def conserved_solvents(a, positions, *, oil_sld, apolar_mean, polar_fraction=0.34,
                       polar_mean=3.34, water_sld=6.33, eo_sld=0.57,
                       po_sld=0.35, outer_width=None):
    """Composition-conserving two-interface profile from the manuscript.

    positions = water_position, core_position, oil_position, outer_position,
    inner_width_position, outer_width_position, each in [0,1]. These bounded
    coordinates enforce fractions and nonoverlapping support. Returned info
    reports the physical radii, widths and solvent fractions, not just coordinates.
    outer_width optionally fixes the full physical outer transition width in A.
    """
    z = np.asarray(positions, float)
    oil, apolar = np.atleast_1d(oil_sld), np.atleast_1d(apolar_mean)
    if z.shape != (6,) or np.any(z < 0) or np.any(z > 1) or not np.all(np.isfinite(z)):
        raise ValueError("Six solvent-profile coordinates in [0,1] are required")
    if oil.shape != apolar.shape or not 0 < polar_fraction < 1 or water_sld == eo_sld or np.any(oil == po_sld):
        raise ValueError("Invalid component SLDs or polar fraction")
    mean_water = (polar_mean-eo_sld)/(water_sld-eo_sld)
    alpha = (apolar-po_sld)/(oil-po_sld)
    if not 0 < mean_water < 1 or np.any(alpha <= 0) or np.any(alpha >= 1):
        raise ValueError("Phase means must be physical solvent/polymer mixtures")
    cell = np.sqrt(3)*a*a/2
    rp, limit = np.sqrt(polar_fraction*cell/np.pi), 0.499*a
    fw = mean_water*z[0]
    cmin = rp*np.sqrt(max(1e-12, (mean_water-fw)/(1-fw)))
    ci = cmin+z[1]*(rp-cmin)
    if outer_width is None:
        co = rp+z[3]*(limit-rp)
    else:
        width = float(outer_width)
        if width < 0 or not np.isfinite(width):
            raise ValueError("outer_width must be finite and nonnegative")
        lower = np.sqrt(rp*rp+rp*width+0.3*width*width)
        upper = np.sqrt(limit*limit-limit*width+0.3*width*width)
        if upper < lower:
            raise ValueError("outer_width is too large for the apolar region")
        co = lower+z[3]*(upper-lower)
    aparea = (1-polar_fraction)*cell
    annulus = np.pi*(co*co-rp*rp)
    fmin = max(0.0, 1-(1-alpha.max())*aparea/max(annulus, 1e-15))
    fo = fmin+z[2]*(alpha.min()-fmin)
    wc = fw+(mean_water-fw)*rp*rp/(ci*ci)
    matrix = (alpha*aparea-fo*annulus)/(aparea-annulus)
    wi = z[4]*maximum_width(ci, 0, rp)
    si = -wi/2+np.sqrt(ci*ci-0.05*wi*wi)
    ei = si+wi
    max_outer = maximum_width(co, max(rp, ei), limit)
    wo = z[5]*max_outer if outer_width is None else float(outer_width)
    if wo < 0 or wo > max_outer+1e-8:
        raise ValueError("Physical outer_width does not fit the current transition support")
    so = -wo/2+np.sqrt(co*co-0.05*wo*wo)
    levels = np.array([[eo_sld+wc*(water_sld-eo_sld), eo_sld+fw*(water_sld-eo_sld),
                        po_sld+fo*(o-po_sld), po_sld+m*(o-po_sld)] for o, m in zip(oil, matrix)])
    if not np.all(np.isfinite(levels)) or min(wc, fw, fo, matrix.min()) < -1e-8 or max(wc, fw, fo, matrix.max()) > 1+1e-8:
        raise ValueError("Solvent fractions left [0,1]")
    profile = SmoothLayers([si, rp, so], [wi, 0.0, wo], levels)
    info = dict(polar_radius=rp, core_equivalent_radius=ci, outer_equivalent_radius=co,
                inner_start=si, inner_end=ei, inner_width=wi, outer_start=so,
                outer_end=so+wo, outer_width=wo, water_in_core=wc, water_in_EO=fw,
                oil_in_PO=fo, matrix_oil_fraction=matrix.tolist(), sld_levels=levels.tolist())
    return profile, info
