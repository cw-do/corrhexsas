# Model and numerical conventions

The periodic reference SLD has hexagonal lattice spacing `a`. A transverse
Gaussian displacement field acts at each point of the composition texture.
Its components have variance `sigma_u**2` and normalized covariance

```text
H(s,z) = exp[-s²/(2 ell_b²) - z²/(2 ell_z²)]
rho(x_perp,z) = rho_0(x_perp - u(x_perp,z) + t)
```

The reference translation `t` is uniform over one lattice cell. It gives a
stationary ensemble. The displacement field is Gaussian. The resulting
SLD field is generally not Gaussian.

The connected density covariance is the sum over nonzero reciprocal modes

```text
C(r) = sum_G |d_G|² exp(i G.r_perp)
       * exp[-G² sigma_u² (1-H(s,z))]
```

For a radial motif, each complete reciprocal shell has
`G = 4 pi sqrt(h²-hk+k²)/(sqrt(3) a)`. Both signs are included once.
The radial coefficient is

```text
d_G = (2 pi/A_cell) integral [rho_rad(r)-rho_medium] J0(G r) r dr
```

Finite grains use normalized overlap of a circular transverse disk and
an axial interval. `grain_radius` is the disk radius, not a radius of
gyration. `grain_length` is the total interval length. Grain overlap
multiplies the covariance before the powder transform. It is not an
independently fitted peak-width function.

The response kernel is the Gaussian radial average of `sinc(q' R)`.
It uses the **one-standard-deviation radial sigma_q**, not FWHM and not
a three-dimensional vector resolution ellipsoid. It is the even-q
extension integrated over the real line. It reduces exactly to sinc
when sigma_q is zero. All radial layers interfere at amplitude level
before taking the shell power.

## Direct and series routes

The direct route evaluates the exponential covariance for every shell.
The series route uses the same function written as Poisson weights
`exp(-v) v**n/n!`, with `v=G² sigma_u²`. The transverse and axial factors
then separate. The current truncation is
`N=ceil(v_max+12 sqrt(v_max+1))`. `last_tail_bound` gives its
absolute intensity remainder bound before scale and background.
It does not bound quadrature, range or shell-cutoff errors, and is
not a relative bound near minima. This expansion and Gaussian phase
averaging are established mathematical tools.

The evaluator sums complete shells through `h²-hk+k² <= max_shell`.
The fixed separation grids allow the instrument kernel to be reused
during fitting. The kernel is cached when it fits within `cache_mb`;
otherwise q blocks are evaluated without storing the full response.
RAM also includes covariance and series work arrays. `cache_mb` limits
the response cache, not total process memory.

## Scope of the implementation

The default model has a sharp uniform disk or an explicitly supplied
radial SLD and a single grain size. No separate Gaussian motif-smoothing
parameter is used. The optional fixed `grain_radius_spread` is retained
only to replay historical calculations. It is a relative halfwidth of
a uniform **volume-weighted** radius distribution, not a number-weighted
distribution or an interface width.

The deformation acts on an Eulerian composition texture. It is not an
exact construction of rigid translated cylinders or of material tubes
with circular locally normal cross sections. It does not enforce mass
conservation under deformation, noncrossing, defect topology or an
elastic Hamiltonian. The covariance lengths are not measured elastic
moduli or wormlike-chain persistence lengths. Finite grain means are
assumed equal, so separate mean-contrast scattering between grains is
absent. This is a bulk powder model, not a GISAS/DWBA calculation.

Rigid-motif form-factor/structure-factor models remain appropriate for
their ensembles. Correlations alone do not remove the common form
factor of identical rigid objects. Here the mode powers `|d_G|²` weight
their own diffuse kernels. A continuous radial-transform zero at the
observed q does not force all those mode contributions to vanish.
This does not guarantee that every observed minimum is filled.

## Prior work

Full-curve hexagonal scattering and radial-profile modeling are established.
See [Frewein et al. (2019)](https://doi.org/10.1107/S1600576719002760).
Correlated displacement scattering in columnar phases also precedes
this package, including [Selinger and Bruinsma (1991)](https://doi.org/10.1103/PhysRevA.43.2910)
and [Welberry (2015)](https://doi.org/10.1515/zkri-2014-1732).
The global-analysis strategy has a lamellar counterpart in
[Pabst et al. (2000)](https://doi.org/10.1103/PhysRevE.62.4000).
This software provides the particular SLD/deformation/grain observation
model described in the associated manuscript. It does not claim that
Gaussian averaging, correlated hexagonal disorder or full-curve fitting
were introduced here.
