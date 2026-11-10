# Radial SLD profiles

All lengths in this page are in Angstrom. SLD numbers use units of
`1e-6 A^-2`. Profiles must become constant before radius `a/2`.
The zero reciprocal mode is omitted because the observation uses
connected density fluctuations.

## Disk and constant shells

```python
from corrhexsas import Disk, Shells

disk = Disk(radius=28, sld=1, medium=0)
layers = Shells(radii=[18, 28, 40], sld=[4, 1, 2], medium=3)
```

Each radius ends a shell. SLD arrays can also have shape
`(contrasts, shells)`. They share geometry and therefore preserve
the amplitude interference between radial layers.

For the CLI, a shell profile can reference declared parameter names
instead of fixed numbers.

```json
"profile": {
  "type": "shells",
  "radii": ["r_core", "r_outer"],
  "sld": [4.0, "rho_shell"],
  "medium": 3.0
}
```

Add `r_core`, `r_outer` and `rho_shell` to `parameters` with values and
optional bounds/vary flags. Bounds must keep radii ordered and below
half the **smallest allowed lattice spacing**. Invalid motifs raise
an error. The optimizer does not silently sort or clip radii.

## Tabulated profiles

```python
from corrhexsas import TabulatedProfile
profile = TabulatedProfile([0, 15, 25, 40], [4, 4, 1, 3], medium=3)
```

This is piecewise-linear interpolation and segmentwise Gaussian radial
quadrature. The table starts at `r=0`, has strictly increasing radii,
and ends at the matrix SLD. Multiple SLD columns provide contrasts.
For sharp jumps use `Shells` or zero-width `SmoothLayers` transitions
instead of duplicated radii. The CLI equivalent is

```json
"profile": {"type": "tabulated", "file": "radial_sld.csv", "medium": 3.0}
```

`radial_sld.csv` contains `r,sld` or `r,sld_contrast_1,sld_contrast_2,...`.
One header row is optional. Tabulated profiles in a CLI fit are fixed
functions of radius. Use the Python API for a custom deformable profile.

## Smooth layers

```python
from corrhexsas import SmoothLayers
profile = SmoothLayers(
    starts=[18, 30, 38], widths=[8, 0, 10], levels=[4, 1, 2, 3])
```

Every transition uses `S(t)=3t²-2t³` on `0<t<1`, with constant endpoints.
The first change spans 18–26 A, the second is a sharp jump at 30 A,
and the last spans 38–48 A. Widths are full support widths, not
Gaussian standard deviations or 10–90% widths. Amplitudes use the
integration-by-parts J1 formula, with 48 quadrature nodes per transition.

## Conserved two-solvent profile

`conserved_solvents` implements the composition-preserving profile in
the manuscript. A sharp polar–apolar boundary separates water/EO
and PO/oil regions. Compact transitions change solvent content within
each region. The outer matrix may contain polymer, and the central
plateau may contain water and polymer.

The default pure-component SLDs are water 6.33, EO 0.57 and PO 0.35.
The default polar area fraction is 0.34 and polar mean SLD is 3.34.
Oil SLD and mean apolar SLD are explicit arrays supplied for each contrast.
These defaults belong to the manuscript system. They are **not universal
material constants**. Supply the appropriate composition and SLDs for
another sample.

The bounded optimization coordinates are

| Coordinate | Meaning |
|---|---|
| `water_position` | EO-side water fraction as a fraction of the polar mean water fraction |
| `core_position` | Inner area-equivalent radius between its composition bound and the polar radius |
| `oil_position` | PO-side oil fraction between its conservation bounds |
| `outer_position` | Outer area-equivalent radius within its allowed interval |
| `inner_width_position` | Inner full support width as a fraction of its maximum allowed width |
| `outer_width_position` | Outer full support width as a fraction of its maximum allowed width |

All lie in `[0,1]`. These are coordinates, not lengths or independently
measured fractions. `decoded_profiles` in the output reports actual
radii, transition starts/ends, widths, solvent fractions and SLD levels.

If `outer_width` is supplied in the profile specification, it fixes the
full physical outer width. `outer_position` is then mapped between
the width-compatible lower/upper radii. `outer_width_position` should
be omitted. The inner transition remains inside the polar region.
The area-equivalent radius obeys

```text
C² = start² + start*width + 3*width²/10
```

It preserves component amounts when the transition width changes.
It is not the transition midpoint. Contrast-dependent matrix fractions
preserve each specified apolar mean, so strictly identical local
composition fractions across contrasts are not imposed.

`examples/chemical_template.json` demonstrates a six-contrast shared
fit. It requires your six data files and is a template, not a released
experimental dataset or the exact archived reference fit. In particular,
`normalize_to_initial` freezes weights from your own starting model.
The paper used an earlier separately frozen normalization model.
