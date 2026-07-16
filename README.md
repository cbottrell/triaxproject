# TriaxProject

`triaxproject` is a compact numerical package for projecting a triaxial
ellipsoid along an arbitrary line of sight. It covers:

1. a uniform-density ellipsoid; and
2. any density stratified on similar ellipsoids, `rho = rho(m)`, within a
   finite outer shell `m_max`, where

   ```text
   m^2 = (x - x0)^T A (x - x0).
   ```

It creates raw line-of-sight surface-density maps, evaluates analytic and Abel
references, extracts contours from raster maps, and fits their centers, axis
ratios, ellipticities, and position angles.

![Worked projections and contour-fit residuals](example_outputs/constant_ellipticity.png)

## Main result demonstrated numerically

For a fixed ellipsoidal quadratic form `A`, a unit line-of-sight vector `n`,
and an orthonormal sky basis `P`, define

```text
alpha = n^T A n
d     = P^T A n
C     = P^T A P
B     = C - d d^T / alpha.
```

At sky position `R = (x,y)`, completing the square gives

```text
m^2 = alpha (s - s_center)^2 + xi^2,
xi^2 = R^T B R.
```

Thus every projected surface density for `rho = rho(m)` is a function only of
`xi`. If the eigenvalues of `B` are `lambda_1 <= lambda_2`, every contour has

```text
a_proj = xi / sqrt(lambda_1)
b_proj = xi / sqrt(lambda_2)
b/a    = sqrt(lambda_1 / lambda_2),
```

independent of contour level and of the radial profile. The package uses the
ellipticity convention `epsilon = 1 - b/a`.

The worked examples evaluate the full three-dimensional ellipsoidal radius at
every line-of-sight quadrature node. They do not construct the map from a
precomputed `Sigma(xi)`. The ellipses are then recovered from interpolated
raster contours rather than simply plotting the analytic contour equations.

## Install

From this directory:

```bash
python -m pip install -e '.[examples]'
```

NumPy and SciPy are runtime dependencies. Matplotlib and ContourPy are included
in the optional `examples` dependency group.

## Quick start

```python
import numpy as np
from triaxproject import (
    PolynomialDensity,
    TriaxialEllipsoid,
    surface_density_los,
)

ellipsoid = TriaxialEllipsoid(axes=(1.8, 1.1, 0.55))
projection = ellipsoid.project(line_of_sight=(0.43, 0.58, 0.69))

x = np.linspace(-1.8, 1.8, 401)
y = np.linspace(-1.0, 1.0, 301)
X, Y = np.meshgrid(x, y)

# rho(m) = rho0 (1-m^2)^2 for m <= 1
profile = PolynomialDensity(rho0=1.0, power=2.0, m_max=1.0)
Sigma = surface_density_los(
    projection, X, Y, profile, m_max=1.0, order=64
)

print(projection.axis_ratio)   # projected b/a
print(projection.ellipticity)  # 1-b/a
print(np.degrees(projection.position_angle))
```

Any vectorized callable can replace the supplied density profiles:

```python
def rho(m):
    return np.exp(-3.0 * np.asarray(m) ** 1.4)
```

`surface_density_los` truncates the integration at the requested `m_max`.
`surface_density_abel` supplies a high-accuracy one-dimensional reference for
the same finite support. If a profile has its own `m_max` attribute, as
`PolynomialDensity` does, the projector infers it when omitted and rejects a
conflicting explicit value. Profiles without an intrinsic cutoff default to
`m_max=1`; pass another value explicitly when needed.

## Run the worked examples

```bash
python examples/worked_examples.py
```

The script uses semiaxes `(1.8, 1.1, 0.55)` and the oblique line of sight
proportional to `(0.43, 0.58, 0.69)`. It produces:

- `example_outputs/constant_ellipticity.png`: maps, numerical contours,
  geometric ellipse overlays, and contour-fit residuals;
- `example_outputs/contour_measurements.csv`: one diagnostic row per contour;
- `example_outputs/example_summary.txt`: a compact human-readable report.

The second example uses the non-polynomial cored profile
`rho(m) = rho0 / (1 + (m/0.25)^2)`, truncated at `m=1`. For the included
501-by-501 calculation, the predicted axis ratio is `0.541139431881`. The
fitted `b/a` range across six contours is about `4.0e-7` for uniform density
and `4.2e-7` for the cored radial profile. The uniform raw map agrees with its
closed form to about `5e-14` of the central surface density; the non-polynomial
cored map agrees with its independent closed form to about `9e-11`.

You can trade resolution for speed:

```bash
python examples/worked_examples.py --grid-size 301 --quadrature-order 32
```

## Tests

The suite uses the standard library's `unittest`, so no test runner is needed.
Install `.[examples]` first because the raster-contour tests use ContourPy:

```bash
python -m unittest discover -s tests -v
```

The tests cover principal and oblique projections, rotated ellipsoids, the
spherical limit, raw LOS quadrature against three closed forms, arbitrary
`rho(m)` against the Abel reduction, constancy around predicted ellipses, and
axis-ratio recovery from raster contours.

## API overview

- `TriaxialEllipsoid`: intrinsic axes, orientation, center, and ellipsoidal
  radius.
- `ProjectionGeometry`: sky basis, projected matrix `B`, chord intersections,
  projected semiaxes, `b/a`, `1-b/a`, and position angle.
- `UniformDensity`, `PolynomialDensity`, `CoredPowerLawDensity`: ready-to-use
  profiles.
- `surface_density_los`: direct finite-chord Gauss-Legendre projection.
- `surface_density_abel`: nonsingular Abel-equivalent reference integral.
- `uniform_surface_density_analytic` and
  `polynomial_surface_density_analytic`: closed-form benchmarks.
- `cored_power_law_surface_density_analytic`: non-polynomial closed-form
  benchmark for outer slope two.
- `surface_density_map`: regular raw-LOS raster generation.
- `measure_contour`, `measure_contours`, and `fit_ellipse`: independent contour
  diagnostics.

## Scope and caveats

The fixed-ellipticity result assumes parallel projection and one fixed `A`:
all intrinsic isodensity shells must be concentric, coaxial, and similar. It is
not generally true when axis ratios or orientation vary with radius, when
density is a function of spherical rather than ellipsoidal radius inside an
ellipsoidal boundary, or when substructure is added.

For a non-monotonic `rho(m)`, one surface-density value can correspond to
multiple nested loops. Each loop is still a similar ellipse, but a level set
need not be a single contour. Position angle is numerically undefined for an
exactly circular projection.

The figure and worked-example scripts are source-tree assets; the minimal wheel
contains the importable library modules but not those auxiliary files.
