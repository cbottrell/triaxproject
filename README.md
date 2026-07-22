# TriaxProject

`triaxproject` is a compact numerical package for projecting a triaxial
ellipsoid along an arbitrary line of sight. It covers:

1. a uniform-density ellipsoid; and
2. any density stratified on similar ellipsoids, $\rho=\rho(m)$, within a
   finite outer shell $m_{\max}$.

The ellipsoidal radius is defined by

$$
m^2=(\boldsymbol{x}-\boldsymbol{x}_0)^{\mathsf T}
    A(\boldsymbol{x}-\boldsymbol{x}_0).
$$

It creates raw line-of-sight surface-density maps, evaluates analytic and Abel
references, extracts contours from raster maps, and fits their centers, axis
ratios, ellipticities, and position angles.

![Worked projections and contour-fit residuals](example_outputs/constant_ellipticity.png)

## Analytic derivation

The result rests on three assumptions:

1. the projection is parallel, so every line of sight has the same direction;
2. one fixed positive-definite matrix $A$ describes the shape and orientation
   at every radius; and
3. the density is a function of ellipsoidal radius alone, $\rho=\rho(m)$,
   with finite support $m\le M$.

Under these assumptions the projected surface density depends on the two sky
coordinates only through one quadratic form. That fact, rather than the
particular choice of $\rho(m)$, is what fixes the ellipticity of every
projected contour.

### 1. Intrinsic ellipsoidal geometry

Let $(a,b,c)$ be the semiaxes in the principal-axis frame. Let $Q$ rotate
vectors from that frame into world coordinates, and define the symmetric
shape matrix

$$
A=Q\mathrm{diag}(a^{-2},b^{-2},c^{-2})Q^{\mathsf T}.
$$

All eigenvalues of $A$ are positive, so $A$ is positive definite. For a point
$\boldsymbol r$ measured from the ellipsoid center $\boldsymbol r_0$, define
the dimensionless ellipsoidal radius by

$$
m^2=(\boldsymbol r-\boldsymbol r_0)^{\mathsf T}
    A(\boldsymbol r-\boldsymbol r_0).
$$

In the principal-axis frame this is simply

$$
m^2=\frac{x'^2}{a^2}+\frac{y'^2}{b^2}+\frac{z'^2}{c^2}.
$$

Consequently, $m=1$ is the ellipsoid with semiaxes $(a,b,c)$ and $m=M$ has
semiaxes $(Ma,Mb,Mc)$. If $\rho=\rho(m)$, changing the density level changes
only $m$: every intrinsic isodensity surface is therefore concentric, coaxial,
and geometrically similar to every other one.

### 2. Coordinates on the sky and along the line of sight

Choose a unit line-of-sight vector $\boldsymbol n$. Choose two unit sky-plane
vectors $\boldsymbol e_x$ and $\boldsymbol e_y$ perpendicular to
$\boldsymbol n$, and collect them as the columns of

$$
P=(\boldsymbol e_x,\boldsymbol e_y).
$$

The basis is orthonormal, so

$$
P^{\mathsf T}P=I_2,
\qquad
P^{\mathsf T}\boldsymbol n=\boldsymbol 0.
$$

Let $\boldsymbol R=(x,y)^{\mathsf T}$ be a position on the sky measured from
the projected center. Every point on the corresponding line of sight can be
written as

$$
\boldsymbol r-\boldsymbol r_0=P\boldsymbol R+\boldsymbol n s,
$$

where $s$ is physical distance along the ray. Substitute this expression into
the definition of $m$:

$$
m^2=(P\boldsymbol R+\boldsymbol n s)^{\mathsf T}
    A(P\boldsymbol R+\boldsymbol n s).
$$

Expanding all four terms gives

$$
m^2=
\boldsymbol R^{\mathsf T}P^{\mathsf T}AP\boldsymbol R
+2s\boldsymbol R^{\mathsf T}P^{\mathsf T}A\boldsymbol n
+s^2\boldsymbol n^{\mathsf T}A\boldsymbol n.
$$

Introduce three quantities that depend only on the ellipsoid and the viewing
direction:

- $\alpha=\boldsymbol n^{\mathsf T}A\boldsymbol n$, a positive scalar;
- $\boldsymbol d=P^{\mathsf T}A\boldsymbol n$, a two-vector; and
- $C=P^{\mathsf T}AP$, a symmetric $2\times2$ matrix.

The ellipsoidal radius along a fixed ray is then the quadratic

$$
m^2=\alpha s^2
    +2s\boldsymbol d^{\mathsf T}\boldsymbol R
    +\boldsymbol R^{\mathsf T}C\boldsymbol R.
$$

### 3. Complete the square

The terms that contain $s$ satisfy

$$
\alpha s^2+2s\boldsymbol d^{\mathsf T}\boldsymbol R
=
\alpha\left(s+\frac{\boldsymbol d^{\mathsf T}\boldsymbol R}{\alpha}\right)^2
-\frac{(\boldsymbol d^{\mathsf T}\boldsymbol R)^2}{\alpha}.
$$

Define the location of closest approach along the ray by

$$
s_c=-\frac{\boldsymbol d^{\mathsf T}\boldsymbol R}{\alpha},
$$

and define the projected shape matrix

$$
B=C-\frac{\boldsymbol d\boldsymbol d^{\mathsf T}}{\alpha}.
$$

Because
$(\boldsymbol d^{\mathsf T}\boldsymbol R)^2
=\boldsymbol R^{\mathsf T}\boldsymbol d\boldsymbol d^{\mathsf T}
\boldsymbol R$, the completed square is

$$
m^2=\alpha(s-s_c)^2+\boldsymbol R^{\mathsf T}B\boldsymbol R.
$$

It is useful to name the remaining sky-plane quadratic form:

$$
\xi^2=\boldsymbol R^{\mathsf T}B\boldsymbol R.
$$

For a fixed sky position, $m^2$ is minimized at $s=s_c$ and its minimum value
is $\xi^2$. Since $A$ is positive definite, this minimum is positive for every
nonzero $\boldsymbol R$; hence $B$ is also positive definite.

The intersections with the outer shell follow by setting $m=M$:

$$
\alpha(s-s_c)^2+\xi^2=M^2.
$$

A real chord exists only when $\xi\le M$. Its two endpoints are

$$
s_{\pm}=s_c\pm\sqrt{\frac{M^2-\xi^2}{\alpha}}.
$$

At $\xi=M$ the two endpoints coincide, so the projected boundary is

$$
\boldsymbol R^{\mathsf T}B\boldsymbol R=M^2.
$$

Because $B$ is a fixed positive-definite $2\times2$ matrix, this boundary is
an ellipse.

For numerical work, the code uses an equivalent expression that is better
conditioned. Let $U=(P,\boldsymbol n)$ be the orthogonal matrix formed from the
sky basis and the line-of-sight vector. In $U^{\mathsf T}AU$, the blocks are
$C$, $\boldsymbol d$, and $\alpha$. The block-inverse formula says that the
sky-plane block of $(U^{\mathsf T}AU)^{-1}$ is $B^{-1}$. On the other hand,

$$
(U^{\mathsf T}AU)^{-1}=U^{\mathsf T}A^{-1}U,
$$

whose sky-plane block is $P^{\mathsf T}A^{-1}P$. Therefore

$$
B=\left(P^{\mathsf T}A^{-1}P\right)^{-1}.
$$

This form avoids subtracting nearly equal matrices when the intrinsic
ellipsoid is very elongated.

### 4. Projection of a general radial profile

For a sky position with $\xi\le M$, the surface density is the density
integrated between the two chord endpoints:

$$
\Sigma(\boldsymbol R)
=\int_{s_-}^{s_+}
\rho\left(\sqrt{\alpha(s-s_c)^2+\xi^2}\right)ds.
$$

Now introduce a centered, dimensionless line-of-sight coordinate

$$
t=\sqrt{\alpha}(s-s_c).
$$

This substitution gives

$$
ds=\frac{dt}{\sqrt{\alpha}},
\qquad
m^2=\xi^2+t^2.
$$

The chord endpoints become

$$
t_-=-\sqrt{M^2-\xi^2},
\qquad
t_+=\sqrt{M^2-\xi^2}.
$$

The integrand depends on $t$ only through $t^2$, so it is even. The integral
from $t_-$ to $t_+$ is therefore twice the integral from zero to $t_+$.

> [!IMPORTANT]
> **General line-of-sight projection**
>
> $$
> \Sigma(\xi)=\frac{2}{\sqrt{\alpha}}
> \int_0^{\sqrt{M^2-\xi^2}}
> \rho\left(\sqrt{\xi^2+t^2}\right)dt
> $$
>
> This expression applies for $0\le\xi\le M$. For $\xi>M$, the ray does not
> intersect the truncated ellipsoid and $\Sigma=0$.

This formula already proves the central geometric claim: the sky coordinates
$x$ and $y$ enter only through
$\xi^2=\boldsymbol R^{\mathsf T}B\boldsymbol R$. The viewing direction affects
the normalization through $\alpha$ and the contour shape through $B$, but it
does not introduce any second sky coordinate into the radial profile.

For an Abel form, change variables from $t$ to $m$. From
$m^2=\xi^2+t^2$ on the positive half of the chord,

$$
t=\sqrt{m^2-\xi^2}.
$$

Differentiating gives

$$
2t dt=2m dm,
$$

and hence

$$
dt=\frac{m}{\sqrt{m^2-\xi^2}}dm.
$$

When $t=0$, $m=\xi$; when $t=\sqrt{M^2-\xi^2}$, $m=M$. Substitution therefore
gives the equivalent Abel integral.

> [!IMPORTANT]
> **Equivalent Abel form**
>
> $$
> \Sigma(\xi)=\frac{2}{\sqrt{\alpha}}
> \int_{\xi}^{M}
> \frac{\rho(m)m}{\sqrt{m^2-\xi^2}}dm
> $$

The two boxed results are mathematically identical. The first is convenient
for direct line-of-sight quadrature because its integrand is nonsingular. The
second makes contact with the standard Abel transform.

### 5. Why every contour has the same ellipticity

Write the projection formula schematically as

$$
\Sigma(\boldsymbol R)=F(\xi),
\qquad
\xi^2=\boldsymbol R^{\mathsf T}B\boldsymbol R.
$$

The scalar function $F$ contains all dependence on the radial density law.
The matrix $B$ contains all dependence on the projected geometry. Crucially,
$B$ does not depend on $\xi$ or on the chosen surface-density level.

Because $B$ is real, symmetric, and positive definite, it has orthonormal
eigenvectors and positive eigenvalues. Order them as

$$
0<\lambda_1\le\lambda_2.
$$

Rotate the sky coordinates into those eigenvectors and call the coordinates
$(u,v)$. In this basis,

$$
\xi^2=\lambda_1u^2+\lambda_2v^2.
$$

For a fixed value $\xi=\xi_0$, divide by $\xi_0^2$:

$$
\frac{u^2}{\xi_0^2/\lambda_1}
+\frac{v^2}{\xi_0^2/\lambda_2}
=1.
$$

This is the standard equation of an ellipse. Since
$\lambda_1\le\lambda_2$, the $u$ direction is the major-axis direction. The
semiaxis lengths are

$$
a_{\mathrm{proj}}=\frac{\xi_0}{\sqrt{\lambda_1}},
\qquad
b_{\mathrm{proj}}=\frac{\xi_0}{\sqrt{\lambda_2}}.
$$

Taking their ratio cancels $\xi_0$.

> [!IMPORTANT]
> **Fixed projected shape**
>
> $$
> q_{\mathrm{proj}}
> =\frac{b_{\mathrm{proj}}}{a_{\mathrm{proj}}}
> =\sqrt{\frac{\lambda_1}{\lambda_2}},
> \qquad
> \epsilon=1-q_{\mathrm{proj}}
> $$
>
> The center, eigenvector directions, axis ratio, and ellipticity are
> independent of contour level and of the function $\rho(m)$.

If $F$ is one-to-one, a specified surface density $\Sigma_0$ selects one value
$\xi_0$ and hence one ellipse. If $F$ is non-monotonic, the equation
$F(\xi)=\Sigma_0$ may have several roots. Those roots produce several nested
ellipses, but every ellipse uses the same matrix $B$, so all have the same
axis ratio and position angle.

Changing the viewing direction changes $\alpha$, $P$, and $B$. Thus the fixed
ellipticity statement applies within one projection; the numerical value of
$q_{\mathrm{proj}}$ generally differs between viewing directions.

### 6. Closed-form benchmark profiles

The package includes three finite-support profiles for which the general
projection integral can be evaluated analytically. These are useful both as
examples and as numerical tests.

#### 6.1 Uniform density

Let $\rho(m)=\rho_0$ for $0\le m\le M$ and zero outside. Define the positive
half-chord length in the transformed coordinate by

$$
L=\sqrt{M^2-\xi^2}.
$$

The density is constant along the chord, so the general projection formula
reduces immediately to

$$
\Sigma(\xi)
=\frac{2}{\sqrt{\alpha}}\int_0^L\rho_0 dt
=\frac{2\rho_0L}{\sqrt{\alpha}}.
$$

> [!IMPORTANT]
> **Uniform-density projection**
>
> $$
> \Sigma(\xi)=\frac{2\rho_0}{\sqrt{\alpha}}
> \sqrt{M^2-\xi^2}
> $$

The square root is simply the half-chord length. For a chosen surface density
$\Sigma_0$, solving the result for $\xi$ gives

$$
\xi^2=M^2-\frac{\alpha\Sigma_0^2}{4\rho_0^2}.
$$

Thus each uniform-density contour is explicitly a constant-$\xi$ ellipse.

#### 6.2 Polynomial density

Consider

$$
\rho(m)=\rho_0
\left[1-\left(\frac{m}{M}\right)^2\right]^p,
\qquad
0\le m\le M,
\qquad
p\ge0.
$$

Along a ray, $m^2=\xi^2+t^2$. With
$L^2=M^2-\xi^2$, the bracket becomes

$$
1-\frac{m^2}{M^2}
=\frac{M^2-\xi^2-t^2}{M^2}
=\frac{L^2-t^2}{M^2}.
$$

The projection integral is therefore

$$
\Sigma(\xi)
=\frac{2\rho_0}{\sqrt{\alpha}M^{2p}}
\int_0^L(L^2-t^2)^pdt.
$$

Set $u=t^2/L^2$. Then

$$
t=L\sqrt{u},
\qquad
dt=\frac{L}{2\sqrt{u}}du,
$$

and the limits remain zero and one. Substitution gives

$$
\Sigma(\xi)
=\frac{\rho_0L^{2p+1}}{\sqrt{\alpha}M^{2p}}
\int_0^1u^{-1/2}(1-u)^pdu.
$$

The remaining integral is the Euler beta function

$$
\mathrm{B}\left(\frac{1}{2},p+1\right)
=\int_0^1u^{-1/2}(1-u)^pdu.
$$

Finally,

$$
\frac{L^{2p+1}}{M^{2p}}
=M\left(1-\frac{\xi^2}{M^2}\right)^{p+1/2}.
$$

> [!IMPORTANT]
> **Polynomial-profile projection**
>
> $$
> \Sigma(\xi)=\frac{\rho_0M}{\sqrt{\alpha}}
> \mathrm{B}\left(\frac{1}{2},p+1\right)
> \left(1-\frac{\xi^2}{M^2}\right)^{p+1/2}
> $$

Here $\mathrm{B}$ is the Euler beta function; it is unrelated to the projected
shape matrix $B$.

#### 6.3 Truncated cored profile

For the cored profile with outer slope two,

$$
\rho(m)=\frac{\rho_0}{1+(m/r_c)^2}
=\frac{\rho_0r_c^2}{r_c^2+m^2},
\qquad
0\le m\le M,
\qquad
r_c>0.
$$

Using $m^2=\xi^2+t^2$ and
$L=\sqrt{M^2-\xi^2}$ gives

$$
\Sigma(\xi)
=\frac{2\rho_0r_c^2}{\sqrt{\alpha}}
\int_0^L\frac{dt}{r_c^2+\xi^2+t^2}.
$$

Let $K=\sqrt{r_c^2+\xi^2}$. The required antiderivative is

$$
\int\frac{dt}{K^2+t^2}
=\frac{1}{K}\tan^{-1}\left(\frac{t}{K}\right).
$$

Evaluating it between zero and $L$ yields the third closed form.

> [!IMPORTANT]
> **Truncated-cored projection**
>
> $$
> \Sigma(\xi)=
> \frac{2\rho_0r_c^2}{\sqrt{\alpha}\sqrt{r_c^2+\xi^2}}
> \tan^{-1}\left(
> \frac{\sqrt{M^2-\xi^2}}{\sqrt{r_c^2+\xi^2}}
> \right)
> $$

#### 6.4 Principal-axis consistency check

As a transparent special case, take an unrotated ellipsoid

$$
\frac{x^2}{a^2}+\frac{y^2}{b^2}+\frac{z^2}{c^2}\le1
$$

and view it along the $z$ axis. Then

$$
\alpha=\frac{1}{c^2},
\qquad
\xi^2=\frac{x^2}{a^2}+\frac{y^2}{b^2}.
$$

Substituting these quantities into the uniform-density result gives

$$
\Sigma(x,y)=2\rho_0c
\sqrt{1-\frac{x^2}{a^2}-\frac{y^2}{b^2}}.
$$

Fixing $\Sigma$ fixes
$x^2/a^2+y^2/b^2$, so every contour is a scaled copy of the projected
boundary ellipse. Its axis ratio is $b/a$ when $a\ge b$.

### 7. Numerical verification

The worked examples evaluate the full three-dimensional ellipsoidal radius at
every line-of-sight quadrature node. They do not construct the map from a
precomputed $\Sigma(\xi)$. The ellipses are then recovered from interpolated
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

The script uses semiaxes $(1.8,1.1,0.55)$ and the oblique line of sight
proportional to $(0.43,0.58,0.69)$. It produces:

- `example_outputs/constant_ellipticity.png`: maps, numerical contours,
  geometric ellipse overlays, and contour-fit residuals;
- `example_outputs/contour_measurements.csv`: one diagnostic row per contour;
- `example_outputs/example_summary.txt`: a compact human-readable report.

The second example uses the non-polynomial cored profile
$\rho(m)=\rho_0/[1+(m/0.25)^2]$, truncated at $m=1$. For the included
$501\times501$ calculation, the predicted axis ratio is $0.541139431881$. The
fitted $b/a$ range across six contours is about $4.0\times10^{-7}$ for uniform
density and $4.2\times10^{-7}$ for the cored radial profile. The uniform raw
map agrees with its closed form to about $5\times10^{-14}$ of the central
surface density; the non-polynomial cored map agrees with its independent
closed form to about $9\times10^{-11}$.

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
$\rho(m)$ against the Abel reduction, constancy around predicted ellipses, and
axis-ratio recovery from raster contours.

## API overview

- `TriaxialEllipsoid`: intrinsic axes, orientation, center, and ellipsoidal
  radius.
- `ProjectionGeometry`: sky basis, projected matrix $B$, chord intersections,
  projected semiaxes, $b/a$, $1-b/a$, and position angle.
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

The fixed-ellipticity result assumes parallel projection and one fixed $A$:
all intrinsic isodensity shells must be concentric, coaxial, and similar. It is
not generally true when axis ratios or orientation vary with radius, when
density is a function of spherical rather than ellipsoidal radius inside an
ellipsoidal boundary, or when substructure is added.

For a non-monotonic $\rho(m)$, one surface-density value can correspond to
multiple nested loops. Each loop is still a similar ellipse, but a level set
need not be a single contour. Position angle is numerically undefined for an
exactly circular projection.

The figure and worked-example scripts are source-tree assets; the minimal wheel
contains the importable library modules but not those auxiliary files.
