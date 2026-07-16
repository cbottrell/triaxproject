"""Analytic and numerical line-of-sight surface-density projections."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import quad
from scipy.special import beta as beta_function

from .geometry import ProjectionGeometry

DensityCallable = Callable[[ArrayLike], ArrayLike]


def _finite_positive(value: float, name: str) -> float:
    converted = float(value)
    if not np.isfinite(converted) or converted <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive")
    return converted


def _broadcast_sky(
    x: ArrayLike, y: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    return np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))


def _resolve_m_max(density: DensityCallable, m_max: float | None) -> float:
    profile_limit = getattr(density, "m_max", None)
    if m_max is None:
        return _finite_positive(1.0 if profile_limit is None else profile_limit, "m_max")
    resolved = _finite_positive(m_max, "m_max")
    if profile_limit is not None and not np.isclose(
        resolved, float(profile_limit), rtol=1.0e-12, atol=0.0
    ):
        raise ValueError(
            "m_max disagrees with density.m_max; use matching support radii"
        )
    return resolved


def uniform_surface_density_analytic(
    projection: ProjectionGeometry,
    x: ArrayLike,
    y: ArrayLike,
    rho0: float = 1.0,
    m_max: float = 1.0,
) -> NDArray[np.float64]:
    r"""Closed-form projection of a homogeneous ellipsoid.

    ``Sigma = 2 rho0 / sqrt(alpha) * sqrt(m_max^2 - xi^2)`` inside
    the projected boundary and zero outside.
    """
    rho0 = _finite_positive(rho0, "rho0")
    m_max = _finite_positive(m_max, "m_max")
    xi_squared = projection.projected_radius_squared(x, y)
    remaining = np.maximum(m_max * m_max - xi_squared, 0.0)
    result = 2.0 * rho0 * np.sqrt(remaining / projection.alpha)
    return np.where(xi_squared <= m_max * m_max, result, 0.0)


def polynomial_surface_density_analytic(
    projection: ProjectionGeometry,
    x: ArrayLike,
    y: ArrayLike,
    rho0: float = 1.0,
    power: float = 2.0,
    m_max: float = 1.0,
) -> NDArray[np.float64]:
    r"""Closed form for ``rho=rho0*(1-(m/m_max)^2)^power``.

    The profile is understood to vanish for ``m > m_max``.
    """
    rho0 = _finite_positive(rho0, "rho0")
    m_max = _finite_positive(m_max, "m_max")
    power = float(power)
    if not np.isfinite(power) or power < 0.0:
        raise ValueError("power must be finite and nonnegative")

    xi_squared = projection.projected_radius_squared(x, y)
    base = np.maximum(1.0 - xi_squared / (m_max * m_max), 0.0)
    normalization = (
        rho0
        * m_max
        * beta_function(0.5, power + 1.0)
        / np.sqrt(projection.alpha)
    )
    result = normalization * np.power(base, power + 0.5)
    return np.where(xi_squared <= m_max * m_max, result, 0.0)


def cored_power_law_surface_density_analytic(
    projection: ProjectionGeometry,
    x: ArrayLike,
    y: ArrayLike,
    rho0: float = 1.0,
    core_radius: float = 0.25,
    m_max: float = 1.0,
) -> NDArray[np.float64]:
    r"""Closed form for a truncated cored profile with outer slope two.

    This projects ``rho(m) = rho0 / (1 + (m/core_radius)^2)`` for
    ``m <= m_max``. It is a useful non-polynomial reference for the raw LOS
    quadrature.
    """
    rho0 = _finite_positive(rho0, "rho0")
    core_radius = _finite_positive(core_radius, "core_radius")
    m_max = _finite_positive(m_max, "m_max")

    xi_squared = projection.projected_radius_squared(x, y)
    inside = xi_squared <= m_max * m_max
    line_limit = np.sqrt(np.maximum(m_max * m_max - xi_squared, 0.0))
    denominator = np.sqrt(core_radius * core_radius + xi_squared)
    result = (
        2.0
        * rho0
        * core_radius
        * core_radius
        / (np.sqrt(projection.alpha) * denominator)
        * np.arctan(line_limit / denominator)
    )
    return np.where(inside, result, 0.0)


def surface_density_los(
    projection: ProjectionGeometry,
    x: ArrayLike,
    y: ArrayLike,
    density: DensityCallable,
    m_max: float | None = None,
    order: int = 64,
) -> NDArray[np.float64]:
    r"""Numerically integrate density along each raw 3D ellipsoid chord.

    Gauss--Legendre nodes span chord roots computed from the raw intrinsic
    quadratic, with a stable projected-radius fallback only when cancellation
    is detected. At every node, ``m^2`` is evaluated as a positive sum in the
    ellipsoid's intrinsic coordinates. The map is not assembled from a
    precomputed one-dimensional ``Sigma(xi)`` profile.
    """
    if not callable(density):
        raise TypeError("density must be callable")
    m_max = _resolve_m_max(density, m_max)
    if not isinstance(order, (int, np.integer)) or order < 2:
        raise ValueError("order must be an integer of at least 2")

    x_array, y_array = _broadcast_sky(x, y)
    intrinsic_transform = (
        np.diag(1.0 / projection.ellipsoid.axes)
        @ projection.ellipsoid.orientation.T
    )
    intrinsic_sky = intrinsic_transform @ projection.sky_basis
    intrinsic_los = intrinsic_transform @ projection.line_of_sight
    sky_intrinsic_0 = (
        intrinsic_sky[0, 0] * x_array + intrinsic_sky[0, 1] * y_array
    )
    sky_intrinsic_1 = (
        intrinsic_sky[1, 0] * x_array + intrinsic_sky[1, 1] * y_array
    )
    sky_intrinsic_2 = (
        intrinsic_sky[2, 0] * x_array + intrinsic_sky[2, 1] * y_array
    )
    beta = (
        sky_intrinsic_0 * intrinsic_los[0]
        + sky_intrinsic_1 * intrinsic_los[1]
        + sky_intrinsic_2 * intrinsic_los[2]
    )
    gamma = (
        sky_intrinsic_0 * sky_intrinsic_0
        + sky_intrinsic_1 * sky_intrinsic_1
        + sky_intrinsic_2 * sky_intrinsic_2
    )

    first_term = beta * beta
    second_term = projection.alpha * (m_max * m_max - gamma)
    raw_discriminant = first_term + second_term
    term_scale = np.abs(first_term) + np.abs(second_term)
    cancellation_limit = 4096.0 * np.finfo(float).eps * np.maximum(
        term_scale, np.finfo(float).tiny
    )
    needs_stable_fallback = (
        ~np.isfinite(raw_discriminant)
        | (np.abs(raw_discriminant) <= cancellation_limit)
    )
    if np.any(needs_stable_fallback):
        xi_squared = projection.projected_radius_squared(x_array, y_array)
        stable_discriminant = projection.alpha * (
            m_max * m_max - xi_squared
        )
        discriminant = np.where(
            needs_stable_fallback, stable_discriminant, raw_discriminant
        )
    else:
        discriminant = raw_discriminant

    inside = discriminant >= 0.0
    chord_center = -beta / projection.alpha
    half_length = np.zeros_like(discriminant, dtype=float)
    half_length[inside] = np.sqrt(np.maximum(discriminant[inside], 0.0)) / (
        projection.alpha
    )

    nodes, weights = leggauss(int(order))
    integral = np.zeros_like(discriminant, dtype=float)
    for node, weight in zip(nodes, weights, strict=True):
        los_coordinate = chord_center + half_length * node
        intrinsic_0 = sky_intrinsic_0 + intrinsic_los[0] * los_coordinate
        intrinsic_1 = sky_intrinsic_1 + intrinsic_los[1] * los_coordinate
        intrinsic_2 = sky_intrinsic_2 + intrinsic_los[2] * los_coordinate
        radius_squared = (
            intrinsic_0 * intrinsic_0
            + intrinsic_1 * intrinsic_1
            + intrinsic_2 * intrinsic_2
        )
        radius = np.sqrt(np.maximum(radius_squared, 0.0))
        # Do not call finite-domain profiles outside the projected boundary.
        evaluation_radius = np.where(inside, radius, 0.0)
        density_values = np.asarray(density(evaluation_radius), dtype=float)
        if density_values.shape != radius.shape:
            density_values = np.broadcast_to(density_values, radius.shape)
        if not np.all(np.isfinite(density_values[inside])):
            raise ValueError("density returned a non-finite value inside the ellipsoid")
        if np.any(density_values[inside] < 0.0):
            raise ValueError("density returned a negative value inside the ellipsoid")
        density_values = np.where(inside, density_values, 0.0)
        integral += weight * density_values

    result = half_length * integral
    return np.where(inside, result, 0.0)


def surface_density_abel(
    projection: ProjectionGeometry,
    x: ArrayLike,
    y: ArrayLike,
    density: DensityCallable,
    m_max: float | None = None,
    epsabs: float = 1.0e-11,
    epsrel: float = 1.0e-11,
) -> NDArray[np.float64]:
    r"""Project an arbitrary ``rho(m)`` through a nonsingular Abel integral.

    For finite support this evaluates

    ``Sigma(xi) = 2/sqrt(alpha) * integral_0^sqrt(m_max^2-xi^2)
    rho(sqrt(xi^2+t^2)) dt``.

    Unlike :func:`surface_density_los`, this routine explicitly uses the
    reduction to projected elliptical radius and is intended for accurate
    one-dimensional reference calculations rather than large raster maps.
    """
    if not callable(density):
        raise TypeError("density must be callable")
    m_max = _resolve_m_max(density, m_max)
    if not np.isfinite(epsabs) or not np.isfinite(epsrel) or epsabs < 0.0 or epsrel < 0.0:
        raise ValueError("integration tolerances must be nonnegative")

    xi = projection.projected_radius(x, y)
    result = np.zeros_like(xi, dtype=float)
    scale = 2.0 / np.sqrt(projection.alpha)

    iterator = np.nditer(
        [xi, result],
        flags=["refs_ok", "zerosize_ok"],
        op_flags=[["readonly"], ["writeonly"]],
    )
    for xi_item, result_item in iterator:
        xi_value = float(xi_item)
        if xi_value > m_max:
            result_item[...] = 0.0
            continue
        upper = np.sqrt(max(m_max * m_max - xi_value * xi_value, 0.0))

        def integrand(t: float) -> float:
            value = density(np.sqrt(xi_value * xi_value + t * t))
            scalar = float(np.asarray(value, dtype=float))
            if not np.isfinite(scalar) or scalar < 0.0:
                raise ValueError("density must return finite, nonnegative values")
            return scalar

        integral, _ = quad(
            integrand, 0.0, upper, epsabs=epsabs, epsrel=epsrel, limit=200
        )
        result_item[...] = scale * integral
    return result


@dataclass(frozen=True)
class SurfaceDensityMap:
    """A regular sky grid and its projected surface density."""

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    sigma: NDArray[np.float64]

    @property
    def mesh(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return 2D ``(X, Y)`` coordinate arrays."""
        return np.meshgrid(self.x, self.y, indexing="xy")


def surface_density_map(
    projection: ProjectionGeometry,
    density: DensityCallable,
    m_max: float | None = None,
    grid_size: int = 401,
    padding: float = 1.05,
    order: int = 64,
) -> SurfaceDensityMap:
    """Build a regular raster map using direct LOS quadrature."""
    if not isinstance(grid_size, (int, np.integer)) or grid_size < 25:
        raise ValueError("grid_size must be an integer of at least 25")
    resolved_m_max = _resolve_m_max(density, m_max)
    x_min, x_max, y_min, y_max = projection.plot_extent(resolved_m_max, padding)
    x = np.linspace(x_min, x_max, int(grid_size))
    y = np.linspace(y_min, y_max, int(grid_size))
    x_grid, y_grid = np.meshgrid(x, y, indexing="xy")
    sigma = surface_density_los(
        projection, x_grid, y_grid, density, m_max=resolved_m_max, order=order
    )
    return SurfaceDensityMap(x=x, y=y, sigma=sigma)
