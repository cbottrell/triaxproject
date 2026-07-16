"""Recover ellipse parameters from contours of a numerical raster map."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class EllipseFit:
    """Parameters of a free-centered ellipse fitted to contour vertices."""

    center: NDArray[np.float64]
    semimajor: float
    semiminor: float
    axis_ratio: float
    ellipticity: float
    position_angle: float
    quadratic_form: NDArray[np.float64]
    rms_conic_residual: float
    n_points: int
    level: float | None = None


def fit_ellipse(points: ArrayLike, level: float | None = None) -> EllipseFit:
    r"""Fit ``(x-c).T Q (x-c)=1`` to a set of contour vertices.

    The fit first solves the free-centered linear conic
    ``q11*x^2 + 2*q12*x*y + q22*y^2 + l1*x + l2*y = 1`` and then converts it
    to centered form. This is accurate for densely sampled, nearly elliptical
    contours and exposes both axis-ratio and concentricity errors.
    """
    vertices = np.asarray(points, dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")
    if vertices.shape[0] < 8:
        raise ValueError("at least eight contour points are required")
    if not np.all(np.isfinite(vertices)):
        raise ValueError("points must contain only finite values")

    # Translating near the data centroid makes the fixed-RHS conic fit valid
    # and well-conditioned even when the ellipse is far from coordinate zero.
    origin = np.mean(vertices, axis=0)
    local_vertices = vertices - origin
    x = local_vertices[:, 0]
    y = local_vertices[:, 1]
    design = np.column_stack((x * x, 2.0 * x * y, y * y, x, y))
    coefficients, _, rank, _ = np.linalg.lstsq(
        design, np.ones(vertices.shape[0]), rcond=None
    )
    if rank < 5:
        raise ValueError("contour points do not determine a free-centered ellipse")

    raw_quadratic = np.array(
        [
            [coefficients[0], coefficients[1]],
            [coefficients[1], coefficients[2]],
        ]
    )
    linear = coefficients[3:5]
    try:
        local_center = -0.5 * np.linalg.solve(raw_quadratic, linear)
    except np.linalg.LinAlgError as error:
        raise ValueError("fitted conic is singular") from error

    scale = 1.0 + local_center @ raw_quadratic @ local_center
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("fitted conic is not a real ellipse")
    quadratic = 0.5 * (
        raw_quadratic / scale + (raw_quadratic / scale).T
    )
    eigenvalues, eigenvectors = np.linalg.eigh(quadratic)
    if np.any(eigenvalues <= 0.0):
        raise ValueError("fitted conic is not a positive-definite ellipse")

    semimajor = float(1.0 / np.sqrt(eigenvalues[0]))
    semiminor = float(1.0 / np.sqrt(eigenvalues[1]))
    axis_ratio = semiminor / semimajor
    major_direction = eigenvectors[:, 0]
    position_angle = float(
        np.mod(np.arctan2(major_direction[1], major_direction[0]), np.pi)
    )
    center = origin + local_center
    centered = vertices - center
    conic_values = np.einsum("...i,ij,...j->...", centered, quadratic, centered)
    rms = float(np.sqrt(np.mean(np.square(conic_values - 1.0))))

    return EllipseFit(
        center=center,
        semimajor=semimajor,
        semiminor=semiminor,
        axis_ratio=axis_ratio,
        ellipticity=1.0 - axis_ratio,
        position_angle=position_angle,
        quadratic_form=quadratic,
        rms_conic_residual=rms,
        n_points=vertices.shape[0],
        level=None if level is None else float(level),
    )


def _polyline_length(points: NDArray[np.float64]) -> float:
    if points.shape[0] < 2:
        return 0.0
    return float(np.sum(np.linalg.norm(np.diff(points, axis=0), axis=1)))


def measure_contour(
    x: ArrayLike, y: ArrayLike, z: ArrayLike, level: float
) -> EllipseFit:
    """Extract the longest contour loop from a raster and fit its ellipse."""
    return measure_contours(x, y, z, [level])[0]


def measure_contours(
    x: ArrayLike, y: ArrayLike, z: ArrayLike, levels: Iterable[float]
) -> list[EllipseFit]:
    """Extract and fit one closed ellipse at each requested contour level."""
    try:
        import contourpy
    except ImportError as error:
        raise ImportError(
            "contour measurement requires contourpy; install the 'examples' extra"
        ) from error

    x_array = np.asarray(x, dtype=float)
    y_array = np.asarray(y, dtype=float)
    z_array = np.asarray(z, dtype=float)
    if x_array.ndim != 1 or y_array.ndim != 1:
        raise ValueError("x and y must be one-dimensional grid coordinates")
    if z_array.shape != (y_array.size, x_array.size):
        raise ValueError("z must have shape (len(y), len(x))")
    if not np.all(np.isfinite(z_array)):
        raise ValueError("z must contain only finite values")

    generator = contourpy.contour_generator(
        x=x_array, y=y_array, z=z_array, name="serial"
    )
    fits: list[EllipseFit] = []
    for level in levels:
        level_value = float(level)
        if not np.isfinite(level_value):
            raise ValueError("contour levels must be finite")
        segments = generator.lines(level_value)
        span = max(float(np.ptp(x_array)), float(np.ptp(y_array)), 1.0)
        closure_tolerance = 1.0e-10 * span
        edge_tolerance = 1.0e-10 * span
        usable = []
        for segment in segments:
            if segment.shape[0] < 8:
                continue
            if np.linalg.norm(segment[0] - segment[-1]) > closure_tolerance:
                continue
            touches_edge = (
                np.any(segment[:, 0] <= x_array[0] + edge_tolerance)
                or np.any(segment[:, 0] >= x_array[-1] - edge_tolerance)
                or np.any(segment[:, 1] <= y_array[0] + edge_tolerance)
                or np.any(segment[:, 1] >= y_array[-1] - edge_tolerance)
            )
            if not touches_edge:
                usable.append(segment)
        if not usable:
            raise ValueError(
                f"no closed, unclipped contour found at level {level_value:g}"
            )
        longest = max(usable, key=_polyline_length)
        fits.append(fit_ellipse(longest, level=level_value))
    return fits
