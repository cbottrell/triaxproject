"""Projection tools for homeoidally stratified triaxial ellipsoids."""

from .diagnostics import EllipseFit, fit_ellipse, measure_contour, measure_contours
from .geometry import ProjectionGeometry, TriaxialEllipsoid
from .profiles import CoredPowerLawDensity, PolynomialDensity, UniformDensity
from .projection import (
    SurfaceDensityMap,
    cored_power_law_surface_density_analytic,
    polynomial_surface_density_analytic,
    surface_density_abel,
    surface_density_los,
    surface_density_map,
    uniform_surface_density_analytic,
)

__all__ = [
    "CoredPowerLawDensity",
    "EllipseFit",
    "PolynomialDensity",
    "ProjectionGeometry",
    "SurfaceDensityMap",
    "TriaxialEllipsoid",
    "UniformDensity",
    "cored_power_law_surface_density_analytic",
    "fit_ellipse",
    "measure_contour",
    "measure_contours",
    "polynomial_surface_density_analytic",
    "surface_density_abel",
    "surface_density_los",
    "surface_density_map",
    "uniform_surface_density_analytic",
]

__version__ = "0.1.0"
