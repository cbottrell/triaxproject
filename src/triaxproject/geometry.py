"""Geometry of a triaxial ellipsoid under an arbitrary parallel projection."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _vector3(value: ArrayLike, name: str) -> NDArray[np.float64]:
    vector = np.asarray(value, dtype=float)
    if vector.shape != (3,):
        raise ValueError(f"{name} must contain exactly three values")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite values")
    return vector


def _sky_basis(
    line_of_sight: NDArray[np.float64], sky_x: ArrayLike | None
) -> NDArray[np.float64]:
    """Return columns (e_x, e_y) with e_x cross e_y = line_of_sight."""
    if sky_x is None:
        reference = np.eye(3)[np.argmin(np.abs(line_of_sight))]
    else:
        reference = _vector3(sky_x, "sky_x")

    e_x = reference - np.dot(reference, line_of_sight) * line_of_sight
    norm = np.linalg.norm(e_x)
    if norm <= 1.0e-14:
        raise ValueError("sky_x must not be parallel to the line of sight")
    e_x /= norm
    e_y = np.cross(line_of_sight, e_x)
    return np.column_stack((e_x, e_y))


@dataclass(frozen=True)
class TriaxialEllipsoid:
    """A family of similar ellipsoids described by ``m^2 = dx.T @ A @ dx``.

    The shell ``m=1`` has principal semiaxes ``axes``.  The columns of
    ``orientation`` are those principal-axis directions in world coordinates.
    Sky coordinates returned by :meth:`project` are centered on ``center``.
    """

    axes: ArrayLike
    orientation: ArrayLike | None = None
    center: ArrayLike = (0.0, 0.0, 0.0)

    def __post_init__(self) -> None:
        axes = _vector3(self.axes, "axes")
        if np.any(axes <= 0.0):
            raise ValueError("all semiaxes must be strictly positive")

        if self.orientation is None:
            orientation = np.eye(3)
        else:
            orientation = np.asarray(self.orientation, dtype=float)
            if orientation.shape != (3, 3):
                raise ValueError("orientation must be a 3 by 3 matrix")
            if not np.all(np.isfinite(orientation)):
                raise ValueError("orientation must contain only finite values")
            if not np.allclose(
                orientation.T @ orientation, np.eye(3), rtol=0.0, atol=1.0e-12
            ):
                raise ValueError("orientation must be orthonormal")
            if not np.isclose(np.linalg.det(orientation), 1.0, atol=1.0e-12):
                raise ValueError("orientation must be a proper rotation (determinant +1)")

        center = _vector3(self.center, "center")
        object.__setattr__(self, "axes", axes.copy())
        object.__setattr__(self, "orientation", orientation.copy())
        object.__setattr__(self, "center", center.copy())

    @property
    def quadratic_form(self) -> NDArray[np.float64]:
        """The positive-definite matrix ``A`` defining ellipsoidal radius."""
        inverse_axes_squared = np.diag(1.0 / np.square(self.axes))
        return self.orientation @ inverse_axes_squared @ self.orientation.T

    @property
    def inverse_quadratic_form(self) -> NDArray[np.float64]:
        """The inverse of ``A``, formed without numerically inverting ``A``."""
        axes_squared = np.diag(np.square(self.axes))
        return self.orientation @ axes_squared @ self.orientation.T

    def ellipsoidal_radius(self, points: ArrayLike) -> NDArray[np.float64]:
        """Evaluate ``m`` at world-coordinate points with final dimension 3."""
        points_array = np.asarray(points, dtype=float)
        if points_array.shape[-1:] != (3,):
            raise ValueError("points must have final dimension 3")
        displacement = points_array - self.center
        radius_squared = np.einsum(
            "...i,ij,...j->...", displacement, self.quadratic_form, displacement
        )
        return np.sqrt(np.maximum(radius_squared, 0.0))

    def project(
        self, line_of_sight: ArrayLike, sky_x: ArrayLike | None = None
    ) -> "ProjectionGeometry":
        """Construct the orthographic projection geometry for a viewing direction."""
        return ProjectionGeometry.from_ellipsoid(self, line_of_sight, sky_x=sky_x)


@dataclass(frozen=True)
class ProjectionGeometry:
    """Precomputed matrices for one orthographic projection.

    The line of sight points toward increasing LOS coordinate. ``sky_basis``
    contains the world-coordinate sky x/y unit vectors as its columns.
    """

    ellipsoid: TriaxialEllipsoid
    line_of_sight: NDArray[np.float64]
    sky_basis: NDArray[np.float64]
    alpha: float
    coupling: NDArray[np.float64]
    sky_matrix: NDArray[np.float64]
    matrix: NDArray[np.float64]

    @classmethod
    def from_ellipsoid(
        cls,
        ellipsoid: TriaxialEllipsoid,
        line_of_sight: ArrayLike,
        sky_x: ArrayLike | None = None,
    ) -> "ProjectionGeometry":
        los = _vector3(line_of_sight, "line_of_sight")
        norm = np.linalg.norm(los)
        if norm <= 0.0:
            raise ValueError("line_of_sight must be nonzero")
        los /= norm
        basis = _sky_basis(los, sky_x)

        quadratic_form = ellipsoid.quadratic_form
        alpha = float(los @ quadratic_form @ los)
        coupling = basis.T @ quadratic_form @ los
        sky_matrix = basis.T @ quadratic_form @ basis
        # H = (P.T A^{-1} P)^{-1} is algebraically the same Schur complement
        # as C - d d.T/alpha, but avoids catastrophic cancellation for highly
        # elongated ellipsoids.
        projected_covariance = basis.T @ ellipsoid.inverse_quadratic_form @ basis
        projected_matrix = np.linalg.solve(projected_covariance, np.eye(2))
        projected_matrix = 0.5 * (projected_matrix + projected_matrix.T)
        if np.any(np.linalg.eigvalsh(projected_matrix) <= 0.0):
            raise ArithmeticError("projected quadratic form is not positive definite")

        return cls(
            ellipsoid=ellipsoid,
            line_of_sight=los.copy(),
            sky_basis=basis.copy(),
            alpha=alpha,
            coupling=coupling.copy(),
            sky_matrix=sky_matrix.copy(),
            matrix=projected_matrix.copy(),
        )

    def projected_radius_squared(
        self, x: ArrayLike, y: ArrayLike
    ) -> NDArray[np.float64]:
        """Evaluate ``xi^2 = (x,y).T @ matrix @ (x,y)`` with broadcasting."""
        x_array, y_array = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        )
        if not np.all(np.isfinite(x_array)) or not np.all(np.isfinite(y_array)):
            raise ValueError("sky coordinates must contain only finite values")
        xi_squared = (
            self.matrix[0, 0] * x_array * x_array
            + 2.0 * self.matrix[0, 1] * x_array * y_array
            + self.matrix[1, 1] * y_array * y_array
        )
        return np.maximum(xi_squared, 0.0)

    def projected_radius(self, x: ArrayLike, y: ArrayLike) -> NDArray[np.float64]:
        """Evaluate the dimensionless projected elliptical radius ``xi``."""
        return np.sqrt(self.projected_radius_squared(x, y))

    def line_of_sight_center(
        self, x: ArrayLike, y: ArrayLike
    ) -> NDArray[np.float64]:
        """LOS coordinate at the midpoint of the ellipsoid chord."""
        x_array, y_array = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        )
        beta = self.coupling[0] * x_array + self.coupling[1] * y_array
        return -beta / self.alpha

    def chord_bounds(
        self, x: ArrayLike, y: ArrayLike, m_max: float = 1.0
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return the two LOS intersections with the shell ``m=m_max``.

        Outside the projected boundary both returned arrays contain ``nan``.
        """
        if not np.isfinite(m_max) or m_max <= 0.0:
            raise ValueError("m_max must be finite and strictly positive")
        xi_squared = self.projected_radius_squared(x, y)
        center = self.line_of_sight_center(x, y)
        inside = xi_squared <= m_max * m_max
        half_length = np.full_like(xi_squared, np.nan, dtype=float)
        half_length[inside] = np.sqrt(
            np.maximum(m_max * m_max - xi_squared[inside], 0.0) / self.alpha
        )
        return center - half_length, center + half_length

    @property
    def eigenvalues(self) -> NDArray[np.float64]:
        """Ascending eigenvalues of the projected quadratic form."""
        return np.linalg.eigvalsh(self.matrix)

    @property
    def principal_directions(self) -> NDArray[np.float64]:
        """Sky-coordinate eigenvectors, with the projected major axis first."""
        _, eigenvectors = np.linalg.eigh(self.matrix)
        return eigenvectors

    @property
    def axis_ratio(self) -> float:
        """Projected minor-to-major axis ratio ``b/a``."""
        eigenvalues = self.eigenvalues
        return float(np.sqrt(eigenvalues[0] / eigenvalues[1]))

    @property
    def ellipticity(self) -> float:
        """Projected ellipticity using the convention ``epsilon = 1 - b/a``."""
        return 1.0 - self.axis_ratio

    @property
    def position_angle(self) -> float:
        """Major-axis angle in radians from sky x, reduced modulo pi."""
        major_direction = self.principal_directions[:, 0]
        return float(np.mod(np.arctan2(major_direction[1], major_direction[0]), np.pi))

    def semiaxes(self, xi: float = 1.0) -> tuple[float, float]:
        """Projected major/minor semiaxes of the contour at radius ``xi``."""
        if not np.isfinite(xi) or xi < 0.0:
            raise ValueError("xi must be finite and nonnegative")
        eigenvalues = self.eigenvalues
        return float(xi / np.sqrt(eigenvalues[0])), float(
            xi / np.sqrt(eigenvalues[1])
        )

    def points_on_ellipse(self, xi: float, count: int = 361) -> NDArray[np.float64]:
        """Sample a projected ``xi=constant`` ellipse in sky coordinates."""
        if count < 4:
            raise ValueError("count must be at least 4")
        major, minor = self.semiaxes(xi)
        angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=True)
        principal_points = np.vstack((major * np.cos(angle), minor * np.sin(angle)))
        return (self.principal_directions @ principal_points).T

    def plot_extent(
        self, m_max: float = 1.0, padding: float = 1.05
    ) -> tuple[float, float, float, float]:
        """Axis-aligned sky bounds enclosing the projected ``m_max`` ellipse."""
        if not np.isfinite(m_max) or m_max <= 0.0:
            raise ValueError("m_max must be finite and strictly positive")
        if not np.isfinite(padding) or padding < 1.0:
            raise ValueError("padding must be finite and at least 1")
        inverse = np.linalg.inv(self.matrix)
        x_limit = padding * m_max * np.sqrt(inverse[0, 0])
        y_limit = padding * m_max * np.sqrt(inverse[1, 1])
        return -x_limit, x_limit, -y_limit, y_limit
