from __future__ import annotations

import unittest
from dataclasses import replace

import numpy as np

from triaxproject import (
    CoredPowerLawDensity,
    PolynomialDensity,
    TriaxialEllipsoid,
    UniformDensity,
    cored_power_law_surface_density_analytic,
    fit_ellipse,
    measure_contours,
    polynomial_surface_density_analytic,
    surface_density_abel,
    surface_density_los,
    surface_density_map,
    uniform_surface_density_analytic,
)


class GeometryTests(unittest.TestCase):
    def test_principal_axis_view_has_known_projected_axes(self) -> None:
        projection = TriaxialEllipsoid((3.0, 2.0, 1.0)).project((0.0, 0.0, 1.0))
        self.assertTrue(np.allclose(projection.semiaxes(), (3.0, 2.0), atol=1e-14))
        self.assertAlmostEqual(projection.axis_ratio, 2.0 / 3.0, places=14)
        self.assertAlmostEqual(projection.ellipticity, 1.0 / 3.0, places=14)

    def test_arbitrary_view_has_positive_projected_form(self) -> None:
        projection = TriaxialEllipsoid((1.8, 1.1, 0.55)).project((0.43, 0.58, 0.69))
        self.assertTrue(np.all(projection.eigenvalues > 0.0))
        self.assertTrue(np.allclose(projection.sky_basis.T @ projection.sky_basis, np.eye(2)))
        self.assertTrue(
            np.allclose(projection.sky_basis.T @ projection.line_of_sight, 0.0)
        )
        self.assertAlmostEqual(np.linalg.norm(projection.line_of_sight), 1.0)

    def test_rotation_covariance(self) -> None:
        angle = 0.47
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        axes = (2.1, 1.2, 0.7)
        los = np.array((0.2, -0.5, 0.84))
        first = TriaxialEllipsoid(axes).project(los)
        second = TriaxialEllipsoid(axes, orientation=rotation).project(rotation @ los)
        self.assertTrue(np.allclose(first.eigenvalues, second.eigenvalues, atol=1e-13))
        self.assertAlmostEqual(first.axis_ratio, second.axis_ratio, places=13)

    def test_sphere_projects_to_circle(self) -> None:
        for los in ((1.0, 0.0, 0.0), (1.0, 2.0, -3.0), (-0.2, 0.9, 0.4)):
            projection = TriaxialEllipsoid((1.7, 1.7, 1.7)).project(los)
            self.assertAlmostEqual(projection.axis_ratio, 1.0, places=13)

    def test_high_dynamic_range_projection_remains_positive(self) -> None:
        projection = TriaxialEllipsoid((1.0e5, 1.0, 1.0e-5)).project(
            (0.43, 0.58, 0.69)
        )
        self.assertTrue(np.all(projection.eigenvalues > 0.0))
        points = projection.points_on_ellipse(0.2, count=65)[:-1]
        sigma = surface_density_los(
            projection,
            points[:, 0],
            points[:, 1],
            PolynomialDensity(power=2.0),
            order=32,
        )
        self.assertLess(np.ptp(sigma) / np.mean(sigma), 1.0e-6)

    def test_invalid_geometry_inputs(self) -> None:
        with self.assertRaises(ValueError):
            TriaxialEllipsoid((1.0, 0.0, 2.0))
        with self.assertRaises(ValueError):
            TriaxialEllipsoid((1.0, 2.0, 3.0)).project((0.0, 0.0, 0.0))
        with self.assertRaises(ValueError):
            TriaxialEllipsoid((1.0, 2.0, 3.0), orientation=np.ones((3, 3)))


class ProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.projection = TriaxialEllipsoid((1.8, 1.1, 0.55)).project(
            (0.43, 0.58, 0.69)
        )

    def _sample_points(self) -> tuple[np.ndarray, np.ndarray]:
        pieces = []
        for xi in (0.0, 0.15, 0.37, 0.61, 0.83, 0.95):
            pieces.append(self.projection.points_on_ellipse(xi, count=19)[:-1])
        points = np.vstack(pieces)
        return points[:, 0], points[:, 1]

    def test_uniform_raw_los_matches_closed_form(self) -> None:
        x, y = self._sample_points()
        numerical = surface_density_los(
            self.projection, x, y, UniformDensity(2.3), order=32
        )
        analytic = uniform_surface_density_analytic(
            self.projection, x, y, rho0=2.3
        )
        self.assertTrue(np.allclose(numerical, analytic, rtol=3e-14, atol=3e-14))

    def test_raw_chords_do_not_depend_on_projected_matrix_for_regular_geometry(self) -> None:
        points = self.projection.points_on_ellipse(0.55, count=17)[:-1]
        corrupted_projection = replace(
            self.projection, matrix=1.7 * self.projection.matrix
        )
        reference = surface_density_los(
            self.projection,
            points[:, 0],
            points[:, 1],
            UniformDensity(),
            order=16,
        )
        corrupted = surface_density_los(
            corrupted_projection,
            points[:, 0],
            points[:, 1],
            UniformDensity(),
            order=16,
        )
        self.assertTrue(np.allclose(corrupted, reference, rtol=1e-14, atol=1e-14))

    def test_polynomial_raw_los_matches_closed_form(self) -> None:
        x, y = self._sample_points()
        profile = PolynomialDensity(rho0=1.7, power=2.0)
        numerical = surface_density_los(
            self.projection, x, y, profile, order=32
        )
        analytic = polynomial_surface_density_analytic(
            self.projection, x, y, rho0=1.7, power=2.0
        )
        self.assertTrue(np.allclose(numerical, analytic, rtol=2e-13, atol=2e-13))

    def test_abel_matches_independent_raw_los(self) -> None:
        profile = lambda m: np.exp(-2.4 * np.asarray(m) ** 2)
        points = np.vstack(
            [
                self.projection.points_on_ellipse(xi, count=11)[:-1]
                for xi in (0.0, 0.2, 0.57, 0.86, 0.97)
            ]
        )
        raw = surface_density_los(
            self.projection, points[:, 0], points[:, 1], profile, order=128
        )
        abel = surface_density_abel(
            self.projection, points[:, 0], points[:, 1], profile
        )
        self.assertTrue(np.allclose(raw, abel, rtol=2e-12, atol=2e-12))

    def test_cored_raw_los_matches_nonpolynomial_closed_form(self) -> None:
        x, y = self._sample_points()
        profile = CoredPowerLawDensity(rho0=1.4, core_radius=0.27, slope=2.0)
        numerical = surface_density_los(
            self.projection, x, y, profile, order=64
        )
        analytic = cored_power_law_surface_density_analytic(
            self.projection, x, y, rho0=1.4, core_radius=0.27
        )
        self.assertTrue(np.allclose(numerical, analytic, rtol=2e-13, atol=2e-13))

    def test_density_is_constant_around_each_predicted_ellipse(self) -> None:
        profiles = (UniformDensity(), PolynomialDensity(power=2.0))
        for profile in profiles:
            for xi in (0.2, 0.43, 0.7, 0.9):
                points = self.projection.points_on_ellipse(xi, count=129)[:-1]
                sigma = surface_density_los(
                    self.projection,
                    points[:, 0],
                    points[:, 1],
                    profile,
                    order=48,
                )
                fractional_range = np.ptp(sigma) / np.mean(sigma)
                self.assertLess(fractional_range, 2e-13)

    def test_random_views_and_orientations_keep_constant_contours(self) -> None:
        rng = np.random.default_rng(1847)
        profile = lambda m: np.exp(-1.7 * np.asarray(m) ** 1.35)
        m_max = 1.3
        for _ in range(8):
            trial_matrix = rng.normal(size=(3, 3))
            orientation, _ = np.linalg.qr(trial_matrix)
            if np.linalg.det(orientation) < 0.0:
                orientation[:, 0] *= -1.0
            line_of_sight = rng.normal(size=3)
            projection = TriaxialEllipsoid(
                (2.0, 1.15, 0.52), orientation=orientation
            ).project(line_of_sight)
            for fraction in (0.18, 0.47, 0.79, 0.93):
                points = projection.points_on_ellipse(
                    fraction * m_max, count=65
                )[:-1]
                sigma = surface_density_los(
                    projection,
                    points[:, 0],
                    points[:, 1],
                    profile,
                    m_max=m_max,
                    order=64,
                )
                self.assertLess(np.ptp(sigma) / np.mean(sigma), 2.0e-11)

    def test_scalar_and_outside_inputs(self) -> None:
        center = surface_density_los(
            self.projection, 0.0, 0.0, UniformDensity(), order=16
        )
        outside = surface_density_los(
            self.projection, 100.0, 100.0, UniformDensity(), order=16
        )
        self.assertGreater(float(center), 0.0)
        self.assertEqual(float(outside), 0.0)
        with self.assertRaises(ValueError):
            uniform_surface_density_analytic(self.projection, np.nan, 0.0)

    def test_abel_rejects_unphysical_density_values(self) -> None:
        with self.assertRaises(ValueError):
            surface_density_abel(
                self.projection, 0.2, 0.1, lambda m: -np.ones_like(m)
            )
        with self.assertRaises(ValueError):
            surface_density_abel(
                self.projection, 0.2, 0.1, lambda m: np.full_like(m, np.nan)
            )

    def test_density_callable_is_not_evaluated_outside_support(self) -> None:
        def finite_domain_profile(m):
            radius = np.asarray(m)
            if np.any(radius > 1.0 + 1.0e-12):
                raise RuntimeError("profile evaluated beyond its support")
            return np.sqrt(np.maximum(1.0 - radius * radius, 0.0))

        values = surface_density_los(
            self.projection,
            np.array((0.0, 100.0)),
            np.array((0.0, 100.0)),
            finite_domain_profile,
            order=24,
        )
        self.assertGreater(values[0], 0.0)
        self.assertEqual(values[1], 0.0)

    def test_profile_support_is_inferred_and_mismatch_is_rejected(self) -> None:
        profile = PolynomialDensity(power=1.0, m_max=1.7)
        inferred = surface_density_los(
            self.projection, 0.0, 0.0, profile, order=24
        )
        explicit = surface_density_los(
            self.projection, 0.0, 0.0, profile, m_max=1.7, order=24
        )
        self.assertAlmostEqual(float(inferred), float(explicit), places=13)
        with self.assertRaises(ValueError):
            surface_density_los(
                self.projection, 0.0, 0.0, profile, m_max=1.0, order=24
            )

    def test_raster_contours_recover_one_axis_ratio(self) -> None:
        profile = PolynomialDensity(power=2.0)
        density_map = surface_density_map(
            self.projection, profile, grid_size=301, padding=1.06, order=32
        )
        xi_values = np.array((0.25, 0.45, 0.65, 0.82, 0.90))
        central = float(
            polynomial_surface_density_analytic(
                self.projection, 0.0, 0.0, power=2.0
            )
        )
        contour_levels = central * np.power(1.0 - xi_values**2, 2.5)
        fits = measure_contours(
            density_map.x, density_map.y, density_map.sigma, contour_levels
        )
        recovered = np.array([fit.axis_ratio for fit in fits])
        centers = np.array([np.linalg.norm(fit.center) for fit in fits])
        angle_errors = np.array(
            [
                abs(
                    0.5
                    * np.arctan2(
                        np.sin(2.0 * (fit.position_angle - self.projection.position_angle)),
                        np.cos(2.0 * (fit.position_angle - self.projection.position_angle)),
                    )
                )
                for fit in fits
            ]
        )
        self.assertLess(np.max(np.abs(recovered - self.projection.axis_ratio)), 5e-4)
        self.assertLess(np.ptp(recovered), 5e-4)
        self.assertLess(np.max(centers), 2e-4)
        self.assertLess(np.max(angle_errors), np.radians(0.01))

    def test_ellipse_fit_supports_a_translated_contour(self) -> None:
        angle = np.linspace(0.0, 2.0 * np.pi, 257)
        rotation_angle = 0.37
        rotation = np.array(
            [
                [np.cos(rotation_angle), -np.sin(rotation_angle)],
                [np.sin(rotation_angle), np.cos(rotation_angle)],
            ]
        )
        points = np.column_stack((2.4 * np.cos(angle), 0.8 * np.sin(angle)))
        points = points @ rotation.T + np.array((3.0, -4.0))
        fit = fit_ellipse(points)
        self.assertTrue(np.allclose(fit.center, (3.0, -4.0), atol=1e-12))
        self.assertAlmostEqual(fit.axis_ratio, 1.0 / 3.0, places=12)


if __name__ == "__main__":
    unittest.main()
