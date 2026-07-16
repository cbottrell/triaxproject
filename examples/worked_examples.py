#!/usr/bin/env python3
"""Worked numerical projections for uniform and ellipsoidally radial density."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

# Some managed environments do not provide a writable default Matplotlib cache.
_mpl_cache = Path(tempfile.gettempdir()) / "triaxproject_matplotlib"
_mpl_cache.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_mpl_cache))
os.environ.setdefault("XDG_CACHE_HOME", str(_mpl_cache))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from triaxproject import (
    CoredPowerLawDensity,
    TriaxialEllipsoid,
    UniformDensity,
    cored_power_law_surface_density_analytic,
    measure_contours,
    surface_density_map,
    uniform_surface_density_analytic,
)


@dataclass
class ExampleResult:
    name: str
    density_map: object
    central_sigma: float
    xi_values: np.ndarray
    levels: np.ndarray
    fits: list
    normalized_map_error: float


def axial_angle_difference(first: float, second: float) -> float:
    """Signed angle difference for unoriented axes, returned in radians."""
    delta = first - second
    return 0.5 * np.arctan2(np.sin(2.0 * delta), np.cos(2.0 * delta))


def run_case(
    name,
    projection,
    profile,
    analytic_function,
    xi_values,
    grid_size,
    quadrature_order,
) -> ExampleResult:
    density_map = surface_density_map(
        projection,
        profile,
        grid_size=grid_size,
        padding=1.06,
        order=quadrature_order,
    )
    x_grid, y_grid = density_map.mesh
    analytic_map = analytic_function(projection, x_grid, y_grid)
    central_sigma = float(analytic_function(projection, 0.0, 0.0))
    levels = np.asarray(
        [
            float(
                analytic_function(
                    projection,
                    *projection.points_on_ellipse(float(xi), count=5)[0],
                )
            )
            for xi in xi_values
        ]
    )
    fits = measure_contours(density_map.x, density_map.y, density_map.sigma, levels)
    normalized_map_error = float(
        np.max(np.abs(density_map.sigma - analytic_map)) / central_sigma
    )
    return ExampleResult(
        name=name,
        density_map=density_map,
        central_sigma=central_sigma,
        xi_values=xi_values,
        levels=levels,
        fits=fits,
        normalized_map_error=normalized_map_error,
    )


def save_measurements(results, projection, output_path: Path) -> None:
    fieldnames = [
        "case",
        "elliptical_radius_xi",
        "sigma_level",
        "sigma_over_sigma0",
        "fitted_axis_ratio_b_over_a",
        "fitted_ellipticity_1_minus_b_over_a",
        "predicted_axis_ratio_b_over_a",
        "axis_ratio_error",
        "fitted_position_angle_deg",
        "predicted_position_angle_deg",
        "position_angle_error_deg",
        "fitted_center_x",
        "fitted_center_y",
        "center_offset",
        "rms_conic_residual",
        "contour_vertex_count",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for result in results:
            for xi, level, fit in zip(
                result.xi_values, result.levels, result.fits, strict=True
            ):
                angle_error = axial_angle_difference(
                    fit.position_angle, projection.position_angle
                )
                writer.writerow(
                    {
                        "case": result.name,
                        "elliptical_radius_xi": f"{xi:.8f}",
                        "sigma_level": f"{level:.12e}",
                        "sigma_over_sigma0": f"{level / result.central_sigma:.12e}",
                        "fitted_axis_ratio_b_over_a": f"{fit.axis_ratio:.12f}",
                        "fitted_ellipticity_1_minus_b_over_a": f"{fit.ellipticity:.12f}",
                        "predicted_axis_ratio_b_over_a": f"{projection.axis_ratio:.12f}",
                        "axis_ratio_error": f"{fit.axis_ratio - projection.axis_ratio:.12e}",
                        "fitted_position_angle_deg": f"{np.degrees(fit.position_angle):.9f}",
                        "predicted_position_angle_deg": (
                            f"{np.degrees(projection.position_angle):.9f}"
                        ),
                        "position_angle_error_deg": f"{np.degrees(angle_error):.12e}",
                        "fitted_center_x": f"{fit.center[0]:.12e}",
                        "fitted_center_y": f"{fit.center[1]:.12e}",
                        "center_offset": f"{np.linalg.norm(fit.center):.12e}",
                        "rms_conic_residual": f"{fit.rms_conic_residual:.12e}",
                        "contour_vertex_count": fit.n_points,
                    }
                )


def save_figure(results, projection, output_path: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12.0, 10.0), constrained_layout=True)
    colors = plt.cm.plasma(np.linspace(0.12, 0.88, len(results[0].xi_values)))

    for axis, result in zip(axes[0], results, strict=True):
        normalized = result.density_map.sigma / result.central_sigma
        mesh = axis.pcolormesh(
            result.density_map.x,
            result.density_map.y,
            normalized,
            shading="auto",
            cmap="magma",
            vmin=0.0,
            vmax=1.0,
        )
        for color, xi, level in zip(
            colors, result.xi_values, result.levels, strict=True
        ):
            # Solid: contour extracted from the raw numerical LOS raster.
            axis.contour(
                result.density_map.x,
                result.density_map.y,
                normalized,
                levels=[level / result.central_sigma],
                colors=[color],
                linewidths=2.0,
            )
            # Dashed black: ellipse predicted solely from the 3D geometry.
            points = projection.points_on_ellipse(float(xi))
            axis.plot(
                points[:, 0],
                points[:, 1],
                color="black",
                linestyle="--",
                linewidth=0.9,
                alpha=0.9,
            )
        axis.set_aspect("equal")
        axis.set_xlabel("sky x")
        axis.set_ylabel("sky y")
        title = result.name
        if result.name == "Cored radial density":
            title = r"Radial density $\rho(m)\propto[1+(m/0.25)^2]^{-1}$"
        axis.set_title(title)
        figure.colorbar(mesh, ax=axis, label=r"$\Sigma/\Sigma_0$", shrink=0.86)
        axis.text(
            0.03,
            0.04,
            "solid: numerical contour\ndashed: geometric ellipse",
            transform=axis.transAxes,
            color="white",
            fontsize=9,
            bbox={"facecolor": "black", "alpha": 0.45, "edgecolor": "none"},
        )

    axis_ratio_axis = axes[1, 0]
    pa_axis = axes[1, 1]
    markers = ("o", "s")
    for result, marker in zip(results, markers, strict=True):
        ratio_errors_ppm = 1.0e6 * np.array(
            [fit.axis_ratio - projection.axis_ratio for fit in result.fits]
        )
        angle_errors_arcsec = np.array(
            [
                3600.0
                * np.degrees(
                    axial_angle_difference(
                        fit.position_angle, projection.position_angle
                    )
                )
                for fit in result.fits
            ]
        )
        axis_ratio_axis.plot(
            result.xi_values,
            ratio_errors_ppm,
            marker=marker,
            linewidth=1.5,
            label=result.name,
        )
        pa_axis.plot(
            result.xi_values,
            angle_errors_arcsec,
            marker=marker,
            linewidth=1.5,
            label=result.name,
        )

    axis_ratio_axis.axhline(0.0, color="black", linestyle="--")
    axis_ratio_axis.set_xlabel(r"projected elliptical radius $\xi$")
    axis_ratio_axis.set_ylabel(r"$10^6\,[(b/a)_{\rm fit}-(b/a)_{\rm geom}]$")
    axis_ratio_axis.set_title(
        f"Axis-ratio residuals (geometric b/a = {projection.axis_ratio:.9f})"
    )
    axis_ratio_axis.grid(alpha=0.25)
    axis_ratio_axis.legend(fontsize=8)

    expected_angle = np.degrees(projection.position_angle)
    pa_axis.axhline(0.0, color="black", linestyle="--")
    pa_axis.set_xlabel(r"projected elliptical radius $\xi$")
    pa_axis.set_ylabel("fitted minus geometric PA (arcsec)")
    pa_axis.set_title(f"Orientation residuals (geometric PA = {expected_angle:.6f} deg)")
    pa_axis.grid(alpha=0.25)
    pa_axis.legend(fontsize=8)

    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_summary(results, ellipsoid, projection, output_path: Path) -> str:
    lines = [
        "TRIAXIAL PROJECTION WORKED EXAMPLES",
        "",
        f"Intrinsic semiaxes (m=1): {tuple(float(v) for v in ellipsoid.axes)}",
        "Normalized LOS: " + np.array2string(projection.line_of_sight, precision=8),
        f"Predicted projected b/a: {projection.axis_ratio:.12f}",
        f"Predicted ellipticity 1-b/a: {projection.ellipticity:.12f}",
        f"Predicted position angle: {np.degrees(projection.position_angle):.9f} deg",
        "",
        "case                       xi      fitted b/a       delta(b/a)       "
        "center offset     PA error (deg)",
    ]
    for result in results:
        for xi, fit in zip(result.xi_values, result.fits, strict=True):
            pa_error = np.degrees(
                axial_angle_difference(fit.position_angle, projection.position_angle)
            )
            lines.append(
                f"{result.name[:25]:25s}  {xi:5.2f}   {fit.axis_ratio:14.10f}  "
                f"{fit.axis_ratio - projection.axis_ratio:+14.3e}  "
                f"{np.linalg.norm(fit.center):14.3e}  {pa_error:+14.3e}"
            )
        recovered = np.array([fit.axis_ratio for fit in result.fits])
        lines.extend(
            [
                f"  {result.name} b/a range: {np.ptp(recovered):.3e}",
                f"  {result.name} max |raw LOS map - analytic map| / Sigma0: "
                f"{result.normalized_map_error:.3e}",
                "",
            ]
        )
    text = "\n".join(lines).rstrip() + "\n"
    output_path.write_text(text, encoding="utf-8")
    return text


def parse_args() -> argparse.Namespace:
    package_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=package_root / "example_outputs",
        help="directory for the PNG, CSV, and text summary",
    )
    parser.add_argument("--grid-size", type=int, default=501)
    parser.add_argument("--quadrature-order", type=int, default=48)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ellipsoid = TriaxialEllipsoid(axes=(1.8, 1.1, 0.55))
    projection = ellipsoid.project(line_of_sight=(0.43, 0.58, 0.69))
    xi_values = np.array((0.25, 0.40, 0.55, 0.70, 0.82, 0.90))

    uniform = run_case(
        "Uniform density",
        projection,
        UniformDensity(rho0=1.0),
        lambda p, x, y: uniform_surface_density_analytic(p, x, y, rho0=1.0),
        xi_values,
        args.grid_size,
        args.quadrature_order,
    )
    radial = run_case(
        "Cored radial density",
        projection,
        CoredPowerLawDensity(rho0=1.0, core_radius=0.25, slope=2.0),
        lambda p, x, y: cored_power_law_surface_density_analytic(
            p, x, y, rho0=1.0, core_radius=0.25, m_max=1.0
        ),
        xi_values,
        args.grid_size,
        args.quadrature_order,
    )
    results = (uniform, radial)

    save_measurements(results, projection, args.output_dir / "contour_measurements.csv")
    save_figure(results, projection, args.output_dir / "constant_ellipticity.png")
    summary = save_summary(
        results, ellipsoid, projection, args.output_dir / "example_summary.txt"
    )
    print(summary, end="")
    print(f"Saved outputs in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
