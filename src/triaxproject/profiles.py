"""Example density profiles expressed as functions of ellipsoidal radius."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


def _positive_finite(value: float, name: str) -> float:
    converted = float(value)
    if not np.isfinite(converted) or converted <= 0.0:
        raise ValueError(f"{name} must be finite and strictly positive")
    return converted


@dataclass(frozen=True)
class UniformDensity:
    """Constant volume density. The projection routine supplies the truncation."""

    rho0: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "rho0", _positive_finite(self.rho0, "rho0"))

    def __call__(self, m: ArrayLike) -> NDArray[np.float64]:
        radius = np.asarray(m, dtype=float)
        return np.full_like(radius, self.rho0, dtype=float)


@dataclass(frozen=True)
class PolynomialDensity:
    r"""Finite profile ``rho0 * (1 - (m/m_max)^2)^power`` inside ``m_max``."""

    rho0: float = 1.0
    power: float = 2.0
    m_max: float = 1.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "rho0", _positive_finite(self.rho0, "rho0"))
        object.__setattr__(self, "m_max", _positive_finite(self.m_max, "m_max"))
        power = float(self.power)
        if not np.isfinite(power) or power < 0.0:
            raise ValueError("power must be finite and nonnegative")
        object.__setattr__(self, "power", power)

    def __call__(self, m: ArrayLike) -> NDArray[np.float64]:
        radius = np.asarray(m, dtype=float)
        base = np.maximum(1.0 - np.square(radius / self.m_max), 0.0)
        density = self.rho0 * np.power(base, self.power)
        return np.where(radius <= self.m_max, density, 0.0)


@dataclass(frozen=True)
class CoredPowerLawDensity:
    r"""Smooth profile ``rho0 * (1 + (m/core_radius)^2)^(-slope/2)``."""

    rho0: float = 1.0
    core_radius: float = 0.2
    slope: float = 2.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "rho0", _positive_finite(self.rho0, "rho0"))
        object.__setattr__(
            self, "core_radius", _positive_finite(self.core_radius, "core_radius")
        )
        slope = float(self.slope)
        if not np.isfinite(slope) or slope < 0.0:
            raise ValueError("slope must be finite and nonnegative")
        object.__setattr__(self, "slope", slope)

    def __call__(self, m: ArrayLike) -> NDArray[np.float64]:
        radius = np.asarray(m, dtype=float)
        return self.rho0 * np.power(
            1.0 + np.square(radius / self.core_radius), -0.5 * self.slope
        )
