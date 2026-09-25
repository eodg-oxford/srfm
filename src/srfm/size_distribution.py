r"""Analytic particle-radius number distributions used by SRFM.

- Name: size_distribution.py
- Parent package: srfm
- Author: Don Grainger
- Contributors: Antonin Knizek
- Date: 3 January 2025
- Changelog:
    - Added Gaussian, multimode log-normal, gamma, modified gamma, inverse
      modified gamma, and regularised power-law distributions. AK. 25 Sep 2026

Units:
    - Particle size is expressed in :math:`\mu`\ m.
    - The distribution (n) is in number per :math:`\mu`\ m per cm\ :sup:`3`.
Integrated values are:
    - Number density is in number per cm\ :sup:`3`.
    - Surface area density is in :math:`\mu`\ m\ :sup:`2` cm\ :sup:`-3`.
    - Volume density is in :math:`\mu`\ m\ :sup:`3` cm\ :sup:`-3`.

The implementations use number distributions per unit radius. For spherical
particles, the surface-area and volume densities are respectively
:math:`4\pi M_2` and :math:`4\pi M_3/3`, where :math:`M_i` is raw moment ``i``.

Formulae follow R. G. Grainger, *Some Useful Formulae for Particle Size
Distributions and their Optical Properties*, revised 27 August 2026.
"""

from abc import ABC, abstractmethod

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaincinv, gammaln, logsumexp, ndtr, ndtri


class SizeDistribution(ABC):
    """Abstract base class for particle-radius number distributions.

    Concrete distributions provide at least an arithmetic :meth:`mean`. The
    implementations in this module also provide ``value(radii)`` and store a
    short machine-readable distribution name in :attr:`type`.
    """

    def __init__(self, type):
        """Initialize a size distribution with its distribution type.

        Args:
            type: Machine-readable name of the particle-size distribution.
        """
        self.type = type

    @abstractmethod
    def mean(self):
        """Return the arithmetic mean particle radius in micrometres.

        Returns:
            Arithmetic mean radius. Heavy-tailed distributions return
            ``np.inf`` when their first moment does not exist.
        """
        pass


def _positive_parameter(name, value):
    """Validate and convert a positive scalar distribution parameter.

    Args:
        name: Parameter name used in the error message.
        value: Candidate scalar value.

    Returns:
        The value converted to ``float``.

    Raises:
        ValueError: If ``value`` is missing, non-scalar, non-finite, or not
            strictly positive.
    """
    if value is None or not np.isscalar(value) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite number greater than 0.")
    return float(value)


def _finite_parameter(name, value):
    """Validate and convert a finite scalar distribution parameter.

    Args:
        name: Parameter name used in the error message.
        value: Candidate scalar value.

    Returns:
        The value converted to ``float``.

    Raises:
        ValueError: If ``value`` is missing, non-scalar, or non-finite.
    """
    if value is None or not np.isscalar(value) or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite number.")
    return float(value)


def _exp(log_value):
    """Exponentiate a scalar stored in the log domain.

    Overflow and underflow are permitted because a mathematically valid moment
    can be outside the floating-point range. The returned value is then
    ``np.inf`` or ``0.0``. Invalid operations are deliberately not suppressed.

    Args:
        log_value: Natural logarithm of the desired value.

    Returns:
        ``exp(log_value)`` converted to ``float``.
    """
    with np.errstate(over="ignore", under="ignore"):
        return float(np.exp(log_value))


def _return_scalar_when_scalar(original, values):
    """Preserve scalar input/output behaviour for vectorized evaluations.

    Args:
        original: Original radius argument supplied by the caller.
        values: NumPy result produced from that argument.

    Returns:
        A Python ``float`` when ``original`` was scalar; otherwise ``values``.
    """
    if np.ndim(original) == 0:
        return float(np.asarray(values))
    return values


def _set_concentration(
    distribution,
    n,
    surface_area_density,
    volume_density,
    second_moment_per_particle,
    third_moment_per_particle,
):
    r"""Set the three equivalent concentration measures on a distribution.

    Inputs use the same deterministic precedence as
    :class:`LogNormalDistribution`: ``n`` first, then surface-area density,
    then volume density. For normalized per-particle raw moments
    :math:`\langle r^i\rangle`, spherical-particle densities are

    .. math::

       A_0 = 4\pi N_0\langle r^2\rangle, \qquad
       V_0 = \frac{4\pi}{3}N_0\langle r^3\rangle.

    Args:
        distribution: Object on which the three densities are set.
        n: Total number concentration in cm\ :sup:`-3`.
        surface_area_density: Total spherical surface area in
            :math:`\mu`\ m\ :sup:`2` cm\ :sup:`-3`.
        volume_density: Total spherical volume in
            :math:`\mu`\ m\ :sup:`3` cm\ :sup:`-3`.
        second_moment_per_particle: Normalized raw second moment
            :math:`\langle r^2\rangle`.
        third_moment_per_particle: Normalized raw third moment
            :math:`\langle r^3\rangle`.

    Raises:
        ValueError: If no concentration measure is supplied, a supplied value
            is non-positive, or the required moment is divergent.
    """
    if n is not None:
        distribution.n = _positive_parameter("n", n)
    elif surface_area_density is not None:
        surface_area_density = _positive_parameter(
            "surface_area_density", surface_area_density
        )
        if not np.isfinite(second_moment_per_particle):
            raise ValueError(
                "surface_area_density cannot parameterise a distribution whose "
                "second moment is not finite."
            )
        distribution.n = surface_area_density / (
            4.0 * np.pi * second_moment_per_particle
        )
    elif volume_density is not None:
        volume_density = _positive_parameter("volume_density", volume_density)
        if not np.isfinite(third_moment_per_particle):
            raise ValueError(
                "volume_density cannot parameterise a distribution whose third "
                "moment is not finite."
            )
        distribution.n = (
            3.0 * volume_density / (4.0 * np.pi * third_moment_per_particle)
        )
    else:
        raise ValueError("Provide one of n, surface_area_density, or volume_density.")

    distribution.surface_area_density = (
        4.0 * np.pi * distribution.n * second_moment_per_particle
    )
    distribution.volume_density = (
        4.0 * np.pi * distribution.n * third_moment_per_particle / 3.0
    )


def _modified_gamma_effective_variance(alpha, gamma):
    r"""Evaluate effective variance for modified-gamma shape parameters.

    This implements

    .. math::

       v_e = \frac{\Gamma(s_2)\Gamma(s_4)}{\Gamma(s_3)^2} - 1,
       \qquad s_i = \frac{\alpha+i+1}{\gamma}.

    Args:
        alpha: Small-radius power exponent, greater than ``-1`` for a
            normalizable distribution.
        gamma: Positive large-radius cutoff exponent.

    Returns:
        Dimensionless effective variance.
    """
    s2 = (alpha + 3.0) / gamma
    s3 = (alpha + 4.0) / gamma
    s4 = (alpha + 5.0) / gamma
    log_ratio = gammaln(s2) + gammaln(s4) - 2.0 * gammaln(s3)
    with np.errstate(over="ignore"):
        return float(np.expm1(log_ratio))


def _modified_gamma_alpha(effective_variance, gamma):
    """Recover modified-gamma ``alpha`` from effective variance.

    Grainger Equation 138 has no general closed-form inverse. At fixed
    ``gamma`` this function brackets the normalizable interval ``alpha > -1``
    and solves it with Brent's method. The ordinary-gamma case ``gamma == 1``
    uses the exact result ``alpha = 1 / effective_variance - 3``.

    Args:
        effective_variance: Positive dimensionless effective variance.
        gamma: Positive modified-gamma cutoff exponent.

    Returns:
        The corresponding exponent ``alpha``.

    Raises:
        ValueError: If the parameters are outside the normalizable range or a
            numerical bracket cannot be found.
    """
    effective_variance = _positive_parameter("effective_variance", effective_variance)
    gamma = _positive_parameter("gamma", gamma)

    if gamma == 1.0:
        alpha = 1.0 / effective_variance - 3.0
        if alpha <= -1.0:
            raise ValueError(
                "effective_variance must be less than 0.5 for a normalisable "
                "gamma distribution."
            )
        return alpha

    maximum_variance = _modified_gamma_effective_variance(-1.0, gamma)
    if effective_variance >= maximum_variance:
        raise ValueError(
            "effective_variance is outside the normalisable modified-gamma "
            f"range for gamma={gamma}; it must be less than {maximum_variance}."
        )

    def residual(alpha):
        """Return the variance mismatch used by Brent's root solver."""
        return _modified_gamma_effective_variance(alpha, gamma) - effective_variance

    lower = np.nextafter(-1.0, 0.0)
    upper = max(1.0, 2.0 / (gamma * effective_variance))
    while residual(upper) > 0.0:
        upper *= 2.0
        if not np.isfinite(upper):
            raise ValueError(
                "Could not bracket alpha for the requested effective variance."
            )

    return float(brentq(residual, lower, upper, xtol=1e-12, rtol=1e-12))


def _positive_vector(name, values):
    """Validate a non-empty vector of positive finite parameters.

    Args:
        name: Parameter name used in error messages.
        values: Scalar or one-dimensional array-like values.

    Returns:
        One-dimensional ``float`` NumPy array. A scalar becomes length one.

    Raises:
        ValueError: If the result is empty, multidimensional, non-finite, or
            contains a non-positive value.
    """
    if values is None:
        raise ValueError(f"{name} must be provided.")
    result = np.atleast_1d(np.asarray(values, dtype=float))
    if result.ndim != 1 or result.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional sequence.")
    if np.any(~np.isfinite(result)) or np.any(result <= 0.0):
        raise ValueError(f"Every value in {name} must be finite and greater than 0.")
    return result


class GaussianDistribution(SizeDistribution):
    r"""Create a Gaussian particle-radius distribution.

    The untruncated number distribution is

    .. math::

       n(r) = \frac{N_0}{\sqrt{2\pi}\sigma_0}
       \exp\left[-\frac{(r-\mu_0)^2}{2\sigma_0^2}\right].

    Args:
        n: Total number concentration in cm\ :sup:`-3`.
        r: Mean ``mu_0`` of the underlying normal distribution in micrometres.
        s: Standard deviation ``sigma_0`` in micrometres.
        surface_area_density: Spherical surface-area density in
            :math:`\mu`\ m\ :sup:`2` cm\ :sup:`-3`, as an alternative to ``n``.
        volume_density: Spherical volume density in
            :math:`\mu`\ m\ :sup:`3` cm\ :sup:`-3`, as an alternative to ``n``.
        truncate: If true, restrict the distribution to ``r > 0`` and
            renormalize it. The default is false, matching the conventional
            mathematical Gaussian and preserving ``mean() == r``.

    Attributes:
        median_radius: Median of the represented distribution.
        mode_radius: Radius of maximum density.
        effective_radius: Area-weighted mean radius ``M3 / M2``.
        effective_variance: Dimensionless area-weighted variance
            ``M2 * M4 / M3**2 - 1``.

    Notes:
        An untruncated Gaussian assigns finite probability to negative radii.
        Its ``n`` and raw moments are defined over the entire real line, so its
        third moment and volume density are mathematical rather than physical
        particle quantities. Use ``truncate=True`` when ``n`` must describe
        only positive particle radii.

    Raises:
        ValueError: If a required scalar or concentration is not positive and
            finite, or no concentration measure is supplied.
        TypeError: If ``truncate`` is not boolean.
    """

    def __init__(
        self,
        n=None,
        r=None,
        s=None,
        surface_area_density=None,
        volume_density=None,
        *,
        truncate=False,
    ):
        super().__init__("gaussian")
        self.r = _positive_parameter("r", r)
        self.s = _positive_parameter("s", s)
        if not isinstance(truncate, (bool, np.bool_)):
            raise TypeError("truncate must be a boolean.")
        self.truncate = bool(truncate)

        mu = self.r
        sigma = self.s
        if self.truncate:
            positive_fraction = float(ndtr(mu / sigma))
            boundary_density = np.exp(-0.5 * (mu / sigma) ** 2) / np.sqrt(2.0 * np.pi)
            unnormalised = np.empty(5, dtype=float)
            unnormalised[0] = positive_fraction
            unnormalised[1] = mu * positive_fraction + sigma * boundary_density
            for order in range(2, 5):
                unnormalised[order] = (
                    mu * unnormalised[order - 1]
                    + (order - 1) * sigma**2 * unnormalised[order - 2]
                )
            self._moments_per_particle = unnormalised / positive_fraction
            self._positive_fraction = positive_fraction
            target_probability = 1.0 - 0.5 * positive_fraction
            self.median_radius = mu + sigma * ndtri(target_probability)
        else:
            self._moments_per_particle = np.array(
                [
                    1.0,
                    mu,
                    mu**2 + sigma**2,
                    mu**3 + 3.0 * mu * sigma**2,
                    mu**4 + 6.0 * mu**2 * sigma**2 + 3.0 * sigma**4,
                ]
            )
            self._positive_fraction = 1.0
            self.median_radius = mu

        _set_concentration(
            self,
            n,
            surface_area_density,
            volume_density,
            self._moments_per_particle[2],
            self._moments_per_particle[3],
        )
        self.effective_radius = (
            self._moments_per_particle[3] / self._moments_per_particle[2]
        )
        self.effective_variance = (
            self._moments_per_particle[2]
            * self._moments_per_particle[4]
            / self._moments_per_particle[3] ** 2
            - 1.0
        )
        self.mode_radius = mu
        self._log_prefactor = (
            np.log(self.n)
            - np.log(sigma)
            - 0.5 * np.log(2.0 * np.pi)
            - np.log(self._positive_fraction)
        )

    def moment(self, order):
        """Return raw moment ``order`` for orders zero through four.

        Args:
            order: Integer moment order in the inclusive interval ``[0, 4]``.

        Returns:
            Raw moment including the total number concentration.

        Raises:
            ValueError: If ``order`` is not an integer from zero through four.
        """
        if not isinstance(order, (int, np.integer)) or not 0 <= order <= 4:
            raise ValueError("order must be an integer from 0 to 4.")
        return self.n * self._moments_per_particle[order]

    def mean(self):
        """Return the arithmetic mean radius.

        Returns:
            ``r`` for the untruncated Gaussian, or the conditional mean over
            positive radii for a truncated Gaussian.
        """
        return self._moments_per_particle[1]

    def value(self, radii):
        """Evaluate differential number density at one or more radii.

        Args:
            radii: Scalar or array-like particle radii in micrometres.

        Returns:
            Number density per micrometre. Scalar input returns ``float`` and
            array-like input returns a NumPy array of the same shape. For a
            truncated distribution, values at non-positive radii are zero.
        """
        radius = np.asarray(radii, dtype=float)
        log_value = self._log_prefactor - 0.5 * ((radius - self.r) / self.s) ** 2
        with np.errstate(over="ignore", under="ignore"):
            values = np.exp(log_value)
        if self.truncate:
            values = np.where(radius > 0.0, values, 0.0)
        return _return_scalar_when_scalar(radii, values)


class LogNormalDistribution(SizeDistribution):
    r"""Create a single-mode log-normal particle-radius distribution.

    The normalized number distribution is

    .. math::

       n(r) = \frac{N_0}{\sqrt{2\pi}\ln(S)r}
       \exp\left[-\frac{(\ln r-\ln r_m)^2}{2\ln^2(S)}\right],

    where ``r`` is the number-median radius :math:`r_m` and ``s`` is the
    geometric standard deviation :math:`S`.

    Args:
        n: Total number concentration in cm\ :sup:`-3`.
        r: Number-median radius in micrometres.
        s: Geometric standard deviation. It must exceed one; ``s == 1``
            represents a delta distribution rather than a finite density.
        surface_area_density: Spherical surface-area density in
            :math:`\mu`\ m\ :sup:`2` cm\ :sup:`-3`, as an alternative to ``n``.
        volume_density: Spherical volume density in
            :math:`\mu`\ m\ :sup:`3` cm\ :sup:`-3`, as an alternative to ``n``.

    Attributes:
        median_radius: Number-median radius, equal to ``r``.
        mode_radius: Radius at which ``dN/dr`` is maximal.
        effective_radius: Area-weighted mean radius ``M3 / M2``.
        effective_variance: Dimensionless area-weighted variance
            ``M2 * M4 / M3**2 - 1``.

    Raises:
        ValueError: If ``r`` is not positive, ``s`` does not exceed one, no
            concentration measure is supplied, or a concentration is invalid.
    """

    def __init__(
        self, n=None, r=None, s=None, surface_area_density=None, volume_density=None
    ):
        super().__init__("log_normal")
        self.r = _positive_parameter("r", r)
        self.s = _positive_parameter("s", s)
        if self.s <= 1.0:
            raise ValueError("s must be greater than 1.")
        self.lnr = np.log(self.r)
        self.lns = np.log(self.s)

        second_moment = self._moment_per_particle(2)
        third_moment = self._moment_per_particle(3)
        _set_concentration(
            self,
            n,
            surface_area_density,
            volume_density,
            second_moment,
            third_moment,
        )
        self.median_radius = self.r
        self.mode_radius = self.r * np.exp(-(self.lns**2))
        self.effective_radius = third_moment / second_moment
        fourth_moment = self._moment_per_particle(4)
        self.effective_variance = second_moment * fourth_moment / third_moment**2 - 1.0
        self._log_prefactor = (
            np.log(self.n) - 0.5 * np.log(2.0 * np.pi) - np.log(self.lns)
        )

    def _moment_per_particle(self, order):
        """Return normalized raw moment ``order`` without number concentration."""
        return _exp(order * self.lnr + 0.5 * order**2 * self.lns**2)

    def moment(self, order):
        r"""Return a non-negative integer raw moment.

        Args:
            order: Non-negative integer moment order.

        Returns:
            :math:`M_i = N_0 r_m^i\exp(i^2\ln^2(S)/2)`.

        Raises:
            ValueError: If ``order`` is not a non-negative integer.
        """
        if not isinstance(order, (int, np.integer)) or order < 0:
            raise ValueError("order must be a non-negative integer.")
        return self.n * self._moment_per_particle(order)

    def mean(self):
        """Return the arithmetic mean radius in micrometres.

        Returns:
            ``r * exp(log(s)**2 / 2)``.
        """
        return self._moment_per_particle(1)

    def value(self, radii):
        """Evaluate differential number density at positive radii.

        Args:
            radii: Scalar or array-like particle radii in micrometres.

        Returns:
            Number density per micrometre. Scalar input returns ``float`` and
            array-like input returns a NumPy array of the same shape. Values at
            non-positive radii are zero.
        """
        radius = np.asarray(radii, dtype=float)
        values = np.zeros_like(radius)
        positive = (radius > 0.0) & np.isfinite(radius)
        log_radius = np.log(radius[positive])
        log_value = (
            self._log_prefactor
            - log_radius
            - 0.5 * ((log_radius - self.lnr) / self.lns) ** 2
        )
        with np.errstate(over="ignore", under="ignore"):
            values[positive] = np.exp(log_value)
        return _return_scalar_when_scalar(radii, values)


class MultimodeLogNormalDistribution(SizeDistribution):
    r"""Create a sum of independent log-normal number-distribution modes.

    If mode ``i`` has number concentration :math:`N_i`, number-median radius
    :math:`r_i`, and geometric standard deviation :math:`S_i`, then

    .. math::

       n(r) = \sum_i \frac{N_i}{\sqrt{2\pi}\ln(S_i)r}
       \exp\left[-\frac{(\ln r-\ln r_i)^2}{2\ln^2(S_i)}\right].

    Raw moments add linearly across modes. Consequently, mixture effective
    radius is ``M3 / M2`` and is generally not the number-weighted average of
    the component effective radii.

    Args:
        n: One-dimensional sequence of per-mode number concentrations in
            cm\ :sup:`-3`.
        r: One-dimensional sequence of per-mode number-median radii in
            micrometres.
        s: One-dimensional sequence of per-mode geometric standard deviations.
            Every value must exceed one.
        surface_area_density: Sequence of per-mode spherical surface-area
            densities, as an alternative to ``n``.
        volume_density: Sequence of per-mode spherical volume densities, as an
            alternative to ``n``.

    Attributes:
        number_of_modes: Number of component modes.
        mode_number_densities: Number concentration of every component.
        mode_surface_area_densities: Surface-area density of every component.
        mode_volume_densities: Volume density of every component.
        component_mode_radii: Locations of the maxima of each ``dN/dr`` mode.
        median_radius: Number median of the complete mixture, found numerically.
        effective_radius: Complete-mixture ``M3 / M2``.
        effective_variance: Complete-mixture ``M2 * M4 / M3**2 - 1``.

    Notes:
        Density alternatives must be supplied per mode. A single mixture total
        cannot determine how concentration should be apportioned among modes.
        If more than one concentration representation is supplied, precedence
        is ``n``, then ``surface_area_density``, then ``volume_density``.

    Raises:
        ValueError: If arrays are empty, non-positive, non-finite, have unequal
            lengths, or if no concentration representation is supplied.
    """

    def __init__(
        self,
        n=None,
        r=None,
        s=None,
        surface_area_density=None,
        volume_density=None,
    ):
        super().__init__("multimode_log_normal")
        self.r = _positive_vector("r", r)
        self.s = _positive_vector("s", s)
        if self.r.shape != self.s.shape:
            raise ValueError("r and s must contain the same number of modes.")
        if np.any(self.s <= 1.0):
            raise ValueError("Every geometric standard deviation in s must exceed 1.")

        self.number_of_modes = self.r.size
        self._log_r = np.log(self.r)
        self._log_s = np.log(self.s)
        area_per_particle = 4.0 * np.pi * self.r**2 * np.exp(2.0 * self._log_s**2)
        volume_per_particle = (
            4.0 * np.pi * self.r**3 * np.exp(4.5 * self._log_s**2) / 3.0
        )

        if n is not None:
            mode_number_densities = _positive_vector("n", n)
        elif surface_area_density is not None:
            mode_surface_area_densities = _positive_vector(
                "surface_area_density", surface_area_density
            )
            self._check_mode_count("surface_area_density", mode_surface_area_densities)
            mode_number_densities = mode_surface_area_densities / area_per_particle
        elif volume_density is not None:
            mode_volume_densities = _positive_vector("volume_density", volume_density)
            self._check_mode_count("volume_density", mode_volume_densities)
            mode_number_densities = mode_volume_densities / volume_per_particle
        else:
            raise ValueError(
                "Provide per-mode n, surface_area_density, or volume_density."
            )
        self._check_mode_count("n", mode_number_densities)

        self.mode_number_densities = mode_number_densities.copy()
        self.mode_surface_area_densities = (
            self.mode_number_densities * area_per_particle
        )
        self.mode_volume_densities = self.mode_number_densities * volume_per_particle
        self.n = float(np.sum(self.mode_number_densities))
        self.surface_area_density = float(np.sum(self.mode_surface_area_densities))
        self.volume_density = float(np.sum(self.mode_volume_densities))
        self._log_mode_number_densities = np.log(self.mode_number_densities)
        self._log_prefactors = (
            self._log_mode_number_densities
            - 0.5 * np.log(2.0 * np.pi)
            - np.log(self._log_s)
        )

        second_moment = self.moment(2)
        third_moment = self.moment(3)
        fourth_moment = self.moment(4)
        self.effective_radius = third_moment / second_moment
        self.effective_variance = second_moment * fourth_moment / third_moment**2 - 1.0
        self.component_mode_radii = self.r * np.exp(-(self._log_s**2))
        self.median_radius = self._mixture_median()

    def _check_mode_count(self, name, values):
        """Require an input vector to contain one value for every mode.

        Args:
            name: Input name used in a possible error message.
            values: One-dimensional NumPy array to compare with ``r``.

        Raises:
            ValueError: If the array length differs from the number of modes.
        """
        if values.shape != self.r.shape:
            raise ValueError(
                f"{name} must contain one value for each of the "
                f"{self.number_of_modes} modes."
            )

    def _mixture_median(self):
        """Calculate the number median of the complete mixture.

        Returns:
            Positive radius whose cumulative number concentration is ``n / 2``.
        """

        def centred_cdf(log_radius):
            """Return normalized mixture CDF minus one half in log-radius."""
            standardised = (log_radius - self._log_r) / self._log_s
            return float(
                np.dot(self.mode_number_densities, ndtr(standardised)) / self.n - 0.5
            )

        lower = float(np.min(self._log_r - 10.0 * self._log_s))
        upper = float(np.max(self._log_r + 10.0 * self._log_s))
        return _exp(brentq(centred_cdf, lower, upper, xtol=1e-12, rtol=1e-12))

    def moment(self, order):
        """Return a non-negative integer raw moment of the mixture.

        Args:
            order: Non-negative integer moment order.

        Returns:
            Sum of component moments
            ``N[i] * r[i]**order * exp(order**2 * log(s[i])**2 / 2)``.

        Raises:
            ValueError: If ``order`` is not a non-negative integer.
        """
        if not isinstance(order, (int, np.integer)) or order < 0:
            raise ValueError("order must be a non-negative integer.")
        log_mode_moments = (
            self._log_mode_number_densities
            + order * self._log_r
            + 0.5 * order**2 * self._log_s**2
        )
        return _exp(logsumexp(log_mode_moments))

    def mean(self):
        """Return the number-weighted arithmetic mean radius.

        Returns:
            Complete-mixture first raw moment divided by total concentration.
        """
        return self.moment(1) / self.n

    def cdf(self, radii):
        r"""Evaluate cumulative number concentration below given radii.

        Args:
            radii: Scalar or array-like particle radii in micrometres.

        Returns:
            Cumulative number concentration in cm\ :sup:`-3`. Scalar input
            returns ``float``; array-like input returns a matching NumPy array.
        """
        radius = np.asarray(radii, dtype=float)
        values = np.zeros_like(radius)
        positive = radius > 0.0
        log_radius = np.log(radius[positive])
        standardised = (
            log_radius[np.newaxis, :] - self._log_r[:, np.newaxis]
        ) / self._log_s[:, np.newaxis]
        values[positive] = np.sum(
            self.mode_number_densities[:, np.newaxis] * ndtr(standardised), axis=0
        )
        return _return_scalar_when_scalar(radii, values)

    def value(self, radii):
        """Evaluate the summed differential number distribution.

        Args:
            radii: Scalar or array-like particle radii in micrometres.

        Returns:
            Number density per micrometre. Scalar input returns ``float`` and
            array-like input returns a matching NumPy array. Values at
            non-positive radii are zero.
        """
        radius = np.asarray(radii, dtype=float)
        values = np.zeros_like(radius)
        positive = (radius > 0.0) & np.isfinite(radius)
        log_radius = np.log(radius[positive])
        standardised = (
            log_radius[np.newaxis, :] - self._log_r[:, np.newaxis]
        ) / self._log_s[:, np.newaxis]
        log_components = (
            self._log_prefactors[:, np.newaxis]
            - log_radius[np.newaxis, :]
            - 0.5 * standardised**2
        )
        values[positive] = np.exp(logsumexp(log_components, axis=0))
        return _return_scalar_when_scalar(radii, values)


class ModifiedGammaDistribution(SizeDistribution):
    r"""Create a normalized modified-gamma particle-radius distribution.

    The distribution is

    .. math::

       n(r) = N_0\frac{\gamma b^{(\alpha+1)/\gamma}}
       {\Gamma((\alpha+1)/\gamma)}r^\alpha\exp(-br^\gamma).

    The caller supplies effective variance and the positive cutoff exponent
    ``gamma``. The small-radius exponent ``alpha`` is obtained numerically from
    effective variance. A radius scale is then supplied either as the true
    number median ``r`` or as ``effective_radius``.

    Args:
        n: Total number concentration in cm\ :sup:`-3`.
        r: Number-median radius in micrometres. Mutually exclusive with
            ``effective_radius``.
        s: Positive dimensionless effective variance.
        surface_area_density: Spherical surface-area density, as an alternative
            to ``n``.
        volume_density: Spherical volume density, as an alternative to ``n``.
        gamma: Positive exponent controlling the large-radius cutoff.
        effective_radius: Area-weighted mean radius in micrometres, as an
            alternative to the median radius ``r``.

    Attributes:
        alpha: Solved small-radius power exponent. Normalization requires
            ``alpha > -1``.
        b: Positive scale coefficient with units micrometres\ :sup:`-gamma`.
        median_radius: True number median, not the modal radius.
        mode_radius: Number-distribution mode, or zero when ``alpha <= 0``.
        effective_radius: Area-weighted mean radius ``M3 / M2``.
        effective_variance: Supplied dimensionless effective variance.

    Notes:
        Grainger denotes the *mode* by :math:`r_m` in the modified-gamma
        mode-radius formula. In this API, ``r`` follows
        :class:`LogNormalDistribution` and always denotes a number median.

    Raises:
        ValueError: If parameters are non-positive or non-finite, both/neither
            radius representations are supplied, or effective variance lies
            outside the normalizable range for the selected ``gamma``.
    """

    def __init__(
        self,
        n=None,
        r=None,
        s=None,
        surface_area_density=None,
        volume_density=None,
        *,
        gamma=None,
        effective_radius=None,
    ):
        super().__init__("modified_gamma")
        self.gamma = _positive_parameter("gamma", gamma)
        self.s = _positive_parameter("s", s)
        self.effective_variance = self.s
        self.alpha = _modified_gamma_alpha(self.effective_variance, self.gamma)

        if (r is None) == (effective_radius is None):
            raise ValueError("Provide exactly one of r or effective_radius.")

        s0 = (self.alpha + 1.0) / self.gamma
        s2 = (self.alpha + 3.0) / self.gamma
        s3 = (self.alpha + 4.0) / self.gamma
        median_gamma_variable = float(gammaincinv(s0, 0.5))

        if r is not None:
            self.r = _positive_parameter("r", r)
            self.median_radius = self.r
            self._log_b = np.log(median_gamma_variable) - self.gamma * np.log(
                self.median_radius
            )
            self.effective_radius = _exp(
                -self._log_b / self.gamma + gammaln(s3) - gammaln(s2)
            )
        else:
            self.effective_radius = _positive_parameter(
                "effective_radius", effective_radius
            )
            self._log_b = self.gamma * (
                gammaln(s3) - gammaln(s2) - np.log(self.effective_radius)
            )
            self.median_radius = _exp(
                (np.log(median_gamma_variable) - self._log_b) / self.gamma
            )
            self.r = self.median_radius

        self.b = _exp(self._log_b)
        if not np.isfinite(self.b) or self.b <= 0.0:
            raise ValueError("The requested parameters produce a non-finite b.")

        second_moment = self._moment_per_particle(2)
        third_moment = self._moment_per_particle(3)
        _set_concentration(
            self,
            n,
            surface_area_density,
            volume_density,
            second_moment,
            third_moment,
        )
        self.mode_radius = (
            _exp((np.log(self.alpha) - np.log(self.gamma) - self._log_b) / self.gamma)
            if self.alpha > 0.0
            else 0.0
        )
        self._log_prefactor = (
            np.log(self.n) + np.log(self.gamma) + s0 * self._log_b - gammaln(s0)
        )

    def _moment_per_particle(self, order):
        """Return normalized raw moment ``order`` without concentration.

        Args:
            order: Non-negative raw-moment order.

        Returns:
            ``b**(-order/gamma) * Gamma(s_order) / Gamma(s_0)``.
        """
        si = (self.alpha + order + 1.0) / self.gamma
        s0 = (self.alpha + 1.0) / self.gamma
        return _exp(-order * self._log_b / self.gamma + gammaln(si) - gammaln(s0))

    def moment(self, order):
        """Return a non-negative integer raw moment.

        Args:
            order: Non-negative integer moment order.

        Returns:
            Raw moment including total number concentration.

        Raises:
            ValueError: If ``order`` is not a non-negative integer.
        """
        if not isinstance(order, (int, np.integer)) or order < 0:
            raise ValueError("order must be a non-negative integer.")
        return self.n * self._moment_per_particle(order)

    def mean(self):
        """Return the finite arithmetic mean radius in micrometres."""
        return self._moment_per_particle(1)

    def value(self, radii):
        """Evaluate differential number density at positive radii.

        Args:
            radii: Scalar or array-like particle radii in micrometres.

        Returns:
            Number density per micrometre. Scalar input returns ``float`` and
            array-like input returns a matching NumPy array. Values at
            non-positive or non-finite radii are zero.
        """
        radius = np.asarray(radii, dtype=float)
        values = np.zeros_like(radius)
        positive = (radius > 0.0) & np.isfinite(radius)
        log_radius = np.log(radius[positive])
        with np.errstate(over="ignore", under="ignore"):
            cutoff = np.exp(self._log_b + self.gamma * log_radius)
            log_value = self._log_prefactor + self.alpha * log_radius - cutoff
            values[positive] = np.exp(log_value)
        return _return_scalar_when_scalar(radii, values)


class GammaDistribution(ModifiedGammaDistribution):
    r"""Create an ordinary gamma distribution in optical-size parameters.

    This is the ``gamma = 1`` special case of
    :class:`ModifiedGammaDistribution`. In terms of effective radius ``r_e``
    and effective variance ``v_e`` its parameters are

    .. math::

       \alpha = \frac{1}{v_e}-3, \qquad b = \frac{1}{r_e v_e}.

    Args:
        n: Total number concentration in cm\ :sup:`-3`.
        r: Number-median radius in micrometres. Mutually exclusive with
            ``effective_radius``.
        s: Effective variance ``v_e``. Normalization requires ``0 < s < 0.5``.
        surface_area_density: Spherical surface-area density, as an alternative
            to ``n``.
        volume_density: Spherical volume density, as an alternative to ``n``.
        effective_radius: Area-weighted mean radius in micrometres, as an
            alternative to the number median ``r``.

    Notes:
        An interior number-distribution mode requires ``s < 1/3``. For
        ``1/3 <= s < 1/2`` the distribution is still normalizable but is
        maximal at zero radius.
    """

    def __init__(
        self,
        n=None,
        r=None,
        s=None,
        surface_area_density=None,
        volume_density=None,
        *,
        effective_radius=None,
    ):
        super().__init__(
            n=n,
            r=r,
            s=s,
            surface_area_density=surface_area_density,
            volume_density=volume_density,
            gamma=1.0,
            effective_radius=effective_radius,
        )
        self.type = "gamma"


class InverseModifiedGammaDistribution(SizeDistribution):
    r"""Create a normalized inverse modified-gamma distribution.

    The distribution is

    .. math::

       n(r) = N_0\frac{\gamma b^{(\alpha-1)/\gamma}}
       {\Gamma((\alpha-1)/\gamma)}r^{-\alpha}\exp(-br^{-\gamma}).

    It is exponentially suppressed at small radius and has a power-law tail at
    large radius. Moment ``i`` exists only when ``i < alpha - 1``.

    Args:
        n: Total number concentration in cm\ :sup:`-3`.
        alpha: Large-radius power exponent. Normalization requires ``alpha > 1``.
        b: Positive scale coefficient with units micrometres\ :sup:`gamma`.
            Mutually exclusive with ``median_radius``.
        gamma: Positive inverse-power cutoff exponent.
        surface_area_density: Spherical surface-area density, as an alternative
            to ``n``. This requires ``alpha > 3``.
        volume_density: Spherical volume density, as an alternative to ``n``.
            This requires ``alpha > 4``.
        median_radius: Number-median radius in micrometres, as an alternative
            to ``b``.

    Attributes:
        median_radius: True number median.
        mode_radius: Number-distribution mode.
        effective_radius: ``M3 / M2`` when ``alpha > 4``, otherwise ``np.inf``.
        effective_variance: ``M2 * M4 / M3**2 - 1`` when ``alpha > 5``,
            otherwise ``np.inf``.

    Notes:
        A distribution can have finite number concentration while its mean,
        area, volume, or optical metrics diverge. Undefined heavy-tail moments
        are represented by ``np.inf`` rather than silently truncated.

    Raises:
        ValueError: If shape or scale parameters are outside their domains,
            both/neither scale representations are supplied, or a density input
            requires a divergent moment.
    """

    def __init__(
        self,
        n=None,
        alpha=None,
        b=None,
        gamma=None,
        surface_area_density=None,
        volume_density=None,
        *,
        median_radius=None,
    ):
        super().__init__("inverse_modified_gamma")
        self.alpha = _finite_parameter("alpha", alpha)
        if self.alpha <= 1.0:
            raise ValueError("alpha must be greater than 1 for normalisation.")
        self.gamma = _positive_parameter("gamma", gamma)
        if (b is None) == (median_radius is None):
            raise ValueError("Provide exactly one of b or median_radius.")

        t0 = (self.alpha - 1.0) / self.gamma
        median_gamma_variable = float(gammaincinv(t0, 0.5))
        if median_radius is not None:
            self.median_radius = _positive_parameter("median_radius", median_radius)
            self._log_b = np.log(median_gamma_variable) + self.gamma * np.log(
                self.median_radius
            )
            self.b = _exp(self._log_b)
        else:
            self.b = _positive_parameter("b", b)
            self._log_b = np.log(self.b)
            self.median_radius = _exp(
                (self._log_b - np.log(median_gamma_variable)) / self.gamma
            )
        if not np.isfinite(self.b) or self.b <= 0.0:
            raise ValueError("The requested parameters produce a non-finite b.")

        self.r = self.median_radius
        self.s = None
        second_moment = self._moment_per_particle(2)
        third_moment = self._moment_per_particle(3)
        _set_concentration(
            self,
            n,
            surface_area_density,
            volume_density,
            second_moment,
            third_moment,
        )
        self.mode_radius = _exp(
            (self._log_b + np.log(self.gamma) - np.log(self.alpha)) / self.gamma
        )
        self.effective_radius = (
            third_moment / second_moment if np.isfinite(third_moment) else np.inf
        )
        fourth_moment = self._moment_per_particle(4)
        self.effective_variance = (
            second_moment * fourth_moment / third_moment**2 - 1.0
            if np.isfinite(fourth_moment)
            else np.inf
        )
        self._log_prefactor = (
            np.log(self.n) + np.log(self.gamma) + t0 * self._log_b - gammaln(t0)
        )

    def _moment_per_particle(self, order):
        """Return a normalized raw moment or ``np.inf`` when divergent.

        Args:
            order: Non-negative raw-moment order.

        Returns:
            ``b**(order/gamma) * Gamma(t_order) / Gamma(t_0)`` when
            ``order < alpha - 1``; otherwise ``np.inf``.
        """
        if order >= self.alpha - 1.0:
            return np.inf
        t0 = (self.alpha - 1.0) / self.gamma
        ti = (self.alpha - 1.0 - order) / self.gamma
        return _exp(order * self._log_b / self.gamma + gammaln(ti) - gammaln(t0))

    def moment(self, order):
        """Return a raw moment, or ``np.inf`` when it diverges.

        Args:
            order: Non-negative integer moment order.

        Returns:
            Raw moment including total number concentration.

        Raises:
            ValueError: If ``order`` is not a non-negative integer.
        """
        if not isinstance(order, (int, np.integer)) or order < 0:
            raise ValueError("order must be a non-negative integer.")
        return self.n * self._moment_per_particle(order)

    def mean(self):
        """Return the arithmetic mean, or ``np.inf`` if ``alpha <= 2``."""
        return self._moment_per_particle(1)

    def value(self, radii):
        """Evaluate differential number density at positive radii.

        Args:
            radii: Scalar or array-like particle radii in micrometres.

        Returns:
            Number density per micrometre. Scalar input returns ``float`` and
            array-like input returns a matching NumPy array. Values at
            non-positive radii are zero.
        """
        radius = np.asarray(radii, dtype=float)
        values = np.zeros_like(radius)
        positive = radius > 0.0
        log_radius = np.log(radius[positive])
        with np.errstate(over="ignore", under="ignore"):
            inverse_cutoff = np.exp(self._log_b - self.gamma * log_radius)
            log_value = self._log_prefactor - self.alpha * log_radius - inverse_cutoff
            values[positive] = np.exp(log_value)
        return _return_scalar_when_scalar(radii, values)


class RegularisedPowerLawDistribution(SizeDistribution):
    r"""Create a normalized regularised power-law distribution.

    The distribution is

    .. math::

       n(r) = N_0\alpha(\gamma-1)b^{-\alpha}r^{\alpha-1}
       [1+(r/b)^\alpha]^{-\gamma}.

    ``alpha`` controls the small-radius slope and transition sharpness, ``b``
    is a radius scale, and ``gamma`` controls the large-radius tail. Moment
    ``i`` exists only when ``gamma > 1 + i / alpha``.

    Args:
        n: Total number concentration in cm\ :sup:`-3`.
        alpha: Positive small-radius and transition exponent.
        b: Positive scale radius in micrometres. Mutually exclusive with
            ``median_radius``.
        gamma: Tail parameter greater than one.
        surface_area_density: Spherical surface-area density, as an alternative
            to ``n``. The second moment must exist.
        volume_density: Spherical volume density, as an alternative to ``n``.
            The third moment must exist.
        median_radius: Number-median radius in micrometres, as an alternative
            to ``b``.

    Attributes:
        median_radius: True number median.
        mode_radius: Interior mode for ``alpha > 1``; otherwise zero.
        effective_radius: ``M3 / M2`` when both moments exist, else ``np.inf``.
        effective_variance: ``M2 * M4 / M3**2 - 1`` when all required moments
            exist, else ``np.inf``.

    Notes:
        This is a scaled Burr type-XII distribution with SciPy parameters
        ``c=alpha`` and ``d=gamma-1``.

    Raises:
        ValueError: If parameters are outside their domains, both/neither scale
            representations are supplied, or a density input requires a
            divergent moment.
    """

    def __init__(
        self,
        n=None,
        alpha=None,
        b=None,
        gamma=None,
        surface_area_density=None,
        volume_density=None,
        *,
        median_radius=None,
    ):
        super().__init__("regularised_power_law")
        self.alpha = _positive_parameter("alpha", alpha)
        self.gamma = _positive_parameter("gamma", gamma)
        if self.gamma <= 1.0:
            raise ValueError("gamma must be greater than 1 for normalisation.")
        if (b is None) == (median_radius is None):
            raise ValueError("Provide exactly one of b or median_radius.")

        log_median_factor = (
            np.log(np.expm1(np.log(2.0) / (self.gamma - 1.0))) / self.alpha
        )
        if median_radius is not None:
            self.median_radius = _positive_parameter("median_radius", median_radius)
            self._log_b = np.log(self.median_radius) - log_median_factor
            self.b = _exp(self._log_b)
        else:
            self.b = _positive_parameter("b", b)
            self._log_b = np.log(self.b)
            self.median_radius = _exp(self._log_b + log_median_factor)
        if not np.isfinite(self.b) or self.b <= 0.0:
            raise ValueError("The requested parameters produce a non-finite b.")

        self.r = self.median_radius
        self.s = None
        second_moment = self._moment_per_particle(2)
        third_moment = self._moment_per_particle(3)
        _set_concentration(
            self,
            n,
            surface_area_density,
            volume_density,
            second_moment,
            third_moment,
        )
        self.mode_radius = (
            _exp(
                self._log_b
                + (
                    np.log(self.alpha - 1.0)
                    - np.log(1.0 + self.alpha * (self.gamma - 1.0))
                )
                / self.alpha
            )
            if self.alpha > 1.0
            else 0.0
        )
        self.effective_radius = (
            third_moment / second_moment if np.isfinite(third_moment) else np.inf
        )
        fourth_moment = self._moment_per_particle(4)
        self.effective_variance = (
            second_moment * fourth_moment / third_moment**2 - 1.0
            if np.isfinite(fourth_moment)
            else np.inf
        )
        self._log_prefactor = (
            np.log(self.n)
            + np.log(self.alpha)
            + np.log(self.gamma - 1.0)
            - self.alpha * self._log_b
        )

    def _moment_per_particle(self, order):
        """Return a normalized raw moment or ``np.inf`` when divergent.

        Args:
            order: Non-negative raw-moment order.

        Returns:
            ``b**order * Gamma(1 + order/alpha)`` multiplied by
            ``Gamma(gamma - 1 - order/alpha) / Gamma(gamma - 1)`` when the
            tail condition is satisfied; otherwise ``np.inf``.
        """
        qi = self.gamma - 1.0 - order / self.alpha
        if qi <= 0.0:
            return np.inf
        pi = 1.0 + order / self.alpha
        return _exp(
            order * self._log_b + gammaln(pi) + gammaln(qi) - gammaln(self.gamma - 1.0)
        )

    def moment(self, order):
        """Return a raw moment, or ``np.inf`` when it diverges.

        Args:
            order: Non-negative integer moment order.

        Returns:
            Raw moment including total number concentration.

        Raises:
            ValueError: If ``order`` is not a non-negative integer.
        """
        if not isinstance(order, (int, np.integer)) or order < 0:
            raise ValueError("order must be a non-negative integer.")
        return self.n * self._moment_per_particle(order)

    def mean(self):
        """Return the arithmetic mean, or ``np.inf`` when it diverges."""
        return self._moment_per_particle(1)

    def value(self, radii):
        """Evaluate differential number density at positive radii.

        Args:
            radii: Scalar or array-like particle radii in micrometres.

        Returns:
            Number density per micrometre. Scalar input returns ``float`` and
            array-like input returns a matching NumPy array. Values at
            non-positive radii are zero.
        """
        radius = np.asarray(radii, dtype=float)
        values = np.zeros_like(radius)
        positive = radius > 0.0
        log_radius = np.log(radius[positive])
        transition = self.alpha * (log_radius - self._log_b)
        log_value = (
            self._log_prefactor
            + (self.alpha - 1.0) * log_radius
            - self.gamma * np.logaddexp(0.0, transition)
        )
        with np.errstate(over="ignore", under="ignore"):
            values[positive] = np.exp(log_value)
        return _return_scalar_when_scalar(radii, values)


# Optional aliases for common alternative spellings/capitalization.
RegularizedPowerLawDistribution = RegularisedPowerLawDistribution
MultiModeLogNormalDistribution = MultimodeLogNormalDistribution


#  Selector
def create_distribution(dist_type, **kwargs):
    """Construct a supported particle-size distribution by name.

    Args:
        dist_type: One of ``"gaussian"``, ``"log_normal"``,
            ``"multimode_log_normal"``, ``"gamma"``, ``"modified_gamma"``,
            ``"inverse_modified_gamma"``, or ``"regularised_power_law"``.
            ``"regularized_power_law"`` is accepted as a spelling alias.
        **kwargs: Constructor arguments forwarded unchanged to the selected
            distribution class.

    Returns:
        Concrete :class:`SizeDistribution` instance selected by ``dist_type``.

    Raises:
        ValueError: If ``dist_type`` is unknown or constructor values are invalid.
        TypeError: If constructor arguments are not accepted by the selected class.
    """
    distribution_types = {
        "gaussian": GaussianDistribution,
        "log_normal": LogNormalDistribution,
        "multimode_log_normal": MultimodeLogNormalDistribution,
        "gamma": GammaDistribution,
        "modified_gamma": ModifiedGammaDistribution,
        "inverse_modified_gamma": InverseModifiedGammaDistribution,
        "regularised_power_law": RegularisedPowerLawDistribution,
        "regularized_power_law": RegularisedPowerLawDistribution,
    }
    try:
        distribution_class = distribution_types[dist_type]
    except KeyError:
        raise ValueError(f"Unknown distribution type: {dist_type}") from None
    return distribution_class(**kwargs)
