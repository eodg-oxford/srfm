Particle size distributions
===========================

The :mod:`srfm.size_distribution` module represents differential particle
number concentration per unit radius. If ``distribution.value(r)`` is
:math:`n(r)`, then

.. math::

   N_0 = \int n(r)\,dr,
   \qquad
   M_i = \int r^i n(r)\,dr.

For spherical particles, SRFM reports total surface-area and volume densities
as

.. math::

   A_0 = 4\pi M_2,
   \qquad
   V_0 = \frac{4\pi}{3}M_3.

Radii are expressed in micrometres, ``n`` in cm\ :sup:`-3`, differential
number density in cm\ :sup:`-3` :math:`\mu`m\ :sup:`-1`, surface-area density
in :math:`\mu`m\ :sup:`2` cm\ :sup:`-3`, and volume density in
:math:`\mu`m\ :sup:`3` cm\ :sup:`-3`.

Common interface
----------------

All concrete classes expose:

* ``value(radii)`` for differential number density;
* ``mean()`` for arithmetic mean radius;
* ``n``, ``surface_area_density`` and ``volume_density``;
* ``median_radius``, ``effective_radius`` and ``effective_variance``;
* ``moment(order)`` for raw moments supported by the distribution.

Except for multimode log-normal distributions, one of ``n``,
``surface_area_density`` or ``volume_density`` may set the concentration. If
more than one is supplied, precedence is number, surface area, then volume.
The multimode class accepts the corresponding quantity separately for every
mode because a mixture total does not determine the component proportions.

Gaussian distribution
---------------------

The untruncated density is

.. math::

   n(r) = \frac{N_0}{\sqrt{2\pi}\sigma_0}
   \exp\left[-\frac{(r-\mu_0)^2}{2\sigma_0^2}\right].

:class:`~srfm.size_distribution.GaussianDistribution` takes ``r`` as the mean
and ``s`` as the standard deviation, both in micrometres::

   from srfm.size_distribution import GaussianDistribution

   gaussian = GaussianDistribution(n=100.0, r=5.0, s=1.5)

The default is the mathematical Gaussian on the complete real line. It assigns
some concentration to negative, unphysical radii. Set ``truncate=True`` to
truncate at zero and renormalize::

   positive_gaussian = GaussianDistribution(
       n=100.0, r=0.5, s=1.0, truncate=True
   )

For the untruncated form, number and moments refer to the complete real line;
in particular, its volume density includes the signed third moment. Use the
truncated form when ``n`` must refer only to physical radii.

Log-normal distributions
------------------------

For one mode,

.. math::

   n(r) = \frac{N_0}{\sqrt{2\pi}\ln(S)r}
   \exp\left[-\frac{(\ln r-\ln r_m)^2}{2\ln^2(S)}\right].

For :class:`~srfm.size_distribution.LogNormalDistribution`, ``r`` is the
number-median radius and ``s`` is the geometric standard deviation::

   from srfm.size_distribution import LogNormalDistribution

   mode = LogNormalDistribution(n=100.0, r=0.3, s=1.7)

``s`` must exceed one. The limiting case ``s == 1`` is a delta distribution and
does not have a finite differential log-normal density.

:class:`~srfm.size_distribution.MultimodeLogNormalDistribution` accepts one
number concentration, median radius and geometric standard deviation for every
mode::

   from srfm.size_distribution import MultimodeLogNormalDistribution

   mixture = MultimodeLogNormalDistribution(
       n=[100.0, 20.0],
       r=[0.1, 1.0],
       s=[1.5, 2.0],
   )

The total density and every raw moment are sums of their component values. The
mixture median is solved numerically, and its effective radius is ``M3 / M2``;
it is not generally a number-weighted average of the component effective
radii.

Gamma distribution
------------------

:class:`~srfm.size_distribution.GammaDistribution` is parameterized by
effective variance ``s`` and either the number median ``r`` or
``effective_radius``::

   from srfm.size_distribution import GammaDistribution

   gamma_distribution = GammaDistribution(
       n=100.0,
       effective_radius=5.0,
       s=0.1,
   )

Its conventional shape and rate are

.. math::

   \alpha = \frac{1}{v_e}-3,
   \qquad
   b = \frac{1}{r_e v_e}.

Normalization requires :math:`0 < v_e < 1/2`. An interior mode requires the
stronger condition :math:`v_e < 1/3`.

Modified gamma distribution
---------------------------

The normalized density is

.. math::

   n(r) = N_0\frac{\gamma b^{(\alpha+1)/\gamma}}
   {\Gamma((\alpha+1)/\gamma)}r^\alpha\exp(-br^\gamma).

:class:`~srfm.size_distribution.ModifiedGammaDistribution` adds the positive
cutoff exponent ``gamma``::

   from srfm.size_distribution import ModifiedGammaDistribution

   modified = ModifiedGammaDistribution(
       n=100.0,
       effective_radius=5.0,
       s=0.05,
       gamma=2.0,
   )

At fixed ``gamma``, SRFM solves the effective-variance equation numerically for
``alpha`` and derives ``b`` from the supplied radius scale. The ``r`` argument
is a true number median, consistent with the log-normal API. In the source
reference, :math:`r_m` in the modified-gamma mode-radius formula denotes a
*mode*, not a median; SRFM exposes that separate value as ``mode_radius``.

Inverse modified gamma distribution
-----------------------------------

The normalized density is

.. math::

   n(r) = N_0\frac{\gamma b^{(\alpha-1)/\gamma}}
   {\Gamma((\alpha-1)/\gamma)}r^{-\alpha}\exp(-br^{-\gamma}).

:class:`~srfm.size_distribution.InverseModifiedGammaDistribution` uses its
natural shape parameters and either ``b`` or ``median_radius``::

   from srfm.size_distribution import InverseModifiedGammaDistribution

   inverse = InverseModifiedGammaDistribution(
       n=100.0,
       alpha=8.0,
       b=0.01,
       gamma=2.0,
   )

Normalization requires ``alpha > 1``. Raw moment ``i`` exists only when
``i < alpha - 1``. Thus finite surface area requires ``alpha > 3``, finite
volume and effective radius require ``alpha > 4``, and finite effective
variance requires ``alpha > 5``. SRFM reports divergent moments and derived
metrics as ``np.inf``.

Regularised power law
---------------------

The normalized density is

.. math::

   n(r) = N_0\alpha(\gamma-1)b^{-\alpha}r^{\alpha-1}
   [1+(r/b)^\alpha]^{-\gamma}.

:class:`~srfm.size_distribution.RegularisedPowerLawDistribution` takes
positive ``alpha``, scale radius ``b`` and ``gamma > 1``::

   from srfm.size_distribution import RegularisedPowerLawDistribution

   power_law = RegularisedPowerLawDistribution(
       n=100.0,
       alpha=3.0,
       b=0.8,
       gamma=4.0,
   )

``median_radius`` may replace ``b``. Moment ``i`` exists only when
``gamma > 1 + i / alpha``; divergent moments are reported as ``np.inf``. The
class is mathematically equivalent to a scaled Burr type-XII distribution with
``c=alpha`` and ``d=gamma-1``.

Factory construction
--------------------

:func:`~srfm.size_distribution.create_distribution` recognizes:

* ``"gaussian"``;
* ``"log_normal"``;
* ``"multimode_log_normal"``;
* ``"gamma"``;
* ``"modified_gamma"``;
* ``"inverse_modified_gamma"``;
* ``"regularised_power_law"``;
* ``"regularized_power_law"`` as a spelling alias.

Keyword arguments are forwarded to the selected class::

   from srfm.size_distribution import create_distribution

   distribution = create_distribution(
       "gamma", n=100.0, effective_radius=5.0, s=0.1
   )

Numerical validation
--------------------

The unit tests compare densities, raw moments zero through four, means,
medians, surface-area densities and volume densities with independent SciPy
distributions. The correspondences are ``norm`` and ``truncnorm`` for the two
Gaussian forms, ``lognorm`` for each log-normal mode, ``gamma`` for the
ordinary gamma distribution, ``gengamma`` for the modified gamma distribution,
``invgamma`` for the ``gamma=1`` inverse modified-gamma case, and ``burr12``
for the regularised power law. Additional quadrature tests check that every
positive-radius density integrates to its requested number concentration.

Mie-layer limitation
--------------------

These classes and the factory are available through the direct Python API. The
current :class:`~srfm.layer.MieLayer` radius-grid code derives integration
limits specifically from the single-mode log-normal ``r`` and ``s``
parameters. Consequently, the driver-table scattering-layer pathway currently
supports only ``"log_normal"`` for complete optical calculations. Supporting
the other distributions there requires distribution-specific or generic
radius-bound selection and additional schema fields for ``alpha``, ``b`` and
``gamma``.

Reference
---------

The formulae follow R. G. Grainger, *Some Useful Formulae for Particle Size
Distributions and their Optical Properties*, revised 27 August 2026,
`doi:10.5281/zenodo.22101058 <https://doi.org/10.5281/zenodo.22101058>`_.
