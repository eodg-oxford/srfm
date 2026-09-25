Driver table
============

SRFM can be executed as a standard radiative transfer model with a driver
table. This page is the complete parameter reference; the
:doc:`worked example <auto_examples/basic_example/run_srfm>` shows these
parameters together in an executable calculation.

The driver table is a Python file containing an ``inputs`` dictionary. Fields
marked **optional** may be omitted, in which case SRFM uses the stated default
or derives the value from the atmosphere. Other fields are required by
``main.run_srfm``. The example driver also defines named grid constants for
readability; the parameters below are the values passed to SRFM. Paths may be
strings or :class:`os.PathLike` objects unless noted otherwise.

General run and output
----------------------

* ``date`` (**optional** :class:`datetime.datetime` or three-item ``tuple``):
  UTC calculation date, or ``(year, month, day)``. It determines the
  Sun--Earth distance correction when solar illumination is enabled. If
  omitted, SRFM uses 23 March 2025, representing approximately the mean
  Sun--Earth distance.

* ``base_plots`` (``bool``): Create one output spectrum plot for every polar
  angle, output level, and azimuthal angle. Permitted values: ``True`` or
  ``False``.
* ``out_mode`` (``str`` or ``None``): Select the saved spectrum format.
  Permitted values: ``"netcdf"``, ``"txt"``, or ``None`` for no spectrum file.
* ``show_plots`` (``bool``): Display generated plots interactively as well as
  saving them. Permitted values: ``True`` or ``False``.
* ``results_fldr`` (path-like): Directory in which SRFM writes RFM products,
  spectra, and plots.
* ``out_fname`` (**optional** ``str`` or ``None``): Filename for the single
  NetCDF output file. Omitting it or setting it to ``None`` selects
  ``srfm.nc``. The value is used as supplied, including its extension, and
  applies only when ``out_mode`` is ``"netcdf"``.
* ``convolve_iasi`` (``bool``): Convolve every output radiance spectrum with
  the IASI instrument line shape. Permitted values: ``True`` or ``False``.
* ``iasi_ils`` (path-like or ``None``): RFM-format IASI instrument-line-shape
  file; required when ``convolve_iasi`` is ``True``.
* ``retain_outputs`` (required sequence of ``str``): Full-grid values kept on
  the returned model. Supported names are ``"rad"``/``"radiance"``/``"uu"``,
  ``"bbt"``, ``"rfldir"``, ``"rfldn"``, ``"flup"``, ``"dfdt"``,
  ``"uavg"``, ``"albmed"``, and ``"trnmed"``. In NetCDF mode, the single
  output file contains exactly these requested SRFM result variables together
  with their coordinates and optical-layer metadata. Layer type, computational
  optical depth, particle SSA and normalized moments, surface albedo, and the
  FBEAM spectrum are stored as NetCDF variables where applicable. Large source
  grids and values are summarized rather than embedded in the JSON global
  attribute. Raw
  DISORT results such as ``"rfldn"`` can be written without retaining BBT or
  radiance. Text mode requires retained BBT or radiance and writes them to
  ``bbt.txt`` and ``rad.txt`` respectively. When ``base_plots`` is true, all
  retained BBT and radiance outputs are plotted; if neither is retained, SRFM
  warns and creates no base plots.
* ``scattering_block_size`` (**optional** ``int``): Positive number of spectral
  points for which particle properties are interpolated at once. The default
  is 10,000.
* ``retain_phase_functions`` (**optional** ``bool``): Keep coarse particle
  phase functions after Legendre expansion for inspection. The default is
  ``False``.

Spectral grids
--------------

* ``fin_wvnmlo`` (``int`` or ``float``): Lower wavenumber of the final output
  grid in cm\ :sup:`-1`.
* ``fin_wvnmhi`` (``int`` or ``float``): Upper wavenumber of the final output
  grid in cm\ :sup:`-1`; must be greater than ``fin_wvnmlo``.
* ``fin_res`` (``int`` or ``float``): Positive final-grid spacing in
  cm\ :sup:`-1`.
* ``spc_wvnmlo`` (``int`` or ``float``): Lower bound of the computational
  spectral grid, in ``spc_units``.
* ``spc_wvnmhi`` (``int`` or ``float``): Upper bound of the computational
  spectral grid, in ``spc_units``; must be greater than ``spc_wvnmlo``.
* ``spc_res`` (``int`` or ``float``): Positive computational-grid spacing, in
  ``spc_units``.
* ``spc_units`` (``str``): Units of the computational grid. Permitted values:
  ``"cm-1"``, ``"um"``, or ``"nm"``.

Output geometry
---------------

* ``out_fmt`` (``str``): Interpretation of ``out``. Permitted values:
  ``"altitude"`` for altitude in kilometres or ``"tau"`` for optical depth.
* ``out`` (``int``, ``float``, or one-dimensional array-like): Non-negative
  output altitudes or optical depths. Altitudes are matched to the nearest
  atmospheric level boundary. Values retain their supplied order.
* ``out_toa`` (``bool``): Add top-of-atmosphere output unless it is already in
  ``out``. Permitted values: ``True`` or ``False``.

Atmospheric calculation levels
------------------------------

* ``levels`` (**optional** one-dimensional sequence): Strictly increasing
  atmospheric calculation levels in kilometres. Supply at least two finite
  numeric values. If omitted, SRFM uses its built-in standard level grid.
  Configured optical-layer boundaries are inserted into the effective grid
  automatically.

RFM execution configuration
---------------------------

The following keys belong to the nested ``rfm_config`` dictionary.

* ``output_mode`` (``str``): How RFM output is obtained. The generic SRFM
  runner requires ``"capture"``; ``"files"`` is supported by the lower-level
  RFM helper only.
* ``driver_path`` (path-like or ``None``): Existing RFM driver-file path, or
  ``None`` when the driver is assembled from ``driver_inputs``.
* ``generate_driver`` (``bool``): Permit generation or replacement of
  ``driver_path`` by the RFM helper. Permitted values: ``True`` or ``False``.
* ``verbose`` (``bool``): Print detailed RFM execution information. Permitted
  values: ``True`` or ``False``.
* ``capture_files_content`` (``bool``): Retain raw generated-file contents in
  the RFM result in addition to processed captured output. Permitted values:
  ``True`` or ``False``.
* ``optical_depth_format`` (``str``): Internal captured representation,
  ``"compact"`` or ``"dataframe"``. The top-level runners default to the
  compact spectral-first array. The lower-level compatibility API continues
  to support the DataFrame representation.

RFM driver sections
-------------------

The following keys belong to ``driver_inputs`` and map to sections of an RFM
driver file.

* ``header`` (``str`` or sequence of ``str``): Descriptive RFM driver header.
* ``flags`` (sequence of ``str``): Three-character RFM feature flags. The SRFM
  optical-depth workflow requires ``"OPT"`` and ``"LEV"``. The example also
  uses ``"NAD"``, ``"SFC"``, ``"PRF"``, ``"DBL"``, ``"CHI"``, and ``"MIX"``.
* ``spectral`` (sequence of ``rfm_helper.SpectralRange`` or
  ``rfm_helper.SpectralFile``): Spectral ranges or grid files supplied to RFM.
  ``run_srfm`` replaces the example value with its generated computational
  grid file.
* ``gases`` (sequence of ``str``): RFM gas identifiers included in gas
  absorption calculations.
* ``atmosphere`` (sequence of path-like): RFM atmosphere profile files. In the
  example, the altitude-grid profile is followed by the physical atmosphere.
* ``hit`` (sequence of path-like): HITRAN line-data files used by RFM.
* ``xsc`` (sequence of path-like): Directories or files containing molecular
  cross-section data used by RFM.

DISORT configuration
--------------------

The basic example exposes the following inputs that are passed directly to
DISORT: ``fisot``, ``albedo``, ``temis``, ``earth_radius``, ``nmom``,
``maxcmu``, ``maxumu``, ``maxphi``, ``usrang``, ``usrtau``,
``ibcnd``, ``onlyfl``, ``prnt``, ``planck``, ``lamber``, ``deltamplus``,
``do_pseudo_sphere``, ``header``, ``btemp``, and ``ttemp``.
``disort_precision`` selects the compiled wrapper and ``adjust_maxcmu`` controls
SRFM's retry behaviour. The geometry fields ``sun``, ``sza``, ``saa``, ``zen``,
and ``azi`` determine DISORT's beam strength and angular inputs. Finally,
``out_fmt``, ``out``, and ``out_toa`` determine its output optical depths, while
the atmosphere and configured layers supply optical depth, single-scatter
albedo, temperature, and phase moments.

SRFM validates the user-controlled values before starting RFM or DISORT. A
failure identifies the field, supplied value, and permitted values, for example
``maxcmu: current value 5; permitted values: an even integer greater than or
equal to 4``. This prevents DISORT's fatal ``VAR in error`` path for malformed
driver values.

* ``fisot`` (``int`` or ``float``): Isotropic illumination incident at the top
  of the atmosphere; must be non-negative.
* ``albedo`` (``int``, ``float``, or spectral-field mapping): Lambertian
  bottom-boundary albedo. A scalar retains the historical behavior. A mapping
  contains either ``grid``/``values`` arrays or the standardized text-file
  fields documented in :doc:`input_files`. Values must be finite and in [0, 1]
  and the source must cover the computational grid.
* ``temis`` (``int`` or ``float``): Boundary emissivity used by DISORT.
  Permitted range: 0 to 1.
* ``earth_radius`` (**optional** ``int`` or ``float``): Positive Earth radius
  in kilometres, used by the pseudo-spherical correction. The default is
  6371 km.
* ``nmom`` (``int``): Non-negative requested number of phase-function moments.
  SRFM raises it when a scattering layer or stream count requires more moments.
* ``maxcmu`` (``int``): Even number of DISORT computational streams, at least
  4. Although the native bounds check admits 2 streams, this bundled DISORT then
  raises a fatal ``2 streams not recommended`` error, so SRFM rejects 2 early.
* ``maxumu`` (``int``): Number of output polar-angle cosines allocated by
  DISORT. It must be positive unless ``onlyfl`` is ``True``. When ``usrang`` is
  ``False`` and angular intensities are requested, it must be at least
  ``maxcmu`` to hold the computational angles.
* ``maxphi`` (``int``): Number of output azimuthal angles allocated by DISORT.
  It must be positive unless ``onlyfl`` is ``True``.
* ``usrang`` (``bool``): Return output at user-defined angles. Permitted values:
  ``True`` or ``False``.
* ``usrtau`` (``bool``): Return output at user-defined optical depths. Every
  top-level runner requires ``True`` when using ``out``.
* ``ibcnd`` (``int``): DISORT boundary-condition mode. Permitted values: ``0``
  for a general atmosphere or ``1`` for the special albedo/transmissivity
  problem. Native DISORT requires ``onlyfl=False`` in mode 1.
* ``onlyfl`` (``bool``): Return fluxes only instead of angular intensities.
  Permitted values: ``True`` or ``False``.
* ``prnt`` (list of five ``bool``): DISORT output-section print switches.
* ``planck`` (``bool``): Include internal Planck thermal emission. Permitted
  values: ``True`` or ``False``.
* ``lamber`` (``bool``): Treat the bottom reflector as Lambertian. Permitted
  values: ``True`` or ``False``.
* ``deltamplus`` (``bool``): Use the Delta-M+ phase-function approximation when
  ``True``, or standard Delta-M when ``False``.
* ``do_pseudo_sphere`` (``bool``): Apply the pseudo-spherical solar-beam
  correction. Permitted values: ``True`` or ``False``.
* ``disort_precision`` (``str``): Compiled DISORT implementation to call.
  Permitted values: ``"single"`` or ``"double"``.
* ``header`` (**optional** ``str``): Header sent to DISORT terminal output.
  ``"NO HEADER"`` suppresses the header and is the default; at most 127
  characters are permitted, including an empty string.
* ``adjust_maxcmu`` (``bool``): Retry with more streams when Delta-M+ produces
  a small negative intensity. Permitted values: ``True`` or ``False``.
* ``btemp`` (**optional** ``int`` or ``float``): Bottom-boundary temperature in
  kelvin; it must be non-negative. If omitted, SRFM derives it from the lowest
  atmospheric temperature.
* ``ttemp`` (**optional** ``int`` or ``float``): Top-boundary temperature in
  kelvin; it must be non-negative. If omitted, SRFM derives it from the highest
  atmospheric temperature.

Scattering-layer configuration
------------------------------

``scat_lyrs_inputs`` is an **optional** dictionary keyed by unique layer name.
Omit it, set it to ``None``, or use an empty dictionary for a clear-sky
calculation. Each value is a dictionary with the
parameters below. The same fields apply to the example's
``Sulphuric_acid_1``, ``Ash_1``, and ``Water_cloud_1`` layers.

* ``name`` (``str``): Layer name, normally identical to its dictionary key.
* ``low_spc`` (``int`` or ``float``): Lower spectral bound for the layer's Mie
  calculation, in ``spec_units``.
* ``upp_spc`` (``int`` or ``float``): Upper spectral bound; must exceed
  ``low_spc``.
* ``res`` (``int`` or ``float``): Positive spacing of the layer optical-property
  grid, in ``spec_units``.
* ``spec_units`` (``str``): Units of the layer optical-property grid. Permitted
  values: ``"cm-1"``, ``"um"``, or ``"nm"``.
* ``rho`` (``int``, ``float``, or ``str``): Particle density in kg m\ :sup:`-3`,
  or a supported named ash density. Named values: ``"pumice"``, ``"glass"``,
  ``"mineral"``, or ``"rock"``.
* ``n`` (``int``, ``float``, or ``None``): Particle number concentration. At
  least one of ``mass_loading``, ``n``, ``s_a_den``, or ``v_den`` must be set.
* ``s`` (``int`` or ``float``): Particle-size distribution spread. A
  log-normal distribution requires a value greater than 1.
* ``s_a_den`` (``int``, ``float``, or ``None``): Particle surface-area density.
* ``v_den`` (``int``, ``float``, or ``None``): Particle volume density.
* ``dist_type`` (``str``): Particle-size distribution. Permitted value for a
  complete calculation: ``"log_normal"``. ``"gaussian"`` is recognized but
  not supported by the full Mie-layer pathway. Additional analytic
  distributions are available through the direct Python API; see
  :doc:`size_distributions`.
* ``comp`` (``str``): Refractive-index composition identifier or refractive-
  index filename, such as ``"sulphuric acid"``, ``"ash"``, or an ``.ri`` file.
* ``center_alt`` (``int``, ``float``, or ``None``): Layer centre altitude in km.
  Supply it with ``thick``, or instead supply both explicit altitude bounds.
* ``thick`` (``int``, ``float``, or ``None``): Positive layer thickness in km.
* ``alt_upp`` (``int``, ``float``, or ``None``): Explicit upper layer boundary
  in km; must exceed ``alt_low``.
* ``alt_low`` (``int``, ``float``, or ``None``): Explicit lower layer boundary
  in km.

* ``radii`` (``int``): Positive number of particle radii used to integrate the
  size distribution.
* ``eta`` (``float``): Relative size-distribution tail cutoff. Permitted range:
  greater than 0 and less than 1.
* ``phase_quad_N`` (``int``): Positive number of phase-function quadrature
  points.
* ``phase_quad_type`` (``str``): Phase-function quadrature rule. Permitted
  values: ``"G"``, ``"R"``, ``"L"``, or ``"T"``.
* ``radii_quad_type`` (``str``): Particle-radius quadrature rule. Permitted
  values: ``"G"``, ``"R"``, ``"L"``, or ``"T"``.
* ``leg_coeffs`` (``bool``): Calculate Legendre expansion coefficients.
  Permitted values: ``True`` or ``False``.
* ``leg_coeffs_type`` (``str``): Legendre coefficient convention. Permitted
  values: ``"regular"`` or ``"normalised"``.
* ``multiprocess`` (``bool``): Enable the experimental multiprocessing path for
  optical-property calculations. Permitted values: ``True`` or ``False``.
* ``mass_loading`` (``int``, ``float``, or ``None``): Column particle mass in
  g m\ :sup:`-2` used to derive layer optical depth.
* ``r`` (``int`` or ``float``): Positive number-median particle radius in
  micrometres for the supported ``"log_normal"`` distribution.

The unused altitude pair may be omitted entirely. If all four altitude values
are supplied, both pairs must describe the same boundaries after SRFM's
0.001-km layer-boundary rounding.

``gbc_lyrs_inputs`` is likewise optional and accepts omission, ``None``, or an
empty dictionary. Its grey-body layers use the same alternative altitude pairs:
``center_alt`` with ``thick``, or ``alt_low`` with ``alt_upp``. When all four
are supplied they must describe the same layer.

Prescribed optical layers
-------------------------

``prescribed_lyrs_inputs`` is a separate optional mapping. It is not part of
``scat_lyrs_inputs`` because its values are already column optical properties,
not Mie microphysics. Each layer requires ``name``, finite ``alt_low`` and
``alt_upp`` in kilometres with ``alt_low < alt_upp``, ``optical_depth``, and
``ssalb``.

An Angstrom optical depth is configured as::

   "optical_depth": {
       "type": "angstrom",
       "reference_value": 0.00058,
       "reference_wavelength_um": 1.0,
       "angstrom_exponent": 1.48,
   }

SRFM evaluates
:math:`\tau(\lambda)=\tau_{ref}(\lambda/\lambda_{ref})^{-\alpha}` with
wavelength in micrometres. The reference value must be non-negative, reference
wavelength finite and positive, and exponent finite. Alternatively,
``type="tabulated"`` accepts a common in-memory or file-backed spectral field.

``ssalb`` accepts a finite scalar in [0, 1], broadcast across the calculation,
or a common spectral field. A bare one-dimensional array is deliberately not
accepted because SRFM cannot infer its coordinate.

Henyey--Greenstein phase data accept scalar or spectral asymmetry::

   "phase_function": {
       "type": "henyey_greenstein",
       "asymmetry": 0.490,
   }

The asymmetry factor lies in [-1, 1]. Normalized moments are calculated as
:math:`\beta_l=g^l` from order zero through ``nmom``. Fully tabulated moments
use ``type="legendre_moments"``, ``convention="normalised"``, a spectral grid,
and values shaped ``(spectral_points, moments)``. File-backed moments use
``value_columns`` as described in :doc:`input_files`. Moment zero equals one
within ``1e-5`` wherever scattering optical depth ``tau * ssalb`` is positive.
An absorption-only prescribed layer may omit phase data.

Mie and prescribed particle objects expose the same blockwise column optical
depth, single-scattering albedo, and normalized-moment contract. Properties are
interpolated in bounded blocks controlled by ``scattering_block_size``.

Vertical optical-layer restrictions
-----------------------------------

Names are unique across ``scat_lyrs_inputs``, ``prescribed_lyrs_inputs``, and
``gbc_lyrs_inputs``. Their resolved altitude intervals must be separated: two
objects may neither overlap by a positive distance nor share an exact boundary.
Validation is independent of dictionary order and occurs before layer optical
calculations, RFM, or DISORT.

Each accepted object replaces every requested internal level inside its bounds
with its lower and upper boundaries. It therefore occupies exactly one model
cell and contributes its complete column optical depth to that cell. SRFM does
not fractionally split, combine, or mix configured optical objects.

Solar and viewing geometry
--------------------------

* ``sun`` (``bool``): Include direct solar illumination and reflection.
  Permitted values: ``True`` or ``False``.
* ``solar_spectrum`` (**optional** spectral-field mapping): Custom beam-normal
  DISORT ``FBEAM`` spectral irradiance. Its required units, conversions, file
  form, and no-scaling convention are documented in :doc:`input_files`.
  Omitting it selects the established Gueymard-2018 spectrum with date-based
  Earth--Sun-distance scaling. With ``sun=False``, FBEAM is zero.
* ``sza`` (``int`` or ``float``): Solar zenith angle in degrees. Permitted
  range: 0 to 180. When ``sun`` is ``True``, DISORT requires a positive incident
  beam cosine, so the permitted range is 0 inclusive to 90 exclusive. Set
  ``sun=False`` for night-side geometries.
* ``saa`` (``int`` or ``float``): Solar azimuth angle in degrees. Permitted
  range: 0 to 360.
* ``zen`` (``int`` or ``float``): Viewing zenith angle in degrees. Permitted
  range: 0 to 180. A value of exactly 90 is not permitted when user-angle
  intensities are requested because DISORT forbids a zero output-angle cosine.
  Mode ``ibcnd=1`` additionally requires a positive cosine, so user angles must
  be less than 90 degrees.
* ``azi`` (``int`` or ``float``): Viewing azimuth angle in degrees. Permitted
  range: 0 to 360.
