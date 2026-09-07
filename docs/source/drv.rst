Driver table
============

As explained above, the SRFM is a package which can be executed as a standard 
radiative transfer model with a driver table.
This section describes all variables present in the driver table. 
The driver table is a Python file containing an ``inputs`` dictionary.  The
example at ``examples/basic_example/driver_table.py`` also defines named grid
constants for readability; the parameters below are the values passed to SRFM.
Paths may be strings or :class:`os.PathLike` objects unless noted otherwise.

General output
--------------

* ``base_plots`` (``bool``): Create one output spectrum plot for every polar
  angle, output level, and azimuthal angle. Permitted values: ``True`` or
  ``False``.
* ``out_mode`` (``str`` or ``None``): Select the saved spectrum format.
  Permitted values: ``"netcdf"``, ``"txt"``, or ``None`` for no spectrum file.
* ``show_plots`` (``bool``): Display generated plots interactively as well as
  saving them. Permitted values: ``True`` or ``False``.
* ``results_fldr`` (path-like): Directory in which SRFM writes RFM products,
  spectra, and plots.
* ``rad`` (``bool``): Save radiance output when file output is enabled.
  Permitted values: ``True`` or ``False``.
* ``bbt`` (``bool``): Save brightness-temperature output when file output is
  enabled. Permitted values: ``True`` or ``False``.
* ``rad_out_fname`` (``str`` or ``None``): Radiance filename stem. ``None``
  selects ``rad``.
* ``bbt_out_fname`` (``str`` or ``None``): Brightness-temperature filename
  stem. ``None`` selects ``bbt``.
* ``plot_type`` (``str``): Quantity shown in base plots. Permitted values:
  ``"rad"`` or ``"bbt"``.
* ``convolve_iasi`` (``bool``): Convolve every output radiance spectrum with
  the IASI instrument line shape. Permitted values: ``True`` or ``False``.
* ``iasi_ils`` (path-like or ``None``): RFM-format IASI instrument-line-shape
  file; required when ``convolve_iasi`` is ``True``.
* ``retain_outputs`` (sequence of ``str`` or ``None``): Full-grid values kept
  on the returned model. Supported names are ``"radiance"``/``"uu"``,
  ``"bbt"``, ``"rfldir"``, ``"rfldn"``, ``"flup"``, ``"dfdt"``,
  ``"uavg"``, ``"albmed"``, and ``"trnmed"``. Values needed for a requested
  file or plot are retained automatically. Omitting the parameter preserves
  the historical all-output behavior.
* ``scattering_block_size`` (``int``): Positive number of spectral points for
  which particle properties are interpolated at once. The default is 10,000.
* ``retain_phase_functions`` (``bool``): Keep coarse particle phase functions
  after Legendre expansion for inspection. The default is ``False`` in the
  runners.

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

* ``fisot`` (``int`` or ``float``): Isotropic illumination incident at the top
  of the atmosphere; must be non-negative.
* ``albedo`` (``int`` or ``float``): Lambertian bottom-boundary albedo.
  Permitted range: 0 to 1.
* ``temis`` (``int`` or ``float``): Boundary emissivity used by DISORT.
  Permitted range: 0 to 1.
* ``earth_radius`` (``int`` or ``float``): Positive Earth radius in kilometres,
  used by the pseudo-spherical correction.
* ``nmom`` (``int``): Requested number of phase-function moments. SRFM raises
  it when a scattering layer or stream count requires more moments.
* ``maxcmu`` (``int``): Positive, even number of DISORT computational streams;
  must be at least 2.
* ``maxumu`` (``int``): Number of output polar-angle cosines allocated by
  DISORT. It must be positive unless ``onlyfl`` is ``True``.
* ``maxphi`` (``int``): Number of output azimuthal angles allocated by DISORT.
  It must be positive unless ``onlyfl`` is ``True``.
* ``maxulv`` (``int``): Compatibility input for DISORT output-level allocation.
  Every top-level runner derives the effective value from ``out`` and
  ``out_toa``.
* ``usrang`` (``bool``): Return output at user-defined angles. Permitted values:
  ``True`` or ``False``.
* ``usrtau`` (``bool``): Return output at user-defined optical depths. Every
  top-level runner requires ``True`` when using ``out``.
* ``ibcnd`` (``int``): DISORT boundary-condition mode. Permitted values: ``0``
  for a general atmosphere or ``1`` for the special flux/albedo problem.
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
* ``header`` (``str``): Header sent to DISORT terminal output. ``"NO HEADER"``
  suppresses the header; fewer than 127 characters are permitted.
* ``adjust_maxcmu`` (``bool``): Retry with more streams when Delta-M+ produces
  a small negative intensity. Permitted values: ``True`` or ``False``.

Scattering-layer configuration
------------------------------

``scat_lyrs_inputs`` is a dictionary keyed by unique layer name. Each value is
a dictionary with the parameters below. The same fields apply to the example's
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
  not implemented by the full Mie-layer pathway.
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
* ``r`` (``int`` or ``float``): Positive mean particle radius in micrometres.

Solar and viewing geometry
--------------------------

* ``sun`` (``bool``): Include direct solar illumination and reflection.
  Permitted values: ``True`` or ``False``.
* ``sza`` (``int`` or ``float``): Solar zenith angle in degrees. Permitted
  range: 0 to 180.
* ``saa`` (``int`` or ``float``): Solar azimuth angle in degrees. Permitted
  range: 0 to 360.
* ``zen`` (``int`` or ``float``): Viewing zenith angle in degrees. Permitted
  range: 0 to 180.
* ``azi`` (``int`` or ``float``): Viewing azimuth angle in degrees. Permitted
  range: 0 to 360.
