Memory-efficient calculations
=============================

The top-level SRFM runners keep particle optical properties on their coarse Mie
grids and interpolate them in spectral blocks. They do not interpolate phase
functions, retain per-wavenumber DISORT dictionaries, or construct the legacy
wide RFM optical-depth DataFrame. These choices make the main temporary storage
scale with ``scattering_block_size`` instead of the complete spectral grid.

Output retention
----------------

Set ``retain_outputs`` to the values needed after the run. For example, a
brightness-temperature-only calculation uses::

   "retain_outputs": ("bbt",),

Radiance is retained temporarily because brightness temperature and IASI
convolution require it, then released when it was not requested.

With NetCDF output, explicitly retained values are saved together in one file,
named ``srfm.nc`` by default. For example,
``"retain_outputs": ("bbt", "flup")`` writes
``flup(wavenumber, output_level)`` alongside BBT. A raw result such as
``"retain_outputs": ("rfldn",)`` is also valid without BBT or radiance.
User-angle radiance uses all four spectral geometry dimensions, while
``albmed`` and ``trnmed`` use ``(wavenumber, output_polar_angle)``.

Compatibility interfaces
------------------------

``MieLayer.regrid`` retains its historical behavior by default. Standalone
callers can pass ``regrid_phase_function=False`` and ``retain_original=False``
to avoid large arrays. ``OpticalDepthGrid.to_dataframe`` and
``get_captured_optical_depths`` provide explicit legacy DataFrame conversion.
Standalone ``DISORT`` objects retain history by default; top-level runners use
``retain_history=False`` and consume the returned ``DisortResult`` directly.
The bundled f2py interface requires the reusable ``pmom`` workspace in C order;
passing a Fortran-contiguous production matrix fails DISORT's value checks, so
the wrapper's input conversion is currently unavoidable.

Benchmarking
------------

The benchmark tool samples aggregate RSS for a command and all of its child
processes and writes machine-readable JSON::

   python tools/benchmark_memory.py --output memory.json -- \
       python examples/basic_example/run_srfm.py

Run the same command and inputs against two revisions, then compare
``peak_rss_bytes`` and ``elapsed_seconds``. A medium grid is appropriate for
routine checks; production-scale grids should be run on a machine with the
corresponding time and memory budget.
