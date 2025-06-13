import numpy as np


class Read_Iasi_L1c:
    """Read IASI L1C nat file

    VERSION:
      - 15JAN24 AD: Fix some bugs associated with renaming.
      - 31OCT23 AD: 
        Renamed from Iasi_L1c to Read_Iasi_L1c.
        Prefix private method names with '_'.
      - 18DEC20 AD: Original version
      - 15APR25 AK: Updated documentation (formatting).
      - 01JUN25 AK: Added saa and azi as parameters to be extracted/read.
      - 11JUN25 AK: Changed self.spec dtype to np.ndarray (previously list).

    USAGE:
      By default this will read the entire IASI L1C file contents, c.90000 spectra.
      However various options are available for the user to subset the data.
      By default spectra are radiances in units of nW/(cm2.sr.cm-1) but calling
      with the option bright=True returns spectra as brightness temperatures.

      Typically call as:
      
        l1c = Read_Iasi_L1c(filename, [options] )
        # check for errors
        if l1c.fail:                      
          print(l1c.errmsg)
          stop

      An alternative method is just to extract the location data and use the
      get_spec method to extract individual spectra:
      
        l1c = Read_Iasi_L1c(filename, loc_only=True, [other options] )
        if l1c.fail ...etc
        for iloc in range(l1c.nloc)
          rad = l1c.get_spec(iloc,[options])
       (options for get_spec method are wnolim, bright)

    OPTIONS:
      avhrr (bool): True = include AVHRR cluster analysis data. Default is False.
      bright (bool): True for brightness temp spectra, False for radiance spectra.
        Default is False.
               
      azilim (tuple): (azimin, azimax), Default is None.
      chkqal (list of bool): Default is [True,True,True].
        If True, exclude data with bad quality flag for each band.
        (bands are 645-1210, 1210-2000, 2000-2760 cm-1).
      cldlim (tuple): (cldmin, cldmax), Default is None.
        min,max cloud cover percentages for inclusion, (0,0)=only cloud-free
      latlim (tuple): (latmin, latmax). Default is None.
        min,max latitudes (range -90:90) for inclusion.
      lndlim (tuple): (lndmin, lndmax), Default is None.
        min,max land (cf.ocean) percentages for inclusion, (0,0)=only ocean.
      loc_only (bool); Default is False.
        Set True if only location data (without spectra) are required
      lonlim (tuple): (lonmin, lonmax), Default is None.
        min,max longitudes (range -180:180) for inclusion.
        lonmin,lonmax are actually east,west boundaries so, for example,
        (-170,170) would include everything apart from the 20deg lon band
        around +/-180 deg, whereas (170,-170) would only include this
        20deg band.
        The same actually applies to all other limits, but is probably
        only useful for longitude.
      mph_only (bool): Default is False.
        Set True if just the Main Product Header data is required (as self.mph).
      saalim (tuple): (saamin, saamax), Default is None.
      szalim (tuple): (szamin,szamax), Default is None.
        min,max solar zenith angle (0,90=daytime, 90-180=nighttime).
      wnolim (tuple): (wnomin, wnomax), Default is (645,2760).
        min,max wavenumber [cm-1] limits for extracted spectra.
      zenlim (tuple): (zenmin,zenmax), Default is None.
        min,max satellite zenith angle (0=sat overhead).

    LIST OF METHODS:

      The following methods are 'private' and called on initialisation:
        - __init__     : Open file and read location data
        - _btspec       : Convert radiance to brightness temperature
        - _outside      : Determine if parameter outside selected limits
        - _read_ipr     : Read Internal Pointer Records
        - _read_loc     : Read pixel location data
        - _read_mdr     : Check Measurement Data Records and get offsets
        - _read_mph     : Read Main Product Header
        - _read_sca     : Read radiance spectra scale factors
        - _read_spec    : Read in spectra at selected locations
      
      The following methods are 'public' - intended for the user:
        - get_spec     : Extract individual spectrum
        - write_stats  : Print statistics on good/bad data

    DATA:
      - errtxt    (str)  : Text message describing error if fail=True          general
      - fail      (bool)  : True=fatal error reading data, False=data read OK.  general
      - file      (str)  : Name of input L1C .nat file                        __init__
      - mph      (dict of str) : Main Product Header                                _read_mph
      - nloc      (int)  : No. of pixel locations selected                    _read_loc
      - scale    (list of float) : Spectral scale factors, size NWNO                  _read_sca
      - stats    (dict of int) : Statistics on good/bad pixels and MDRs             _read_loc
      - mdruse   (list of bool) : True = usable MDR, size TOTAL_MDR                  _read_mdr
      - wno      (list of float) : Wavenumber axis for spectra                        _read_spec
      
      A set of lists each of size nloc:                                   _read_loc
        - day    [int] : Day# (1=1Jan2000)
        - msc    [int] : Milliseconds into day
        - daylin [int] : Day# at start of scan-line
        - msclin [int] : Milliseconds at start of scan-line
        - lin    [int] : Scan line/MDR within file (1:c.760)
        - stp    [i2]  : Mirror  Step# across scan (1:30)
        - pix    [i1]  : Pixel# with FOV (1:4)
        - iof    [ui4] : Byte offset of associated spectrum
        - qal    [3i1] : Quality flags for or each band (0=OK,1=bad)
        - lat    [flt] : Latitude (-90:90)
        - lon    [flt] : Longitude (-180:180)
        - sza    [flt] : Solar Zenith Angle (0:180)
        - zen    [flt] : Satellite Zenith Angle (0:90)
        - saa    [flt] : Solar Azimuth Angle (0:360)
        - azo    [flt] : Satellite Azimuth Angle (0:360)
        - cld    [flt] : % cloud fraction (0:100)
        - lnd    [flt] : % land fraction (0:100)
        - spec   [nwno*flt] : spectrum                                      read_spec
    """

    # Define constants required throughout Class
    NAVC = 6  # Number of AVHRR channels
    NBND = 3  # Number of spectral bands
    NCLS = 7  # Max number of AVHRR clusters
    NPIX = 4  # Number of pixels within FOV (=spectra per scan step)
    NSTP = 30  # Number of scan steps per line across-track
    NWNO = 8461  # Size of spectra in L1C file (8461 = 645-2760@0.25)
    RADFAC = 1.0e7  # Conversion from [W/(m2.sr.m-1)] to [nW/(cm2.sr.cm-1)]
    WNO1 = 645.0  # Start of spectral range [cm-1]
    WNOD = 0.25  # Sampling interval of spectra [cm-1]

    # Define numpy data types required throughout Class
    # Short CDS time (Day0 = 1Jan2000)
    CDS = np.dtype([("day", ">u2"), ("millisec", ">u4")])  #   2b  #   4b  # = 6b total
    # Variable Scale Factor Integer
    VI4 = np.dtype([("ibexp", "b"), ("i4man", ">i4")])  #   1b  #   4b  # = 5b total
    # Generic Record Header
    GRH = np.dtype(
        [
            ("record_class", "b"),  #     1b
            ("instrument_group", "b"),  #     1b
            ("record_subclass", "b"),  #     1b
            ("record_subclass_version", "b"),  #     1b
            ("record_size", ">u4"),  #     4b
            ("record_start_time", CDS),  #     6b
            ("record_stop_time", CDS),  #     6b
        ]
    )  #  = 20b total
    GRH_SIZE = 20  # size of GRH in bytes

    # Measurement Data Record v5 (current version). Differences from v4 = ***
    MDR_A = np.dtype(
        [  #   field   offset
            ("GRH", GRH),  #     20b       0b
            # Generic Quality Indicators
            ("DEGRADED_INST_MDR", "b"),  #      1b      20b
            ("DEGRADED_PROC_MDR", "b"),  #      1b      21b
            # Level 1 Data
            ("GEPSIasiMode", ">i4"),  #      4b      22b
            ("GEPSOPSProcessingMode", ">i4"),  #      4b      26b
            ("GEPSIdConf", "b", 32),  #     32b      30b
            ("GEPSLocIasiAvhrr_IASI", "b", 1200),  #   1200b      62b
            ("GEPSLocIasiAvhrr_IIS", "b", 7500),  #   7500b    1262b
            ("OBT", "b", 180),  #    180b    8762b
            ("OnboardUTC", CDS, NSTP),  #    180b    8942b
            ("GEPSDatIasi", "b", 180),  #    180b    9122b
            ("GisfLinOrigin", ">i4", 2),  #      8b    9302b
            ("GisfColOrigin", ">i4", 2),  #      8b    9310b
            ("GisfPds1", ">i4", 2),  #      8b    9318b
            ("GisfPds2", ">i4", 2),  #      8b    9326b
            ("GisfPds3", ">i4", 2),  #      8b    9334b
            ("GisfPds4", ">i4", 2),  #      8b    9342b
            ("GEPS_CCD", "b", 30),  #     30b    9350b
            ("GEPS_SP", ">i4", 30),  #    120b    9380b
            ("GircImage", ">u2", (NSTP, 64, 64)),  # 245760b    9500b
            ("GQisFlagQual", "b", (NSTP, NPIX, NBND)),  #    360b  255260b  ***
            ("GQisFlagQualDetailed", ">u2", (NSTP, NPIX)),  #    240b  255620b  ***
            ("GQisFlagQualIndex", "b", 5),  #      5b  255860b
            ("GQisFlagQualIndexIIS", "b", 5),  #      5b  255865b
            ("GQisFlagQualIndexLoc", "b", 5),  #      5b  255870b
            ("GQisFlagQualIndexRad", "b", 5),  #      5b  255875b
            ("GQisFlagQualIndexSpect", "b", 5),  #      5b  255880b
            ("GQisSysTecIISQual", ">u4"),  #      4b  255885b
            ("GQisSysTecSondQual", ">u4"),  #      4b  255889b
            ("GGeoSondLoc", ">i4", (NSTP, NPIX, 2)),  #    960b  255893b
            ("GGeoSondAnglesMETOP", ">i4", (NSTP, NPIX, 2)),  #    960b  256853b
            ("GGeoIISAnglesMETOP", ">i4", (NSTP, 25, 2)),  #   6000b  257813b
            ("GGeoSondAnglesSun", ">i4", (NSTP, NPIX, 2)),  #    960b  263813b
            ("GGeoIISAnglesSun", ">i4", (30, 25, 2)),  #   6000b  264773b
            ("GGeoIISLoc", ">i4", (30, 25, 2)),  #   6000b  270773b
            ("EARTH_SATELLITE_DISTANCE", ">u4"),  #      4b  276773b
            # Level 1C Specific Data
            ("IDefSpectDWn1b", VI4),  #      5b  276777b
            ("IDefNsfirst1b", ">i4"),  #      4b  276782b
            ("IDefNslast1b", ">i4"),  #      4b  276786b
        ]
    )
    IOFSPC = 276790  # offset of spectra within MDR record
    # Level 1C spectra inserted here                     2088000b  276790b
    IOFMDRB = 2364790  # offset of MDR_B within MDR record
    MDR_B = np.dtype(
        [  #   field   offset
            ("IDefCovarMatEigenVal1c", "b", 1000),  #   1000b 2364790b
            ("IDEFCcsChannelID", ">i4", 6),  #     24b 2365790b
            ("GCcsRadAnalNbClass", ">i4", (NSTP, NPIX)),  #    480b 2365814b
            ("GCcsRadAnalWgt", VI4, (NSTP, NPIX, NCLS)),  #   4200b 2366294b
            ("GCcsRadAnalY", ">i4", (NSTP, NPIX, NCLS)),  #   3360b 2370494b
            ("GCcsRadAnalZ", ">i4", (NSTP, NPIX, NCLS)),  #   3360b 2373854b
            ("GCcsRadAnalMean", VI4, (NSTP, NPIX, NCLS, NAVC)),  #  25200b 2377213b
            ("GCcsRadAnalStd", VI4, (NSTP, NPIX, NCLS, NAVC)),  #  25200b 2402414b
            ("GCcsImageClassified", "b", 300000),  # 300000b 2427614b
            ("IDefCcsMode", ">u4"),  #      4b 2727614b
            ("GCcsImageClassifiedNbLin", ">i2", NSTP),  #     60b 2727618b
            ("GCcsImageClassifiedNbCol", ">i2", NSTP),  #     60b 2727678b
            ("GCcsImageClassifiedFirstLin", VI4, NSTP),  #    150b 2727738b
            ("GCcsImageClassifiedFirstCol", VI4, NSTP),  #    150b 2727888b
            ("GCcsRadAnalType", "b", (NSTP, NCLS)),  #    210b 2728038b
            ("GlacVarImagIIS", VI4, NSTP),  #    150b 2728248b ***
            ("GlacAvgImagIIS", VI4, NSTP),  #    150b 2728398b ***
            ("GEUMAvhrr1BCldFrac", "b", (NSTP, NPIX)),  #    120b 2728548b ***
            ("GEUMAvhrr1BLandFrac", "b", (NSTP, NPIX)),  #    120b 2728668b ***
            ("GEUMAvhrr1BQual", "b", (NSTP, NPIX)),  #    120b 2728788b ***
        ]
    )
    MDR_SIZE = 2728908  #  Total  2728908b

    def __init__(
        self,
        iasifile,
        chkqal=(True, True, True),
        latlim=None,
        lonlim=None,
        szalim=None,
        zenlim=None,
        saalim=None,
        azilim=None,
        cldlim=None,
        lndlim=None,
        wnolim=(645, 2760),
        avhrr=False,
        mph_only=False,
        loc_only=False,
        bright=False,
    ):
        """Open file and read data.

        Opens file, performs a variety of checks, and read data
        On exit, check self.fail - set True for any error conditions, when
        self.errmsg contains text describing the error.

        Args:
          iasifile (str): Name of .nat file to be read.
          chkqal (list of bool): If True, use band quality flags to screen. Default is
            [True, True, True].
          latlim (array-like): min,max latitude (-90,90) for inclusion. Default is None.
          lonlim (array-like): min,max longitude (-180, 180) for inclusion. Default is None.
          szalim (array-like): min,max solar zenith angle (0,180) for inclusion. Default is None.
          zenlim (array-like): min,max satellite zenith angle (0,90) for inclusion. Default is None.
          saalim (array-like): min,max solar azimuth angle (0,360) for inclusion. Default is None.
          azilim (array-like): min,max satellite azimuth angle (0,360) for inclusion. Default is None.
          cldlim (array-like): min,max %cloud cover (0:100) for inclusion. Default is None.
          lndlim (array-like): min,max %land cover (0:100) for inclusion. Default is None.
          wnolim (array-like): (wnomin, wnomax), Default is (645,2760).
          avhrr (bool): True=include AVHRR cluster analysis in loc structure. Default is None.
          mph_only (bool): True=just extract MPH from file. Default is None.
          loc_only (bool): True=just extract location data, no spectra. Default is None.
          bright (bool): True=convert radiance spectra to brightness temperature. Default is None.

        """

        self.fail = False  # Flag to indicate if an error has occurred
        self.errtxt = ""  # Description of error if self.fail is True
        self.file = iasifile  # Save filename for subsequent get_spec

        with open(iasifile, "rb") as f:  # Open binary file, read-only

            self.mph = self._read_mph(f)  # Read Main Product Header
            if mph_only:
                return

            ptr = self._read_ipr(f)  # Read Internal Pointer Records
            if self.fail:
                return

            self.scale = self._read_sca(f, ptr["giadr_s"])  # Radiance Scale Factors
            if self.fail:
                return

            version, mdruse, mdroff = self._read_mdr(f, ptr["any_mdr"])
            if version is None:
                self.fail = True
                self.errtxt = "No valid Measurement Data Records in file"
                return
            elif version != 5:
                self.fail = True
                self.errtxt = "Measurement Data Records are not v5, instead=v" + str(
                    version
                )
                return
            self.mdruse = (
                mdruse  # Save since this may be useful diagnostic for the user
            )

            self._read_loc(
                f,
                mdruse,
                mdroff,
                avhrr,
                chkqal,  # Read Location data
                latlim,
                lonlim,
                szalim,
                zenlim,
                saalim,
                azilim,
                cldlim,
                lndlim,
            )

            if not loc_only:
                self._read_spec(f, wnolim, bright)  # Read spectral data

        return

    def _read_mph(self, f):
        """Read Main Product Header

        ARGUMENTS
          f obj : file object

        RETURNS
          mph {str} : dictionary of MPH parameter=value pairs

        DESCRIPTION
          The MPH is an arbitrary series of value = parameter text strings.
          These are read and returned as a dictionary mph
          This method assumes the file pointer is initially located at the start of
          the MPH record (which should be the start of the file), and leaves the
          pointer at the end of the MPH record.
        """

        # First read GRH to establish total size of MPH (including GRH part)
        grh = np.fromfile(f, dtype=self.GRH, count=1)[0]
        mphsiz = grh["record_size"] - self.GRH_SIZE  # subtract GRH component
        rec = np.fromfile(f, dtype="b", count=mphsiz)  # Read remaining MPH as bytes
        mphrec = rec.tostring().decode("utf-8")  # convert bytes to characters
        ipt = 0
        mph = {}  # create dictionary
        while ipt < mphsiz:
            ieq = mphrec.find("=", ipt)  # locate next '=' sign
            icr = mphrec.find("\n", ipt)  # locate next end-of-line marker
            param = mphrec[ipt : ieq - 1].strip()  # use left of '=' for key
            value = mphrec[ieq + 1 : icr].strip()  # use right of '=' for value
            mph[param] = value
            ipt = icr + 1  # move to next line
        return mph

    def _read_ipr(self, f):
        """Read Internal Pointer Records

        ARGUMENTS
          f obj : file object

        RETURNS
          ptr {int} : dictionary of IPR record types and associated offsets

        DESCRIPTION
          IPRs are index records which point to the location within the file of
          'target' records of specific types.
            The IPRs follow immediately after the MPH, and the number of IPRs is
          given in the MPH as TOTAL_IPR.
            Usually there are 3 IPR records, sometimes more, and a normal file
          contains IPRs pointing to the following record types
            giadr_q : quality data (expected, but not actually used here)
            giadr_s : scale factors for scaling spectral radiance data
            mdr     : measurement data, ie records corresponding to each scan line
          If there is a loss of data, a specific type of 'lost data' record is
          inserted with an IPR pointing to the start of this section, followed by a
          second MDR IPR pointing to the next record of nominal measurements, etc.
            This method reads the IPRs, checks all the expected IPRs are present,
          and just returns the associated offsets of the first occurrence of each
          type of IPR.
        """

        # Internal Pointer Record structure
        IPR = np.dtype(
            [
                ("grh", self.GRH),  #   20b
                ("target_record_class", "b"),  #    1b 5=GIADR, 8=MDR
                ("target_instrument_group", "b"),  #    1b 8=IASI
                ("target_record_subclass", "b"),  #    1b 0=Qual, 1=ScaleF, 2=MDR-1C
                ("target_record_offset", ">u4"),  #    4b
            ]
        )  # = 27b total

        nipr = int(self.mph["TOTAL_IPR"])
        ptr = {}
        for i in range(nipr):
            ipr = np.fromfile(f, dtype=IPR, count=1)[0]
            target = (
                ipr["target_record_class"],
                ipr["target_instrument_group"],
                ipr["target_record_subclass"],
            )
            offset = ipr["target_record_offset"]
            if target == (5, 8, 0):  # GIADR - quality
                if "giadr_q" not in ptr:
                    ptr["giadr_q"] = offset
            elif target == (5, 8, 1):  # GIADR - scale factors
                if "giadr_s" not in ptr:
                    ptr["giadr_s"] = offset
            elif target == (8, 8, 2):  # MDR
                if "mdr" not in ptr:
                    ptr["mdr"] = offset
            elif target == (8, 13, 1):  # lost data MDR
                if "lost" not in ptr:
                    ptr["lost"] = offset
        # Check all three expected IPR records have been found
        if "giadr_q" not in ptr:
            self.fail = True
            self.errtxt = "No IPR pointer found for GIADR - Quality record"
            return
        if "giadr_s" not in ptr:
            self.fail = True
            self.errtxt = "No IPR pointer found for GIADR - Scale Factors record"
            return
        if "mdr" not in ptr:
            self.fail = True
            self.errtxt = "No IPR pointer found for Measurement Data record"
            return
        # When reading through all MDRs it is necessary to start with the first
        # MDR which may in fact be a 'lost' data MDR, so define as 'any_mdr'
        ptr["any_mdr"] = min(ptr["mdr"], ptr["lost"]) if "lost" in ptr else ptr["mdr"]

        return ptr

    def _read_sca(self, f, offset):
        """Read radiance spectra scale factors

        ARGUMENTS
          f      obj : file object
          offset int : byte offset for scale factors record

        RETURNS
          scale [flt] : numpy 1D array of spectral scale factors

        DESCRIPTION
          Radiance values in the .nat file are stored as 16bit integers, with
          spectrally varying scale factors stored in one of the Global Internal
          Auxiliary Data Records (GIADR - scale factors) to convert these
          to floating values in units of W/(m2.sr.m-1) (NB m-1 rather than cm-1).
            Although the code will handle different specifications of scale factors
          and spectral ranges, it seems that these are always fixed so any variation
          from expected values is treated as an error.
            The scale factors are stored as self.scale(NWNO) - same size as the full
          spectral data.
            Note extra factor RADFAC in method self.spec which changes the units.
        """

        # GIADR Scale Factor Record
        SCA = np.dtype(
            [
                ("grh", self.GRH),  #   20b
                ("IDefScaleSondNbScale", ">i2"),  #    2b No. bands (max=10)
                ("IDefScaleSondNsfirst", ">i2", 10),  #   20b Start Chan# for band
                ("IDefScaleSondNslast", ">i2", 10),  #   20b Last  Chan# for band
                ("IDefScaleSondScaleFactor", ">i2", 10),  #   20b Scale factor 10^i
                ("IDefScaleIISScaleFactor", ">i2"),  #    2b IIS SF as 10^i
            ]
        )  # = 84b Total
        # Locate and load GIADR scale factors record
        f.seek(offset)
        sca = np.fromfile(f, dtype=SCA, count=1)[0]
        # Check number of different ranges is as expected
        if sca["IDefScaleSondNbScale"] != 5:
            self.fail = True
            self.errtxt = "F-READ_L1C: No.band scale factors not set = 5, but =" + str(
                sca["IDefScaleSondNbScale"]
            )
            return None
        # Set up radiance scale factor array, 8421 elements set to zero
        scale = np.zeros((self.NWNO))
        # Scale factors are exponents, eg '7'= multiply spc * 1.0E7 to get radiance
        # in units of W/(m2.sr.m-1) (NB m-1 rather than cm-1).
        # Indices start with 1=0cm-1, 2=0.25cm-1 etc so wno = (i-1)/4 cm-1
        nomlist = (
            2581,
            5921,
            9009,
            9541,
            10721,
            11042,
        )  # 645,1480,2242,2385,2680,2760.25
        ioff = nomlist[0]
        for ibnd in range(5):
            ista = sca["IDefScaleSondNsfirst"][ibnd]
            iend = sca["IDefScaleSondNslast"][ibnd]
            iexp = sca["IDefScaleSondScaleFactor"][ibnd]
            istanom = nomlist[ibnd]
            iendnom = nomlist[ibnd + 1] - 1
            if ista != istanom or iend != iendnom:
                self.fail = True
                self.errtxt = (
                    "F-READ_L1C: Band#"
                    + str(ibnd)
                    + " scale factor not expected range "
                    + str(istanom)
                    + ":"
                    + str(iendnom)
                    + " but "
                    + str(ista)
                    + ":"
                    + str(iend)
                )
                return None
            scale[ista - ioff : iend + 1 - ioff] = 10.0 ** (
                -iexp
            )  # for rad in (W/(m2.sr.m-1))
        return scale

    def _read_mdr(self, f, offset):
        """Check Measurement Data Records and get offsets

        ARGUMENTS
          f      obj : file object
          offset int : byte offset of start of first MDR within file

        RETURNS
          version   int  : MDR format version (5 is expected)
          mdruse   [boo] : True = usable MDR
          mdroff [u4int] : Byte offset for each MDR within file

        DESCRIPTION
          Check MDR is expected size and determine offsets
          Returns lists of size matching 'TOTAL_MDR' in MPH, typically c760.
          Determines version from the first usable MDR, so if initial value
          of None is returned this indicates no usable MDRs were found.
        """

        nmdr = int(self.mph["TOTAL_MDR"])  # No of MDR (scan lines)
        iof = offset
        version = None
        mdruse = []
        mdroff = []
        for imdr in range(nmdr):
            mdroff.append(iof)
            f.seek(iof)
            grh = np.fromfile(f, dtype=self.GRH, count=1)[0]
            record_size = grh["record_size"]
            valid = record_size == self.MDR_SIZE
            if valid and version is None:
                version = grh["record_subclass_version"]
            mdruse.append(valid)
            iof += record_size
        return version, mdruse, mdroff

    def _read_loc(
        self,
        f,
        mdruse,
        mdroff,
        avhrr,
        chkqal,
        latlim,
        lonlim,
        szalim,
        zenlim,
        saalim,
        azilim,
        cldlim,
        lndlim,
    ):
        """Read pixel location data

        ARGUMENTS
          f        obj : file object
          mdruse  [boo]: True = useable MDR
          mdroff  [int]: Byte offset for each MDR
          avhrr    boo : True=include AVHRR cluster analysis in loc structure
          chkqal  (boo): True=use band quality flags to screen
          latlim  (flt): min/max latitude
          lonlim  (flt): min/max longitude
          szalim  (flt): min/max solar zenith angle
          zenlim  (flt): min/max satellite zenith angle
          saalim  (flt): min/max solar azimuth angle
          azilim  (flt): min/max satellite azimuth angle
          cldlim  (flt): min/max cloud percentage
          lndlim  (flt): min/max land percentage

        DESCRIPTION
          Load location data for every pixel/spectrum in file
          Initially creates local lists of maximum possible size but saves
          shorter versions reduced by selection criteria
        """

        # Create local arrays for location data
        NAVC = self.NAVC
        NBND = self.NBND
        NCLS = self.NCLS
        nloc = mdruse.count(True) * self.NSTP * self.NPIX  # Max no expected
        iof = np.zeros(nloc, dtype="u4")  # Byte Offset for start of spectrum
        lin = np.zeros(nloc, dtype="i2")  # Along-track Scan line# (1- c.761)
        stp = np.zeros(nloc, dtype="b")  # Across-track step# (1-30)
        pix = np.zeros(nloc, dtype="b")  # Pixel# (1-4)
        day = np.zeros(nloc, dtype="i4")  # Day# since 1 Jan 2000
        daylin = np.zeros(nloc, dtype="i4")  # Day# at start of MDR/Scan line
        msc = np.zeros(nloc, dtype="i4")  # Millisecs since start of day
        msclin = np.zeros(nloc, dtype="i4")  # Millisecs at start of MDR/Scan line
        qal = np.zeros((nloc, NBND), "b")  # Quality flags for each band (0=OK)
        lat = np.zeros(nloc, "f")  # Latitude [deg N]
        lon = np.zeros(nloc, "f")  # Longitude [deg E]
        sza = np.zeros(nloc, "f")  # Solar Zenith Angle [deg] at sfc (<90=day)
        zen = np.zeros(nloc, "f")  # Zenith angle [deg] at surface (0=nadir)
        saa = np.zeros(nloc, "f")  # Solar Azimuth Angle [deg] at sfc
        azi = np.zeros(nloc, "f")  # Azimuth angle [deg] at surface
        cld = np.zeros(nloc, "f")  # Cloud fraction [%]
        lnd = np.zeros(nloc, "f")  # Land fraction [%]
        if avhrr:
            ncl = np.zeros(nloc, "b")
            pct = np.zeros((nloc, NCLS), "f")
            rav = np.zeros((nloc, NCLS, NAVC), "f")
            rsd = np.zeros((nloc, NCLS, NAVC), "f")
        nbadqal = 0
        nbadgeo = 0
        nmdr = len(mdroff)
        iloc = 0
        for imdr in range(nmdr):  # loop over scan lines/MDRs
            if mdruse[imdr]:
                iofmdr = mdroff[imdr]  # byte offset of start of MDR
                f.seek(iofmdr)
                mdr_a = np.fromfile(f, dtype=self.MDR_A, count=1)[0]
                daymdr = mdr_a["GRH"]["record_start_time"]["day"]
                mscmdr = mdr_a["GRH"]["record_start_time"]["millisec"]
                f.seek(iofmdr + self.IOFMDRB)
                mdr_b = np.fromfile(f, dtype=self.MDR_B, count=1)[0]
                iofspc = iofmdr + self.IOFSPC
                for istp in range(self.NSTP):
                    for ipix in range(self.NPIX):
                        iof[iloc] = (
                            iofspc  # Byte offset for start of associated spectrum
                        )
                        iofspc += 17400  # Spectra 8700 points of I*2 = 17400 bytes each

                        qal[iloc] = mdr_a["GQisFlagQual"][istp, ipix, :]
                        if chkqal:
                            if (
                                chkqal[0] * qal[iloc][0]
                                + chkqal[1] * qal[iloc][1]
                                + chkqal[2] * qal[iloc][2]
                                > 1
                            ):
                                nbadqal += 1
                                continue
                        # Extract location data
                        lat[iloc] = mdr_a["GGeoSondLoc"][istp, ipix, 1] * 1.0e-6
                        lon[iloc] = mdr_a["GGeoSondLoc"][istp, ipix, 0] * 1.0e-6
                        sza[iloc] = mdr_a["GGeoSondAnglesSun"][istp, ipix, 0] * 1.0e-6
                        zen[iloc] = mdr_a["GGeoSondAnglesMETOP"][istp, ipix, 0] * 1.0e-6
                        saa[iloc] = mdr_a["GGeoSondAnglesSun"][istp, ipix, 1] * 1.0e-6
                        azi[iloc] = mdr_a["GGeoSondAnglesMETOP"][istp, ipix, 1] * 1.0e-6
                        cld[iloc] = mdr_b["GEUMAvhrr1BCldFrac"][istp, ipix]
                        lnd[iloc] = mdr_b["GEUMAvhrr1BLandFrac"][istp, ipix]
                        # Bad location set to -2147483648, converted to lat=-2147.48 deg
                        if lat[iloc] < -100.0:
                            nbadgeo += 1
                            continue
                        if (
                            self._outside(lat[iloc], latlim)
                            or self._outside(lon[iloc], lonlim)
                            or self._outside(sza[iloc], szalim)
                            or self._outside(zen[iloc], zenlim)
                            or self._outside(saa[iloc], saalim)
                            or self._outside(azi[iloc], azilim)
                            or self._outside(cld[iloc], cldlim)
                            or self._outside(lnd[iloc], lndlim)
                        ):
                            continue
                        # At this point the location will be included
                        daylin[iloc] = daymdr
                        msclin[iloc] = mscmdr
                        lin[iloc] = imdr + 1
                        stp[iloc] = istp + 1
                        pix[iloc] = ipix + 1
                        day[iloc] = mdr_a["OnboardUTC"][istp]["day"]
                        msc[iloc] = mdr_a["OnboardUTC"][istp]["millisec"]
                        if avhrr:
                            ncl[iloc] = mdr_b["GCcsRadAnalNbClass"][istp, ipix]
                            pct[iloc] = (
                                100.0
                                * mdr_b["GCcsRadAnalWgt"][istp, ipix, :]["i4man"]
                                * 0.1 ** mdr_b["GCcsRadAnalWgt"][istp, ipix, :]["ibexp"]
                            )
                            rav[iloc] = (
                                100.0
                                * mdr_b["GCcsRadAnalMean"][istp, ipix, :, :]["i4man"]
                                * 0.1
                                ** mdr_b["GCcsRadAnalMean"][istp, ipix, :, :]["ibexp"]
                            )
                            rsd[iloc] = (
                                100.0
                                * mdr_b["GCcsRadAnalStd"][istp, ipix, :, :]["i4man"]
                                * 0.1
                                ** mdr_b["GCcsRadAnalStd"][istp, ipix, :, :]["ibexp"]
                            )

                        iloc += 1
        self.nloc = iloc
        self.iof = iof[0:iloc]
        self.lin = lin[0:iloc]
        self.stp = stp[0:iloc]
        self.pix = pix[0:iloc]
        self.day = day[0:iloc]
        self.daylin = daylin[0:iloc]
        self.msc = msc[0:iloc]
        self.msclin = msclin[0:iloc]
        self.qal = qal[0:iloc]
        self.lat = lat[0:iloc]
        self.lon = lon[0:iloc]
        self.sza = sza[0:iloc]
        self.zen = zen[0:iloc]
        self.saa = saa[0:iloc]
        self.azi = azi[0:iloc]
        self.cld = cld[0:iloc]
        self.lnd = lnd[0:iloc]
        self.stats = {
            "tot_mdr": nmdr,  # Tot.No of Meas.Data records
            "bad_mdr": mdruse.count(False),  # No. bad/lost MDRs
            "use_mdr": mdruse.count(True),  # No. used MDRs
            "tot_loc": nloc,  # Tot.No pixels in valid MDRs
            "use_loc": iloc,  # No.used  pixels in valid MDRs
            "qal_loc": nbadqal,  # No.excluded for bad qual flg
            "geo_loc": nbadgeo,
        }  # No.excluded for bad geo. data

    def _outside(self, value, limits):
        """True if value outside limits

        ARGUMENTS
          value   flt  : Value to be tested
          limits (flt) : min,max values for test, or None if not specified

        RETURNS
           boo : True if outside min,max or None; False if inside or unspecified

        DESCRIPTION
           Normally called with min<max, but if max>min this returns True if
           value is *outside* range min:max (usefule for longitude wrap-around)
           If called with limits=None this will always return False
        """

        if limits is None:  # no limits specified
            return False
        else:
            if limits[1] >= limits[0]:  # normal min < max
                return not limits[0] <= value <= limits[1]
            else:  # abnormal min > max
                return limits[1] <= value <= limits[0]

    def _btspec(self, wno, rad):
        """Convert radiance spectrum to brightness temperature

        ARGUMENTS
          wno [flt] : Wavenumber [cm-1] values
          rad [flt] : Radiance spectrum

        RETURNS
          [flt] : spectrum of brightness temperatures [K]

        DESCRIPTION
          Use the inverse of the Planck function to convert radiance to
          brightness temperature.
            Also include RADFAC here since C2 depends on radiance units
        """

        # Radiation constants for brightness temperature conversion
        C1 = 1.4387686  # =hc/k [K/cm-1]
        # RADFAC value will be 1e7 if units match C2
        C2 = 1.191044e-3 * self.RADFAC / 1.0e7  # [nW/(cm2.sr.(cm-1))]
        w3 = wno**3
        return C1 * wno / np.log(1.0 + C2 * w3 / rad)

    def _read_spec(self, f, wnolim, bright):
        """Read in spectra at selected locations

        ARGUMENTS
          f       obj : file object
          wnolim (flt): min, max wavenumbers for spectra
          bright  boo : True = output as brightness temperature spectra

        DEFINES
          self.spec, self.wno

        DESCRIPTION
          The byte offsets of spectra within the file are saved as self.iof, so this
          method can simply move to the correct location.
            L1C spectral data are stored as big-endian 16-byte signed integers, '>i2'
          which are multiplied by self.scale to covert to radiances in
          W/(m2.sr.m-1) (NB m-1 rather than cm-1), and then by RADFAC to convert to
          nW/(cm2.sr.cm-1)
            Also sets self.wno containing wavenumber axis for spectra
        """

        # identify start,end points - would be 0,8460 for full range
        i = round((wnolim[0] - self.WNO1) / self.WNOD)
        j = round((wnolim[1] - self.WNO1) / self.WNOD)
        if i < 0:
            print("F-read_l1c: requested lower wavenumber < 645cm-1")
            exit()
        if j > self.NWNO:
            print("F-read_l1c: requested upper wavenumber > 2760 cm-1")
            exit()
        n = j - i + 1
        w = wnolim[0] + np.arange(n) * self.WNOD
        self.wno = w
        self.spec = [np.array(n)] * self.nloc
        for iloc in range(self.nloc):
            iof = self.iof[iloc]
            f.seek(iof + 2 * i)
            rad = (
                np.fromfile(f, dtype=">i2", count=n)
                * self.scale[i : j + 1]
                * self.RADFAC
            )
            if bright:
                self.spec[iloc] = self._btspec(w, rad)
            else:
                self.spec[iloc] = rad
        if len(self.spec) > 0:
            self.spec = np.vstack(self.spec)

    def get_spec(self, iloc, wnolim=(645, 2760), bright=False):
        """Read in spectra at selected locations.

        The byte offsets of spectra within the file are saved in the self.iof list
        so this method can simply move to the correct location.
        Assumes self.file contains name of required file which is reopened.
        Otherwise works as read_spec method.

        Args:
          iloc (int): Index of pixel for spectrum
          wnolim (tuple): (Min,Max) wavenumber limits [cm-1] for extracted spectra.
            Default is (645.2760)
          bbt (bool): If True, return spec as brightness temps, if False, return
            radiances. Default is False.

        Returns:
          rad (array_like): radiances or brightness temperatures
          wno (array_like): wavenumbers
          
          

        """
        # identify start,end points - would be 0,8460 for full range
        i = round((wnolim[0] - self.WNO1) / self.WNOD)
        j = round((wnolim[1] - self.WNO1) / self.WNOD)
        if i < 0:
            print("F-read_l1c: requested lower wavenumber < 645cm-1")
            exit()
        if j > self.NWNO:
            print("F-read_l1c: requested upper wavenumber > 2760 cm-1")
            exit()
        n = j - i + 1

        wno = wnolim[0] + np.arange(n) * self.WNOD
        with open(self.file, "rb") as f:  # Open binary file, read-only
            iof = self.iof[iloc]
            f.seek(iof + 2 * i)
            rad = (
                np.fromfile(f, dtype=">i2", count=n)
                * self.scale[i : j + 1]
                * self.RADFAC
            )
            if bright:
                return self._btspec(wno, rad), wno
            else:
                return rad, wno

    def write_stats(self):
        """Write summary statistics to terminal."""

        print(
            "{:6}".format(self.stats["tot_mdr"])
            + "          = Tot.No of Meas.Data records"
        )
        print(
            "{:6}".format(self.stats["use_mdr"])
            + " ({:5.1f}%)".format(100 * self.stats["use_mdr"] / self.stats["tot_mdr"])
            + " = No. used MDRs"
        )
        print(
            "{:6}".format(self.stats["bad_mdr"])
            + " ({:5.1f}%)".format(100 * self.stats["bad_mdr"] / self.stats["tot_mdr"])
            + " = No. bad/lost MDRs"
        )
        print(
            "{:6}".format(self.stats["tot_loc"])
            + "          = Tot.No pixels in valid MDRs"
        )
        print(
            "{:6}".format(self.stats["use_loc"])
            + " ({:5.1f}%)".format(100 * self.stats["use_loc"] / self.stats["tot_loc"])
            + " = No used pixels in valid MDRs"
        )
        print(
            "{:6}".format(self.stats["qal_loc"])
            + " ({:5.1f}%)".format(100 * self.stats["qal_loc"] / self.stats["tot_loc"])
            + " = No.excluded pixels for bad quality flags"
        )
        print(
            "{:6}".format(self.stats["geo_loc"])
            + " ({:5.1f}%)".format(100 * self.stats["geo_loc"] / self.stats["tot_loc"])
            + " = No.excluded pixels for bad geo. data"
        )
