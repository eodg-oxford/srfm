"""Reader for IASI L2 data.

- Name: read_iasi_l2.py
- Parent package: srfm
- Author: Anu Dudhia
- Contributors: Antonin Knizek
- Date: 25 July 2025

"""

import numpy as np

class Iasi_L2: 
  """ Self-contained Python module to read IASI v6 L2 product

  Version
    05AUG25 AK Added solar and satellite azimuth angle readings.
    04AUG25 AD Corrections to reading Forli data
    15MAR25 AD Modify _loclim to exclude invalid fg_tem profiles
    11FEB25 AD Rename from read_iasi_l2.py to iasi_l2_class.py
    13JAN25 AD Correction to loclim implementation.
    11OCT24 AD Renamed from Iasi_L2.
    05SEP24 AD In _read_ipr only take first ipr_mdr if multiple in file.
    28AUG24 AD Renamed from read_iasi_l2. Add loclim parameter.
    04MAR24 AD Converted from IDL

  Parameters (for initialisation)
    l2fil     str  : Name of L2 file 
   (The following are optional, defaulting to False or None)
    mph_only  boo  : True=Only read Main Product Header 
    forli     boo  : True=include Forli/Brescia data 
    latlim   [flt] : [South,North] latitude limits  
    lonlim   [flt] : [West,East] longitude limits 
    loclim    obj  : Lists of .lat[] and .lon[] values of reqd locations
    kwargs    {}   : Dictionary of mdr fields for each fov and required values
                     eg: flg_retcheck=0 will only return 'good' retrievals
  
  Returns
    structure, containing
      .errtxt str  : Text message describing error if fail is True
      .fail   boo  : True=fatal error reading data, False=data read OK
      .mph   {dic} : Main Product Header - see _read_mph
      .grd   {dic} : Profile grid sizes and levels
      .nloc   int  : Total number of pixel locations in output
    then L2 data as lists, each of size nloc. See _inivar for definitions.

  Description
    Read Eumetsat v6 .nat format IASI L2 data file for one orbit.
    The profile data are returned as a set of lists, with each element
      corresponding to a pixel location. In addition there are two
      dictionaries containing the main product header and the profile grids.
    Default is to read entire orbit but exclude the forli/brescia parts.
    Set mph_only=True to just read the main product header.
    Set forli=True to include the forli (co,hno3,o3) and brescia (so2) 
      products. Note that forli o3 is in addition to regular o3 profile.
    Set latlim and lonlim to restrict lat and/or lon limits. The first 
      lonlim value is interpreted as the western boundary so, for example
      lonlim=[170,-170] would span 20 deg longitude across the 180 meridian

  Reference
    PDF_TEN_980760-EPS-IASI-L2.pdf

  Usage
    Assuming 'l2fil.nat' is the name of an IASI L2 .nat file, 
    to select profiles within a lat/lon box

      from iasi_l2_class import Iasi_L2

      l2 = Iasi_L2 ( l2fil.nat, latlim=[30,60], lonlim=[0,90], \
                              flg_retcheck=0 )
      if l2.fail:             # check for error during read
        print(l2.errtxt)
        exit()
      
      print(l2.nloc)               # total no. of locations
      print(l2.grd['pre_t'])       # press.levels for T profile
      print(l2.lat[0], l2.lon[0])  # lat,lon of first profile
      print(l2.tem[0])             # first temperature profile


  Private Methods
    __init__    : Initialise
    _check_mph  : Perform checks on file structure as described in MPH
    _get_grd    : Convert GIADR contents into dictionary
    _get_mdroff : Determine offsets for each MDR
    _inivar     : Initialise L2 list variables
    _latlim     : Limit pixels to within latitude range
    _lonlim     : Limit pixels to within longitude range
    _loclim     : Limit pixels to (lat,lon) locations
    _read_err   : Read error covariance part of MDR record
    _read_giadr : Read GIADR 
    _read_ipr   : Read Internal Pointer Records
    _read_mph   : Read Main Product Header
    _set_loc    : Set output location data

  """
  # List of constants and fixed data structures

  # Scale factor for converting L2 data to physical units in output profiles
  TSCALE = 0.01             # Temp values are K*100, convert to K
  QSCALE = 0.1*28.96/18.0   # H2O values are kg/kg*1e7, convert to ppmv 
  OSCALE = 0.01*28.96/48.0  # O3  values are kg/kg*1e8, convert to ppmv 

  # Array sizes given in the 8.6 Parameter Table (p69 EUM/OPS-ES/MAN/04/0033)
  NFOV      = 120 # No. pixels in scan line
  NLT       = 101 # No press levels for temperature  (increased from 90 in v5)
  NLQ       = 101 # No press levels for water vapour (increased from 90 in v5)
  NLO       = 101 # No press levels for ozone        (increased from  4 in v5)
  NEW       = 12  # No wavelengths fo emissivities
  NL_CO     = 19  # No. of partial layers for CO
  NL_HNO3   = 19  # No. of partial layers for HNO3
  NL_O3     = 40  # No. of partial layers for O3
  NL_SO2    = 5   # No. of estimate SO2 plume heights
  NPCT      = 28  # No. of PCs for temperature in error data 
  NPCW      = 18  # No. of PCs for water-vapour in error data 
  NPCO      = 10  # No. of PCs for ozone in error data 
  NEVA_CO   = 10  # Max No EVs for CO sensitivity matrix
  NEVE_CO   = 190 # Max No elements in the CO sensitivity matrix EVs
  NEVA_HNO3 = 10  # Max No EVs for HNO3 sensitivity matrix
  NEVE_HNO3 = 190 # Max No elements in the HNO3 sensitivity matrix EVs
  NEVA_O3   = 10  # Max No EVs for O3 sensitivity matrix
  NEVE_O3   = 190 # Max No elements in the O3 sensitivity matrix EVs
  NERRT     = 406 # No. elements in temperature error record
  NERRW     = 171 # No. elements in water-vapour error record
  NERRO     = 55  # No. elements in ozone error record
  NERR      = 30  # No. err data records in scan line
  CO_NBR    = 50  # No. CO retrievals in scan line
  HNO3_NBR  = 50  # No. HNO3 retrievals in scan line
  O3_NBR    = 50  # No. O3 retrievals in scan line
  NCLD      = 3   # No. diffent cloud types per fov (local definition)

  # Scale factors for FORLI data
  SCALEF_CO   = 1.0e-13
  SCALEF_HNO3 = 1.0e-11
  SCALEF_O3   = 1.0e-14

  # Bad data values
  BADI2 = np.int16(32767)               # 2^15 - 1,   16 '1's bit pattern
  BADU2 = np.uint16(65535)              # 2^16 - 1,   16 '1's bit pattern
  BADU4 = np.uint32(4294967295)         # 2^32 - 1,   32 '1's bit pattern

  # Data structure definitions
  # Compound data types (section 2.3, p23 of EPS.GGS.SPE.96167)
  # vu-integer2 (1 byte exp, 2 byte mant.)
  VUI2 = np.dtype ( [ ('ibexp', 'b'), ('i2man', '>u2') ] )

  # v-integer4 (1 byte exp, 4 byte mant.)
  VI4 = np.dtype ( [ ('ibexp', 'b'), ('i4man', '>i4') ] )

  # Short CDS time (Day0 = 1Jan2000)
  CDS = np.dtype ( [ ('day', '>u2'),         # 2b   
                     ('millisec', '>u4') ] ) # 4b = 6b total

  # Generic Record Header (p27 of EPS.GGS.SPE.96167)
  GRH = np.dtype ( [ ('record_class', 'b'),            # 1b
                     ('instrument_group', 'b'),        # 1b
                     ('record_subclass', 'b'),         # 1b 
                     ('record_subclass_version', 'b'), # 1b
                     ('record_size', '>u4'),           # 4b
                     ('record_start_time', CDS),       # 6b
                     ('record_stop_time', CDS)   ] )   # 6b = 20b total
  GRH_SIZE = 20              # size = 20 bytes, first 20 bytes of every record

  # Fixed length part of Measurement Data Record (p61 EUM/OPS-EPS/MAN/04/0033)
  MDR = np.dtype ( [( 'grh', GRH ),                  #     20b     20
    # Generic quality indicators
               ('degraded_inst_mdr',   'b'),              #      1b     21
               ('degraded_proc_mdr',   'b'),              #      1b     22
    # First guess profiles
      ('fg_atmospheric_temperature', '>u2', (NFOV,NLT)),  #  24240b  24262
     ('fg_atmospheric_water_vapour', '>u4', (NFOV,NLQ)),  #  48480b  72742
            ('fg_atmospheric_ozone', '>u2', (NFOV,NLO)),  #  24240b  96982
          ('fg_surface_temperature', '>u2', NFOV),        #    240b  97222
   ('fg_qi_atmospheric_temperature',  'u1', NFOV),        #    120b  97342
  ('fg_qi_atmospheric_water_vapour',  'u1', NFOV),        #    120b  97462
         ('fg_qi_atmospheric_ozone',  'u1', NFOV),        #    120b  97582
       ('fg_qi_surface_temperature',  'u1', NFOV),        #    120b  97702
    # Measurement data
         ('atmospheric_temperature', '>u2', (NFOV,NLT)),  #  24240b 121942
        ('atmospheric_water_vapour', '>u4', (NFOV,NLQ)),  #  48480b 170422
               ('atmospheric_ozone', '>u2', (NFOV,NLO)),  #  24240b 194662
             ('surface_temperature', '>u2', NFOV),        #    240b 194902
         ('integrated_water_vapour', '>u2', NFOV),        #    240b 195142
                ('integrated_ozone', '>u2', NFOV),        #    240b 195382
                  ('integrated_n2o', '>u2', NFOV),        #    240b 195622
                   ('integrated_co', '>u2', NFOV),        #    240b 195862
                  ('integrated_ch4', '>u2', NFOV),        #    240b 196102
                  ('integrated_co2', '>u2', NFOV),        #    240b 196342
              ('surface_emissivity', '>u2', (NFOV,NEW)),  #   2880b 199222 
         ('number_cloud_formations',  'u1', NFOV),        #    120b 199342
          ('fractional_cloud_cover', '>u2', (NFOV,NCLD)),   #    720b 200062
           ('cloud_top_temperature', '>u2', (NFOV,NCLD)),   #    720b 200782
              ('cloud_top_pressure', '>u4', (NFOV,NCLD)),   #   1440b 202222
                     ('cloud_phase',   'b', (NFOV,NCLD)),   #    360b 202582
                ('surface_pressure', '>u4', NFOV),        #    480b 203062
    # Instrument 
                 ('instrument_mode',   'b'),              #      1b 203063
    # Navigation data at scan line
             ('spacecraft_altitude', '>u4'),              #      4b 203067
    # Navigation data at IFOV
                ('angular_relation', '>i2', (NFOV,4)),    #    960b 204027
                  ('earth_location', '>i4', (NFOV,2)),    #    960b 204987
    # Processing and quality flags
                    ('flg_amsubad',   'b', (NFOV)),      #    120b 205107
                   ('flg_avhrrbad',   'b', (NFOV)),      #    120b 205227
                     ('flg_cldfrm',   'b', (NFOV)),      #    120b 205347
                     ('flg_cldnes',   'b', (NFOV)),      #    120b 205467
                     ('flg_cldtst', '>u2', (NFOV)),      #    240b 205707
                     ('flg_daynit',   'b', (NFOV)),      #    120b 205827
                    ('flg_dustcld',   'b', (NFOV)),      #    120b 205947
                    ('flg_fgcheck', '>u2', (NFOV)),      #    240b 206187
                    ('flg_iasibad',   'b', (NFOV)),      #    120b 206307
                     ('flg_initia',   'b', (NFOV)),      #    120b 206427
                     ('flg_itconv',   'b', (NFOV)),      #    120b 206547
                     ('flg_lansea',   'b', (NFOV)),      #    120b 206667
                     ('flg_mhsbad',   'b', (NFOV)),      #    120b 206787
                      ('flg_numit',   'b', (NFOV)),      #    120b 206907
                     ('flg_nwpbad',   'b', (NFOV)),      #    120b 207027
                  ('flg_physcheck',   'b', (NFOV)),      #    120b 207147
                   ('flg_retcheck', '>u2', (NFOV)),      #    240b 207387
                     ('flg_satman',   'b', (NFOV)),      #    120b 207507
                    ('flg_sunglnt',   'b', (NFOV)),      #    120b 207627
                     ('flg_thicir',   'b', (NFOV)) ])    #    120b 207747

  # EXECUTABLE CODE-------------------------------------------------------------

  def __init__ ( self, l2fil, mph_only=False, forli=False, 
                 latlim=None, lonlim=None, loclim=None, **kwargs ):
    """ Initialise 

    Parameters (for initialisation)
      l2fil     str  : Name of L2 file 
    The following are optional, defaulting to False or None
      mph_only  boo  : True=Only read Main Product Header 
      forli     boo  : True=include Forli/Brescia data 
      latlim   [flt] : [South,North] latitude limits  
      lonlim   [flt] : [West,East] longitude limits 
      loclim    obj  : .lat[] and .lon[] locations
    """

    # Define error indicators and set to 'OK' status
    self.fail = False
    self.errtxt = None

    # Create flags for points which match loclim
    if loclim:
      locmsc = []
      for loc in loclim:
        if loc['msclin'] not in locmsc: locmsc.append(loc['msclin'])
      self.useloc = np.zeros(len(loclim),dtype='bool')

    # Open L2 .nat file 
    with open ( l2fil, 'rb' ) as f:     # open binary file, readonly

      self.mph = self._read_mph ( f )   # Read Main Product Header
      if mph_only: return

      self._check_mph()                 # Check file structure
      if self.fail: return

      ipr_giadr, ipr_mdr = self._read_ipr ( f ) # Read internal pointer records
      if self.fail: return

      # Extract grid levels for profiles from GIADR
      giadr = self._read_giadr ( f, ipr_giadr )
      self.grd = self._get_grd ( giadr, forli )   # save as dictionary 

      # No of Measurement Data Records = No.scan lines 
      # (expect around 770 per orbit)
      total_mdr = int(self.mph['TOTAL_MDR'])
      # Get byte offsets for each MDR
      mdroff = self._get_mdroff(f, ipr_mdr, total_mdr) 
      if self.fail: return

      self._inivar ( forli )             # Initialise data lists for output

      for imdr in range(total_mdr):      # Loop over scan lines
        if mdroff[imdr] is None: continue
        f.seek ( mdroff[imdr] )          # Byte offset for start of MDR
        # Read fixed part of MDR
        mdr = np.fromfile ( f, dtype=self.MDR, count=1)[0]  

        # usefov is size 120 array of booleans, True=use data from this location
        usefov = np.ones(self.NFOV,dtype=bool)
        if latlim: usefov &= self._latlim ( mdr, latlim )
        if lonlim: usefov &= self._lonlim ( mdr, lonlim )
        if loclim: usefov &= self._loclim ( mdr, loclim, locmsc )
        if not any ( usefov ) : continue  # nothing with lat,lon, so skip to next MDR
        for keyword in kwargs:
          usefov &= mdr[keyword] == kwargs[keyword]

        # surface altitude sfz might be useful, although need to read into the 
        # variable part of the MDR to extract it
        err = self._read_err ( f )  # read error data, but don't pass to output
        mdr_sfz = np.fromfile ( f, dtype='>i2', count=self.NFOV ) # 240b

        # idxfov is list of indices in range 0:119 of locations in MDR to use
        idxfov = np.arange(self.NFOV)[ np.where ( usefov ) ]
        self._set_loc ( mdr, imdr, idxfov )    # extract time/location data
        self._set_flags ( mdr, idxfov )        # extract flags
        self._set_profs ( mdr, idxfov )        # extract profile data
        self._set_cols  ( mdr, idxfov )        # extract column data
        self._set_sfc ( mdr, mdr_sfz, idxfov ) # extract surface data 
        self._set_cloud ( mdr, idxfov )        # extract cloud data
 
        if not forli: continue

        # Continue reading Forli/Brescia parts of MDR - variable sizes
        forli = self._read_forli ( f, self.NL_CO, self.NEVE_CO, self.NEVA_CO ) 
        self._set_forli ( forli, idxfov, self.SCALEF_CO, 'co' )

        forli = self._read_forli ( f, self.NL_HNO3, self.NEVE_HNO3, self.NEVA_HNO3 ) 
        self._set_forli ( forli, idxfov, self.SCALEF_HNO3, 'hno3' )

        forli = self._read_forli ( f, self.NL_O3, self.NEVE_O3, self.NEVA_O3 ) 
        self._set_forli ( forli, idxfov, self.SCALEF_O3, 'o3' )

        brescia = self._read_brescia ( f, self.NL_SO2 )
        self._set_brescia ( brescia, idxfov )

    self.nloc = len ( self.pix )  
    return 


  def _check_mph ( self ):
    """ Perform checks on file structure as described in MPH

    Description
      Checks file structure as defined in the Main Product Header and
      sets self.fail = True and sets self.errmsg if anything unexpected
      is found.
      
    Returns
      None
    """

    # Check data version 
    if self.mph['PROCESSOR_MAJOR_VERSION'] != '6':
      self.fail = True
      self.errtxt = 'F-Iasi_L2_Class: Expected PROCESSOR_MAJOR_VERSION=6' + \
                    ' but found =' + self.mph['PROCESSOR_MAJOR_VERSION'] 

    # Expect 1 Global Internal Auxiliary Data Record: Pressure levels etc
    elif self.mph['TOTAL_GIADR'] != '1':
      self.fail = True
      self.errtxt = 'F-Iasi_L2_Class: Expected TOTAL_GIADR=1' + \
                    ' but found =' + self.mph['TOTAL_GIADR']

    # Expect 2 Internal Pointer Records, 1 GIADR and 1 MDR:
#    elif self.mph['TOTAL_IPR'] != '2':
#      self.fail = True
#      self.errtxt = 'F-Iasi_L2_Class: Expected TOTAL_IPT=2' + \
#                    ' but found =' + self.mph['TOTAL_IPR']
#    return



  def _get_grd ( self, giadr, forli ):
    """ Convert GIADR contents into dictionary

    Parameters
      giadr  obj : GIADR record contents
      forli  boo : True = include Forli/Brescia grids

    Returns
             {}  : Dictionary of data grid sizes and levels
    """

    grd = { 'nlt':giadr['num_pressure_levels_temp'], 
          'pre_t':giadr['pressure_levels_temp'] * 1.0e-4,          # hPa
            'nlq':giadr['num_pressure_levels_humidity'], 
          'pre_q':giadr['pressure_levels_humidity'] * 1.0e-4,      # hPa
            'nlo':giadr['num_pressure_levels_ozone'], 
          'pre_o':giadr['pressure_levels_ozone'] * 1.0e-4,         # hPa
            'new':giadr['num_surface_emissivity_wavelengths'], 
          'wno_s':1.0e8/giadr['surface_emissivity_wavelengths'] }  # cm-1
    if forli:
      grd = { **grd, 
          'nl_co':giadr['forli_num_layers_co'], 
         'hgt_co':giadr['forli_layer_heights_co']*1.0,             # m
        'nl_hno3':giadr['forli_num_layers_hno3'], 
       'hgt_hno3':giadr['forli_layer_heights_hno3']*1.0,           # m
          'nl_o3':giadr['forli_num_layers_o3'], 
         'hgt_o3':giadr['forli_layer_heights_o3']*1.0,             # m
         'nl_so2':giadr['brescia_num_altitudes_so2'], 
        'hgt_so2':giadr['brescia_altitudes_so2']*1.0       }       # m

      # lowest level of Forli grids always seems to be set to BADU2
      # interpret this as -1 m.
      if giadr['forli_layer_heights_co'][0]   == self.BADU2: 
        grd['hgt_co'][0]   = -1.0
      if giadr['forli_layer_heights_hno3'][0] == self.BADU2: 
        grd['hgt_hno3'][0] = -1.0
      if giadr['forli_layer_heights_o3'][0]   == self.BADU2: 
        grd['hgt_o3'][0]   = -1.0
    return grd


  def _get_mdroff ( self, f, offset, nmdr ):
    """ Determine offsets for each MDR

    Parameters
      f       obj : file object
      offset  int : Byte offset of first MDR
      nmdr    int : No. MDRs to find

    Returns
      [int] : list of offsets for each MDR

    Description
      Starting from the first MDR, reads the record size from the GRH 
      header of each MDR and adds this to the previous offset
      Skips any dummy MDR records
    """

    mdroff = []
    #print('nmdr=', nmdr)
    for imdr in range(nmdr): 
      f.seek(offset)
      #print('offset=', offset, imdr)
      grh = np.fromfile(f,dtype=self.GRH,count=1)[0]  
      #print(grh['instrument_group'],grh['record_size'])
      if grh['instrument_group'] == 15:      # standard MDR, 15=IASI L2
        mdroff.append(offset)
      elif grh['instrument_group'] == 13:    # missing record
        print('W-Iasi_L2_Class: Ignoring missing MDR#' + str(imdr+1))
        mdroff.append(None)
      else: 
        self.fail = True
        self.errtxt = 'F-Iasi_L2_Class: Unexpected MDR.INSTRUMENT_GROUP =' + \
               str(grh['instrument_group']) + ', expected: 15 or 13'
        return None
      offset += grh['record_size']
    return mdroff


  def _inivar ( self, forli ): 
    """ Initialise L2 list variables

    Parameters
      forli     boo  : True=include Forli/Brescia data 

    Returns
      None
    """
    self.fov = []     # FOV# within scan line 1, 2, ... 120
    self.lin = []     # Scan line number (equivalent to MDR#) 1, 2, ... c770 
    self.stp = []     # Cross-scan step number 1, 2, ... 30
    self.pix = []     # Pixel number at each step 1, 2, 3, 4
    self.daylin = []  # Julian day at start of scan line
    self.msclin = []  # Millisec within day at start of scan line
    self.lat = []     # Latitude [degN] (-90:90_
    self.lon = []     # Longitude [degE] (-180:180)
    self.sza = []     # Solar Zenith Angle [deg] (0:180) 0=sun at zenith
    self.zen = []     # Satellite Zenith Angle [deg] (0:c.60) 0=sat at zenith
    self.saa = []     # Solar Azimuth Angle [deg] (0:360)
    self.azi = []     # Satellite Azimuth Angle [deg] (0:360) 
    self.flg_amsubad   = [] # (byte) Availability and quality of AMSU measurements
    self.flg_avhrrbad  = [] # (byte) Availability and quality of AVHRR measurements
    self.flg_cldfrm    = [] # (byte) Origin of characterisation of the cloud formations
    self.flg_cldnes    = [] # (byte) Cloudiness assessment summary
    self.flg_cldtst    = [] # (i2)   Details of cloud tests exectued and their results
    self.flg_daynit    = [] # (byte) Discrimination between day(0) and night(1)
    self.flg_dustcld   = [] # (byte) Indicates presence of dust clouds in the IFOV
    self.flg_fgcheck   = [] # (i2)   Check geophys.params from the first guess within bounds
    self.flg_iasibad   = [] # (byte) Availability and quality of IASI L1 measurements
    self.flg_initia    = [] # (byte) Indicates measurements used in the first guess retrieval
    self.flg_itconv    = [] # (byte) Convergence and acceptance of the OEM result
    self.flg_lansea    = [] # (byte) Specifies surface type
    self.flg_mhsbad    = [] # (byte) Availability and quality of MHS measurements
    self.flg_numit     = [] # (byte) Number of iterations in the OEM
    self.flg_nwpbad    = [] # (byte) Availability and quality of NWP data
    self.flg_physcheck = [] # (byte) Potential correct. for superadiabatic/supersat. conditions
    self.flg_retcheck  = [] # (i2)   Check that geophys.params from the OEM within bounds
    self.flg_satman    = [] # (byte) Indication of satellite manouevre
    self.flg_sunglnt   = [] # (byte) Identification of sun glint
    self.flg_thicir    = [] # (byte) Thin cirrus cloud test

    self.tem_a = []   # First guess temperature profile [K]
    self.h2o_a = []   # First guess water vapour profile [ppmv]
    self.o3_a  = []   # First guess ozone profile [ppmv]
    self.sft_a = []   # First guess surface skin temperature [K]
    self.tem_aq = []  # Quality indicator for tem_a
    self.h2o_aq = []  # Quality indicator for h2o_a
    self.o3_aq  = []  # Quality indicator for o3_a
    self.sft_aq = []  # Quality indicator for sft_a
    self.sfz = []     # Surface elevation [m] (from Forli data)
    self.tem = []     # Retrieved temperature profile [K]
    self.h2o = []     # Retrieved water vapour profile [ppmv] 
    self.o3  = []     # Retrieved ozone profile [ppmv]
    self.sft = []     # Retrieved surface temperature [K]
    self.h2o_col = [] # Integrated water vapour column [kg/m2]
    self.o3_col  = [] # Integrated ozone column [kg/m2]
    self.n2o_col = [] # Integrated nitrous oxide column [kg/m2]
    self.co_col  = [] # Integrated carbon monoxide column [kg/m2]
    self.ch4_col = [] # Integrated methane column [kg/m2]
    self.co2_col = [] # Integrated carbon dioxide column [kg/m2]
    self.sfe = []     # Surface emissivity at various wavelengths
    self.ncld = []    # Number of cloud formations
    self.cld = []     # Fractional cloud cover [%]
    self.ctt = []     # Cloud top temperature [K]
    self.ctp = []     # Cloud top pressure [hPa]
    self.cph = []     # Cloud phase (0=nocloud,1=liq,2=ice,3=mix,4=undefined)
    self.sfp = []     # Surface pressure [hPa]

    if forli: 
      # For CO
      self.co_qflag = []       # General retrieval quality flag
      self.co_bdiv = []        # Retrieval flags
      self.co_nfitlayers = []  # No.layers actually retrieved
      self.co_cp_air = []      # Air partial columns in each layer [molec/cm2]
      self.co_cp_a = []        # A Priori partial columns [molec/cm2]
      self.co_cp   = []        # Retrieved partial columns [molec/cm2]
      # Repeat for HNO3
      self.hno3_qflag = []
      self.hno3_bdiv = []
      self.hno3_nfitlayers = []
      self.hno3_cp_air = []
      self.hno3_cp_a = []
      self.hno3_cp   = []
      # Repeat for O3
      self.o3_qflag = []
      self.o3_bdiv = []
      self.o3_nfitlayers = []
      self.o3_cp_air = []
      self.o3_cp_a = []
      self.o3_cp   = []
      # BRESCIA SO2 data, different structure to FORLI
      self.so2_qflag = []            # General retrieval quality flag
      self.so2_altitude = []         # Retrieved plume altitude [m]
      self.so2_col_at_altitudes = [] # SO2 column for plume at diff.est. altitudes [DU]
      self.so2_col = []              # SO2 col.at rtv.plume alt from OEM [DU]
      self.so2_bt_difference = []    # Indicative bright.temp.difference [K]



  def _latlim ( self, mdr, latlim ):
    """ Limit pixels to within latitude range

    Parameters
      mdr    obj : Measurement data record
      latlim [2] : Min,Max latitudes [deg N] for inclusion

    Returns
      [boo]  : List of 120 flags, True=fov within lat.limits

    """
    imin = round ( latlim[0]*1e4 ) # scale to match MDR units 
    imax = round ( latlim[1]*1e4 ) 
    ilat = np.array ( mdr['earth_location'][:,0] )
    return ( ilat >= imin ) & ( ilat <= imax )

  def _loclim ( self, mdr, loclst, locmsc ):
    """ Limit pixels to match L1C times

    Parameters
      mdr    obj  : Measurement data record
      loclst [{}] : List of L1C locations

    Returns
      [boo]  : List of 120 flags, True=fov within loc list

    """
    imatch = np.zeros(self.NFOV,dtype=bool)
    msclin = mdr['grh']['record_start_time'][1] 
    if msclin not in locmsc: return imatch
    iloc = 0
    for loc in loclst:
      if loc['msclin'] == msclin:
        ifov = loc['stp']*4 + loc['pix'] - 5
        if mdr['fg_atmospheric_temperature'][ifov,0] == self.BADU2:
          print('W-iasi_l2_class._loclim: Invalid FG Tem - exclude profile')
        else:
          imatch[ifov] = True
          self.useloc[iloc] = True
      iloc += 1
    return imatch

  def _lonlim ( self, mdr, lonlim ):
    """ Limit pixels to within longitude range

    Parameters
      mdr    obj : Measurement data record
      lonlim [2] : Min,Max longitudes [deg E] for inclusion

    Returns
      [boo]  : List of 120 flags, True=fov within lon.limits

    """
    imin = round ( lonlim[0]*1e6 ) # scale to match MDR units 
    imax = round ( lonlim[1]*1e6 ) 
    ilon = np.array ( mdr['earth_location'][:,1] )
    if imin < imax:      # conventional lon. coverage
      return ( ilon >= imin ) & ( ilon <= imax )
    else:                # reqd range spans dateline
      return ( ilon >= imin ) | ( ilon <= imax )


  def _read_err ( self, f, skip=True ):
    """ Read error covariance part of MDR record

    Parameters
      f     obj : file object
      skip  boo : True=skip over this part

    Returns
      {}  : dictionary of contents of err part of record (if not skip)

    """
    # Error data - variable size, depends on nerr, so read separately
    nerr = int ( np.fromfile(f,dtype='b',count=1)[0] )
    if skip: 
      sizerr = 120 + nerr * ( self.NERRT + self.NERRW + self.NERRO ) * 4
      f.seek(sizerr,1)  
      return None
    err = {}
    err['nerr'] = nerr                                             # 1b
    # The following error data occupies 120 + nerr*(NERRT+NERRW+NERRO)*4 bytes
    err['error_data_index']   = np.fromfile ( f, dtype='b',count=self.NFOV )  # 120b
    if nerr > 0:
      err['temperature_error']  = np.fromfile ( f, dtype='>u4',count=nerr*self.NERRT ) # 
      err['water_vapour_error'] = np.fromfile ( f, dtype='>u4',count=nerr*self.NERRW ) # 
      err['ozone_error']        = np.fromfile ( f, dtype='>u4',count=nerr*self.NERRO ) # 
    return err

  def _read_giadr ( self, f, offset ): 
    """ Read GIADR 

    Parameters
      f      obj : file object
      offset int : byte offset for record

    Returns
      giadr  obj : structure of GIADR, as defined below
    """  

    # L2 Global Internal Auxiliary Data Record (p59 EUM/OPS-ES/MAN/04/0033)
    GIADR = np.dtype([           ('grh', self.GRH),             #  20b    20
            ('num_pressure_levels_temp',   'b'),                #   1b    21
                ('pressure_levels_temp', '>u4', self.NLT),      # 404b   425
        ('num_pressure_levels_humidity',   'b'),                #   1b   426
            ('pressure_levels_humidity', '>u4', self.NLQ),      # 404b   830
           ('num_pressure_levels_ozone',   'b'),                #   1b   831
               ('pressure_levels_ozone', '>u4', self.NLO),      # 404b  1235  
  ('num_surface_emissivity_wavelengths',   'b'),                #   1b  1236
      ('surface_emissivity_wavelengths', '>u4', self.NEW),      #  48b  1284
                 ('num_temperature_pcs',   'b'),                #   1b  1285
                ('num_water_vapour_pcs',   'b'),                #   1b  1286
                       ('num_ozone_pcs',   'b'),                #   1b  1287
                 ('forli_num_layers_co',   'b'),                #   1b  1288
              ('forli_layer_heights_co', '>u2', self.NL_CO),    #  38b  1326
               ('forli_num_layers_hno3',   'b'),                #   1b  1327
            ('forli_layer_heights_hno3', '>u2', self.NL_HNO3),  #  38b  1365
                 ('forli_num_layers_o3',   'b'),                #   1b  1366
              ('forli_layer_heights_o3', '>u2', self.NL_O3),    #  80b  1446
           ('brescia_num_altitudes_so2',   'b'),                #   1b  1447
               ('brescia_altitudes_so2', '>u2', self.NL_SO2) ]) #  10b  1457
    GIADR_SIZE = 1457                                           #  was 824 in v5

    f.seek(offset) 
    giadr = np.fromfile ( f, dtype=GIADR, count=1 )[0]
    return giadr

  def _read_ipr ( self, f ):
    """ Read Internal Pointer Records

    Parameters
      f obj : file object

    Returns
      int : byte offset GIADR
      int : byte offset for first MDR

    Description
      IPRs are index records which point to the location within the file of
        'target' records of specific types.
      The IPRs follow immediately after the MPH, and the number of IPRs is
        given in the MPH as TOTAL_IPR.
      This returns the byte offsets for the first occurrence of the
        following target records
        giadr : global internal auxiliary data record (1 expected)
        mdr   : measurement data, ie records corresponding to each scan line
      If there is a loss of data, a specific type of 'lost data' record is
        inserted with an IPR pointing to the start of this section, followed 
        by a second MDR IPR pointing to the next record of nominal 
        measurements, etc.
    """

    # L2 Internal Pointer Record (p54 of EPS.GGS.SPE.96167)
    IPR = np.dtype([ 
            ('grh', self.GRH ),                   # 20b  20b
            ('target_record_class',     'b'),     #  1b  21b  (5=GIADR, 8=MDR)
            ('target_instrument_group', 'b'),     #  1b  22b  (15=IASI_L2, 13=Missing data)
            ('target_record_subclass',  'b'),     #  1b  23b
            ('target_record_offset',  '>u4')  ])  #  4b  27b total

    # Identify internal pointer records associated with GIADR and first MDR
    # See Table 14, p28 of  EPS.GGS.SPE.96167
    nipr = int ( self.mph['TOTAL_IPR'] )          # nipr=2 expected, unless missing data
    ipr_mdr = None
    for i in range(nipr):
      ipr = np.fromfile(f,dtype=IPR,count=1)[0]
      #print(ipr['target_record_class'])
      if ipr['target_record_class'] == 5: 
        ipr_giadr = ipr['target_record_offset']
      elif ipr['target_record_class'] == 8:
        if ipr_mdr:
          print('W-Iasi_L2_Class: multiple MDR IPR records, only taking first')
        else:
          ipr_mdr   = ipr['target_record_offset']
      else: 
        self.fail = True
        self.errtxt = 'F-Iasi_L2_Class: Unexpected Internal Pointer Record,' +\
                      ' Class=' + ipr['target_record_class']
        return None, None
    return ipr_giadr, ipr_mdr

  def _read_mph ( self, f ):
    """ Read Main Product Header 

    Parameters
      f obj : file object

    Returns
      {str} : dictionary of MPH parameter=value pairs

    Description
      The MPH is an arbitrary series of parameter=value text strings.
      These are read and returned as a dictionary mph (parameters as keys)
      This method assumes the file pointer is initially located at the start of 
      the MPH record (which should be the start of the file), and leaves the 
      pointer at the end of the MPH record.
    """

    # First read GRH to establish total size of MPH (including GRH part)      
    grh = np.fromfile(f,dtype=self.GRH,count=1)[0]
    mphsiz = grh['record_size'] - self.GRH_SIZE  # subtract GRH component
    rec = np.fromfile(f,dtype='b',count=mphsiz)  # Read remaining MPH as bytes
    mphrec = rec.tostring().decode('utf-8')      # convert bytes to characters
    ipt = 0
    mph = {}                                     # create dictionary
    while ipt < mphsiz:
      ieq = mphrec.find('=',ipt)                 # locate next '=' sign
      icr = mphrec.find('\n',ipt)                # locate next end-of-line mark
      param = mphrec[ipt:ieq-1].strip()          # use left of '=' for key 
      value = mphrec[ieq+1:icr].strip()          # use right of '=' for value
      mph[param] = value
      ipt = icr + 1                              # move to next line 
    return mph

  def _read_forli ( self, f, nl, neve, neva ):
    """ Read FORLI data for one species
    """

    forli = {}
    forli['nl'] = nl  # required for unpacking 2D arrays later
    # FORLI - fixed length part
    forli['qflag']      = np.fromfile ( f, dtype='b',  count=self.NFOV )   # 120b
    forli['bdiv']       = np.fromfile ( f, dtype='>u4',count=self.NFOV )   # 480b
    forli['npca']       = np.fromfile ( f, dtype='b',  count=self.NFOV )   # 120b
    forli['nfitlayers'] = np.fromfile ( f, dtype='b',  count=self.NFOV )   # 120b
    forli['nbr']        = np.fromfile ( f, dtype='b',  count=1)[0]         #   1b
    nbr = int(forli['nbr'])     # no profiles actually retrieved in scan line 
    # FORLI - variable length part
    if nbr > 0: 
      forli['cp_air']   = np.fromfile ( f, dtype='>u2',count=nl*nbr ) # nl*nbr*2
      forli['cp_a']     = np.fromfile ( f, dtype='>u2',count=nl*nbr ) # nl*nbr*2
      forli['x']        = np.fromfile ( f, dtype=self.VUI2, count=nl*nbr ) # nl*nbr*3
      forli['h_eigenvalues']  = np.fromfile ( f, dtype=self.VI4, count=neva*nbr ) 
      forli['h_eigenvectors'] = np.fromfile ( f, dtype=self.VI4, count=neve*nbr ) 
    # Also useful to construct list of indices for profile data from each FOV
    forli['idxdat'] = [-1]*self.NFOV
    idat = 0
    for ifov in range(self.NFOV):
      if forli['qflag'][ifov] != 0:
        forli['idxdat'][ifov] = idat
        idat += 1
    return forli


  def _read_brescia ( self, f, nl ):
    """ Read BRESCIA SO2 data 
    """
    brescia = {}
    brescia['nl'] = nl
    brescia['qflag']            = np.fromfile ( f, dtype='b', count=self.NFOV )  # 120b
    brescia['col_at_altitudes'] = np.fromfile ( f, dtype='>u2', count=self.NFOV*nl )     
    brescia['altitude']         = np.fromfile ( f, dtype='>u2', count=self.NFOV ) 
    brescia['col']              = np.fromfile ( f, dtype='>u2', count=self.NFOV ) 
    brescia['bt_difference']    = np.fromfile ( f, dtype='>u2', count=self.NFOV ) 
    return brescia

  def _set_loc ( self, mdr, imdr, idxfov ):
    """ Set output location data 

    Parameters
      mdr     obj  : Measurement data record (ie IASI scan line)
      imdr    int  : Index of MDR within file (starting at 0, up to c.760)
      idxfov [int] : List of indices of FOV locations to use (0:199)

    Returns
      None - sets class attributes 

    """     
    self.stp += ( 1 + idxfov//4 ).tolist()   #  1, 1, 1, 1, 2, 2, ... scan step# (1:30)
    self.lin += [ 1 + imdr ] * len(idxfov)   #  1, 2, ...             MDR# within file/orbit
    self.pix += ( 1 + idxfov%4 ).tolist()    #  1, 2, 3, 4, 1, 2 ...  FOV# with line (1:120)
    self.fov += ( 1 + idxfov ).tolist()      #  1, 2, 3, 4, 5, 6, ... Pixel# within step (1:4)
    self.lat += ( mdr['earth_location'][idxfov,0] * 1.0E-4 ).tolist()    # latitutde in deg
    self.lon += ( mdr['earth_location'][idxfov,1] * 1.0E-4 ).tolist()    # longitude in deg
    self.sza += ( mdr['angular_relation'][idxfov,0] * 1.0E-2 ).tolist()  # Solar.zen. in deg
    self.zen += ( mdr['angular_relation'][idxfov,1] * 1.0E-2 ).tolist()  # Satel.zen. in deg
    self.saa += ( mdr['angular_relation'][idxfov,2] * 1.0E-2 ).tolist()  # Soalr.azi. in deg
    self.azi += ( mdr['angular_relation'][idxfov,3] * 1.0E-2 ).tolist()  # Satel.azi. in deg
    self.daylin += [ mdr['grh']['record_start_time'][0] ] * len(idxfov)  # Jul.Day at start of MDR
    self.msclin += [ mdr['grh']['record_start_time'][1] ] * len(idxfov)  # Millsec at start of MDR

  def _set_flags ( self, mdr, idxfov ):
    """ Set flags
    """     
    self.flg_amsubad   += ( mdr['flg_amsubad'][idxfov] ).tolist()
    self.flg_avhrrbad  += ( mdr['flg_avhrrbad'][idxfov] ).tolist()
    self.flg_cldfrm    += ( mdr['flg_cldfrm'][idxfov] ).tolist()
    self.flg_cldnes    += ( mdr['flg_cldnes'][idxfov] ).tolist()
    self.flg_cldtst    += ( mdr['flg_cldnes'][idxfov] ).tolist()
    self.flg_daynit    += ( mdr['flg_daynit'][idxfov] ).tolist()
    self.flg_dustcld   += ( mdr['flg_dustcld'][idxfov] ).tolist()
    self.flg_fgcheck   += ( mdr['flg_fgcheck'][idxfov] ).tolist()
    self.flg_iasibad   += ( mdr['flg_iasibad'][idxfov] ).tolist()
    self.flg_initia    += ( mdr['flg_initia'][idxfov] ).tolist()
    self.flg_itconv    += ( mdr['flg_itconv'][idxfov] ).tolist()
    self.flg_lansea    += ( mdr['flg_lansea'][idxfov] ).tolist()
    self.flg_mhsbad    += ( mdr['flg_mhsbad'][idxfov] ).tolist()
    self.flg_numit     += ( mdr['flg_numit'][idxfov] ).tolist()
    self.flg_nwpbad    += ( mdr['flg_nwpbad'][idxfov] ).tolist()
    self.flg_physcheck += ( mdr['flg_physcheck'][idxfov] ).tolist()
    self.flg_retcheck  += ( mdr['flg_retcheck'][idxfov] ).tolist()
    self.flg_satman    += ( mdr['flg_satman'][idxfov] ).tolist()
    self.flg_sunglnt   += ( mdr['flg_sunglnt'][idxfov] ).tolist()
    self.flg_thicir    += ( mdr['flg_thicir'][idxfov] ).tolist()

  def _set_profs ( self, mdr, idxfov ):
    """ Set profile data
    """
    for ifov in idxfov:
      # A priori profiles
      prf = mdr['fg_atmospheric_temperature'][ifov,:]
      self.tem_a.append( np.where ( prf == self.BADU2, np.nan, prf*self.TSCALE ) )
      prf = mdr['fg_atmospheric_water_vapour'][ifov,:]
      self.h2o_a.append( np.where ( prf == self.BADU4, np.nan, prf*self.QSCALE ) )
      prf = mdr['fg_atmospheric_ozone'][ifov,:] 
      self.o3_a.append  ( np.where ( prf == self.BADU2, np.nan, prf*self.OSCALE ) )

      # First guess/a priori quality indicators
      self.tem_aq.append( mdr['fg_qi_atmospheric_temperature'][ifov] ) 
      self.h2o_aq.append( mdr['fg_qi_atmospheric_water_vapour'][ifov] ) 
      self.o3_aq.append( mdr['fg_qi_atmospheric_ozone'][ifov] ) 

      # Retrieved profiles
      prf = mdr['atmospheric_temperature'][ifov,:]
      self.tem.append( np.where ( prf == self.BADU2, np.nan, prf*self.TSCALE ) )
      prf = mdr['atmospheric_water_vapour'][ifov,:]
      self.h2o.append( np.where ( prf == self.BADU4, np.nan, prf*self.QSCALE ) )
      prf = mdr['atmospheric_ozone'][ifov,:] 
      self.o3.append( np.where ( prf == self.BADU2, np.nan, prf*self.OSCALE ) )

  def _set_cols ( self, mdr, idxfov ):
    """ Set column data
    """
    # Integrated column amounts are all kg/m^2
    col = mdr['integrated_water_vapour'][idxfov]
    self.h2o_col += ( np.where ( col == self.BADU2, np.nan, col*1.0E-2 ) ).tolist()
    col = mdr['integrated_ozone'][idxfov]
    self.o3_col  += ( np.where ( col == self.BADU2, np.nan, col*1.0E-6 ) ).tolist()
    col = mdr['integrated_n2o'][idxfov]
    self.n2o_col += ( np.where ( col == self.BADU2, np.nan, col*1.0E-6 ) ).tolist()
    col = mdr['integrated_co'][idxfov]
    self.co_col  += ( np.where ( col == self.BADU2, np.nan, col*1.0E-7 ) ).tolist()
    col = mdr['integrated_ch4'][idxfov]
    self.ch4_col += ( np.where ( col == self.BADU2, np.nan, col*1.0E-6 ) ).tolist()
    col = mdr['integrated_co2'][idxfov]
    self.co2_col += ( np.where ( col == self.BADU2, np.nan, col*1.0E-3 ) ).tolist()

  def _set_sfc ( self, mdr, mdr_sfz, idxfov ):
    """ Set surface data
    """
    # Surface Temperature
    sft = mdr['fg_surface_temperature'][idxfov]
    self.sft_a  += ( np.where ( sft == self.BADU2, np.nan, sft*self.TSCALE ) ).tolist()
    sft = mdr['fg_qi_surface_temperature'][idxfov] 
    self.sft_aq += ( np.where ( sft == self.BADU2, np.nan, sft*self.TSCALE ) ).tolist() 
    sft = mdr['surface_temperature'][idxfov]
    self.sft    += ( np.where ( sft == self.BADU2, np.nan, sft*self.TSCALE ) ).tolist()

    # Surface pressure
    sfp = mdr['surface_pressure'][idxfov]
    self.sfp += ( np.where ( sfp == self.BADU4, np.nan, sfp*0.01 ) ).tolist()      # hPa

    # Surface altitude
    sfz = mdr_sfz[idxfov]
    self.sfz += ( np.where ( sfz == self.BADI2, np.nan, sfz ) ).tolist()   # m

    # Surface emissivity 
    for ifov in idxfov:
      sfe = mdr['surface_emissivity'][ifov,:] 
      self.sfe.append( np.where ( sfe == self.BADU2, np.nan, sfe*1.0E-4 ) ) 

  def _set_cloud ( self, mdr, idxfov ):
    """ Set cloud data 
    """
    # Cloud parameters
    for ifov in idxfov:
      ncld = mdr['number_cloud_formations'][ifov]
      self.ncld.append( ncld ) 
      self.cld.append ( mdr['fractional_cloud_cover'][ifov,0:ncld]*0.01 ) # %
      self.ctt.append ( mdr['cloud_top_temperature'][ifov,0:ncld]*0.01 )  # K
      self.ctp.append ( mdr['cloud_top_pressure'][ifov,0:ncld]*0.01 )     # hPa
      self.cph.append ( mdr['cloud_phase'][ifov,0:ncld] )     

  def _set_forli ( self, forli, idxfov, scalef, molec ):
    """ Convert FORLI to profile data
    """
    qflag = forli['qflag'][idxfov].tolist()
    bdiv  = forli['bdiv'][idxfov].tolist()
    nfitlayers = forli['nfitlayers'][idxfov].tolist()
    if forli['nbr'] > 0:
      forli_cp_air = np.reshape ( forli['cp_air'], ( -1, forli['nbr'] ) )
      forli_cp_a   = np.reshape ( forli['cp_a'],   ( -1, forli['nbr'] ) )
      xf = forli['x'][:]['i2man'] / 10.0**forli['x'][:]['ibexp']
      xf = np.reshape ( xf, ( -1, forli['nbr'] ) )
    cp_air = []
    cp_a   = []
    cp     = []
    for ifov in idxfov:
      idat = forli['idxdat'][ifov]
      if idat < 0: 
        cp_air.append(np.nan)
        cp_a.append(np.nan)
        cp.append(np.nan)
      else:
        cp_air.append( forli_cp_air[:,idat]*1.0E20 ) 
        cp_a.append(     forli_cp_a[:,idat]*scalef ) 
        cp.append(       forli_cp_a[:,idat]*scalef * xf[:,idat] ) 
    if molec == 'co': 
      self.co_qflag      += qflag
      self.co_bdiv       += bdiv
      self.co_cp_air     += cp_air
      self.co_cp_a       += cp_a
      self.co_cp         += cp
      self.co_nfitlayers += nfitlayers
    elif molec == 'o3': 
      self.o3_qflag      += qflag
      self.o3_bdiv       += bdiv
      self.o3_cp_air     += cp_air
      self.o3_cp_a       += cp_a
      self.o3_cp         += cp
      self.o3_nfitlayers += nfitlayers
    elif molec == 'hno3': 
      self.hno3_qflag      += qflag
      self.hno3_bdiv       += bdiv
      self.hno3_cp_air     += cp_air
      self.hno3_cp_a       += cp_a
      self.hno3_cp         += cp
      self.hno3_nfitlayers += nfitlayers
    else:
      print('F-_set_forli: unexpected value of "molec": ' + molec )
      exit()      

  def _set_brescia ( self, brescia, idxfov ):
    """ Convert BRESCIA to profile data
    """
    qflag = brescia['qflag'][idxfov].tolist()

    balt = brescia['altitude'][idxfov]
    altitude = ( np.where ( balt == self.BADU2, np.nan, balt*1.0 ) ).tolist()  # already in m

    bcol = brescia['col'][idxfov]
    col = ( np.where ( bcol == self.BADU2, np.nan, bcol*0.1 ) ).tolist()  # 0.1 for DU

    btd = brescia['bt_difference'][idxfov]
    bt_difference = ( np.where ( btd == self.BADI2, np.nan, btd*0.01 ) ).tolist()  # 0.01 for K

    bcaa = np.reshape ( brescia['col_at_altitudes'], (-1,self.NFOV) )
    col_at_altitudes = []
    for ifov in idxfov:
      col_at_altitudes.append( np.where ( bcaa[:,ifov] == self.BADU2, np.nan, 
                                          bcaa[:,ifov] * 0.1 ) )   # 0.1 to conver to DU
    self.so2_qflag    += qflag
    self.so2_col      += col
    self.so2_altitude += altitude
    self.so2_col_at_altitudes += col_at_altitudes


#l2dir = '/network/group/aopp/sat_data/SAT015_METOP_IASI_L1C_GLOBAL/IASI-A_L2/2018/06/'
#l2fil = l2dir + \
'IASI_SND_02_M02_20180618041455Z_20180618055655Z_N_O_20180618055859Z.nat'

#'IASI_SND_02_M02_20180618155959Z_20180618174159Z_N_O_20180618174149Z.nat'


#l2 = Iasi_L2 ( l2fil, forli=True )
#if l2.fail: print(l2.errtxt)

# print total number of locations in output
#print(l2.nloc)
