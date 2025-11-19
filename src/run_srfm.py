"""If running the SRFM with a driver table, run this code.
"""

import driver_table as dt
from srfm import *

srfm_inps = iasi_main.Srfm_iasi_inputs()
srfm_inps.read_srfm_drv("driver_table.py")

srfm_res = iasi_main.run_srfm(srfm_inps)
