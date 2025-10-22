"""If running the SRFM with a driver table, run this code.
"""

import driver_table as dt
from srfm import *

iasi_main.run_srfm(dt.x,dt.b)

