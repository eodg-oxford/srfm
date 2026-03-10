"""If running the SRFM with a driver table, run this code."""
from srfm import *

srfm_inps = inputs.Inputs()
srfm_inps.read_srfm_drv("driver_table.py")

srfm_res = main.run_srfm(srfm_inps)
