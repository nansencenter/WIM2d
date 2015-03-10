import fns_grid_setup as gs
import os
import sys
import matplotlib.pyplot as plt

GRID_OPT = 1 # Grid configuration
TEST     = 0 # if 1: don't save grid files to ../run/inputs 
             # if 0: save grid files to test/out_py
             # *Both options make a plot of the LAND MASK
LAND_OPT = 2
gs.grid_setup(GRID_OPT=GRID_OPT,TEST=TEST,LAND_OPT=LAND_OPT)
