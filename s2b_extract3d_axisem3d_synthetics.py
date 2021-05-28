# ----------------------------------------------------
# INPUT
# ----------------------------------------------------
# Paths to axisem3d stuff:
# Path to simulation results:
indir = "/home/lermert/Dropbox/Kristiina_simulations/output_1D_pressureSource/stations/oceanfloor_stations/"
# Path to station list from axisem3d
stafile = "/home/lermert/Dropbox/Kristiina_simulations/STATIONS_0"

# Path to sourcegrid defined in step 1:
sourcegrid_file = "sourcegrid.npy"
# The result will be saved to the output file:
f_out = "pressure1d.displacement..MXZ.h5"
physical_quantity = "DIS"  # put P for pressure, DIS for displacement
# define a name for the "reference station",
# i.e. the location where the source was: will be saved just for reference
reference_station = "azim_0.station_1"
nt_keep = 2400  # throw away whatever comes after this timestep
# source coordinates in lon lat
src_x = -16.5
src_y = 37.5
earth_r = 6371000.0
# Force used for simulation in N:
input_N = 1.e10
# filter the synthetics with Butterworth with these parameters:
fmin=0.05
fmax=0.25
filtord=4

# ----------------------------------------------------
# END INPUT
# ----------------------------------------------------


# read in from Claudia's runs
import numpy as np
from netCDF4 import Dataset
import os
import h5py
from glob import glob
from obspy import Trace
from obspy.geodetics import gps2dist_azimuth, degrees2kilometers
from obspy.signal.invsim import cosine_taper
import matplotlib.pyplot as plt

# To do