# ----------------------------------------------------
# INPUT
# ----------------------------------------------------
# Paths to axisem3d stuff:
# Path to simulation results:
indir = "/home/lermert/Dropbox/Kristiina_simulations/output_1D_pressureSource/stations/oceanfloor_stations/"
# Path to station list from axisem3d
stafile = "/home/lermert/Dropbox/Kristiina_simulations/STATIONS_0"

# Path to sourcegrid defined in step 1:
sourcegrid_file = "grids/sourcegrid.npy"
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
outfile_plot = "moveout.pressuresource.png"
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


def run_s2(indir, stafile, sourcegrid_file, f_out,
           physical_quantity, nt_keep, 
           src_x, src_y, input_N, fmin, fmax, filtord=4,
           earth_r=6371000.0, outfile_plot=None):

    azs = np.arange(0., 360., 15)
    sourcegrid = np.load(sourcegrid_file)
    # read in stations -> rank and array index file
    fstas = open(os.path.join(indir, "rank_station.info"), "r")
    fstas = fstas.read().split("\n")
    sta_to_rank = {}
    for sta in fstas:
        if sta == "":
            continue
        inf = sta.split()
        sta_to_rank[inf[1]] = [inf[0], inf[2]]  # station string identifies rank, index in rank

    # get station locations from STATIONS file
    rec_az = []
    rec_dist = []
    rec_x = []
    rec_y = []
    fstas = open(stafile, "r")
    fstas = fstas.read().split("\n")
    for sta in fstas:
        if sta == "":
            continue
        inf = sta.split()
        rec_az.append(np.rad2deg(float(inf[3])))
        if float(inf[3]) == 0:
            rec_dist.append(float(inf[2]))   # distance in radians
    rec_az = np.array(list(set(rec_az)))
    rec_dist = np.array(list(set(rec_dist)))
    # print(rec_dist)
    rec_dist.sort()
    # set up the hdf output file
    f_out = h5py.File(f_out, "w")

    # open one file to get time vector
    testfile = Dataset(glob(os.path.join(indir, "*.nc*"))[0])
    t = testfile["data_time"][:]

    # DATASET NR 1: STATS
    stats = f_out.create_dataset('stats', data=(0,))
    stats.attrs['reference_station'] = reference_station
    stats.attrs['data_quantity'] = physical_quantity
    stats.attrs['ntraces'] = sourcegrid.shape[-1]
    fs = round(1./ np.diff(t).mean(), 5)
    stats.attrs['Fs'] = fs
    stats.attrs['nt'] = nt_keep
    stats.attrs['npad'] = 2 * len(t)  # TODO
    stats.attrs['fdomain'] = False

    # DATASET NR 2: Source grid
    f_out.create_dataset('sourcegrid', data=sourcegrid)

    # DATASET Nr 3: Seismograms itself
    data_out = f_out.create_dataset('data', (sourcegrid.shape[1], nt_keep),
                                    dtype=np.float)

    taper = cosine_taper(nt_keep, p=0.05)
    taper[0: len(taper) // 2] = 1.
    # for all noise source locations:
    ix_az = np.argmin(np.abs(rec_az))

    for i in range(sourcegrid.shape[-1]):
        x = sourcegrid[0, i]
        y = sourcegrid[1, i]

        # find the station(s) to interpolate from, 
        # start with nearest neighbour
        # I have no idea if this is correct, maybe it will be easier
        # to define a latlon grid of stations beforehand
        # ang_to_source = gps2dist_azimuth(0., 0., x, y)[1] / earth_r

        # use only distance to source for 1-D output
        if i == 0:
            print("Warning: using only azimuth=0 results")
        src_dist_km = gps2dist_azimuth(src_y, src_x, y, x)[0] / 1000.

        #print(ix_az)
        rec_dist_km = degrees2kilometers(np.rad2deg(rec_dist))
        #src_dist_km = degrees2kilometers(np.rad2deg(dist_to_source))
        ix_d1 = max(np.argsort((rec_dist_km - src_dist_km)**2)[0], 1)
        #print(ix_d1, rec_dist_km[ix_d1], src_dist_km)
        
        dd1 = np.abs(src_dist_km - rec_dist_km[ix_d1])
        stastr1 = "azim_0.station_{}".format(ix_d1)
        # get the time series
        rank_ix1 = sta_to_rank[stastr1]
        #print(rank_ix)
        dd = Dataset(os.path.join(indir, "axisem3d_synthetics.nc.rank{}".format(rank_ix1[0])))
        seis1 = dd["data_wave"][rank_ix1[1], -1, 0: nt_keep] / input_N
        seis1 *= taper
        trace1 = Trace(data=seis1)
        trace1.stats.sampling_rate = fs
        trace1.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=filtord)

        data_out[i, :] = trace1.data

        if i % 20 == 0 and outfile_plot is not None:
            traceplot = trace1.copy()
            plt.plot(np.arange(nt_keep) / traceplot.stats.sampling_rate, 
                    10*traceplot.data / traceplot.data.max() + rec_dist_km[ix_d1], 
                    color="0.7", alpha=0.5)

    if outfile_plot is not None:
        plt.xlabel("Seconds after source")
        plt.ylabel("Norm. waveforms at distance from source in km")
        plt.savefig(outfile_plot)

if __name__ == "__main__":
    run_s2(indir, stafile, sourcegrid_file, f_out,
           physical_quantity, nt_keep, 
           src_x, src_y, input_N, fmin, fmax, filtord,
           earth_r, outfile_plot)

