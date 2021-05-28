from s1_makenoisegrid import run_s1
from s2a_extract1d_axisem3d_synthetics import run_s2
from s3_makenoisesource import run_s3
from s4_generate_timeseries import run_s4
import os

# static input
latmin = 35.5
latmax = 39.5
lonmin = -19.0
lonmax = -14.0
# source coordinates in lon lat
src_x = -16.5
src_y = 37.5
earth_r = 6371000.0
# temporary:
grid_filename = 'grids/sourcegrid.npy'
# temporary:
surfel_filename = "grids/surfel_temp.npy"

# Path to station list from axisem3d
stafile = "/home/lermert/Dropbox/Kristiina_simulations/STATIONS_0"
# define a name for the "reference station",
# i.e. the location where the source was: will be saved just for reference
reference_station = "azim_0.station_1"
nt_keep = 2400  # throw away whatever comes after this timestep
# Force used for simulation in N:
input_N = 1.e10
# filter the synthetics with Butterworth with these parameters:
fmin=0.05
fmax=0.25
filtord=4
greens_function_file = "greensfcts/src.greens..MXZ.h5"

# input related to noise source parametrization
noise_source_output_file = "noisesources/source1.h5"
# Number of basis functions for representing spectra:
n_basis = 100
# represent spectra between fmin and fmax (noisi source needs all freq. between 
# 0 and Nyquist)
fmin_interpolatespectrum = 0.0

# Variable input
grid_dx_in_ms = [5000, 10000, 15000]
# Path to Ocean source model:
src_model_files = ["oceanmodels/pressure_PSD_2007-12-04-00.npz",
                   "oceanmodels/pressure_PSD_2007-12-05-00.npz",
                   "oceanmodels/pressure_PSD_2007-12-06-00.npz",]
# Path to simulation results:
indirs = ["/home/lermert/Dropbox/Kristiina_simulations/output_1D_pressureSource/stations/oceanfloor_stations/",
         "/home/lermert/Dropbox/Kristiina_simulations/output_1D_forceSource/stations/oceansurface_stations/"]

f_out = "greensfcts/green..MXZ.h5"
# duration of the time series to be constructed
durations_seconds = [3600., 7200., 14400., 6. * 3600., 12. * 3600.]
# nr of sources for time domain dirac comb:
ns_sources = [10, 100, 1000, 10000]

for grid_dx_in_m in grid_dx_in_ms:
    run_s1(grid_dx_in_m, latmin, latmax, lonmin, lonmax,
                           grid_filename, surfel_filename)

    for indir in indirs:
        if "pressureSource" in indir:
            physical_quantity = "DIS"
            sourcename = "pressure1d"
            quantname = "displacement"
        else:
            physical_quantity = "P"
            sourcename = "force1d"
            quantname = "pressure"
        outfile_plot = "moveouts/moveout.{}.{}.{}km.png".format(sourcename, quantname, grid_dx_in_m//1000)

        run_s2(indir, stafile, grid_filename, greens_function_file,
                               physical_quantity, nt_keep, 
                               src_x, src_y, input_N, fmin, fmax, filtord,
                               earth_r, outfile_plot)

        for src_model_file in src_model_files:
            run_s3(noise_source_output_file, grid_filename, surfel_filename,
                           src_model_file, greens_function_file, n_basis, fmin_interpolatespectrum)

            for duration_seconds in durations_seconds:
                for n_sources in ns_sources:
                    src_model_name = "sm_" + os.path.splitext(os.path.basename(src_model_file))[0].split("_")[-1]
                    # source type, recording type, noise source name, grid dx in km, duration in seconds
                    outname = "noisetimeseries/{}.{}.{}.{}km.{}s".format(sourcename, quantname, src_model_name,
                                                         int(grid_dx_in_m/1000), int(duration_seconds))

                    run_s4(greens_function_file, noise_source_output_file, duration_seconds, n_sources, outname)