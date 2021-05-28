# INPUT
noise_source_output_file = "source1.h5"
grid_file = "sourcegrid.npy"
surf_file = "surfel_temp.npy"
src_model_file = "oceanmodels/pressure_PSD_2007-12-11-00.npz"
greens_function_file = "pressure1d.displacement..MXZ.h5"
# Number of basis functions for representing spectra:
n_basis = 100
# represent spectra between fmin and fmax (noisi source needs all freq. between 
# 0 and Nyquist)
fmin = 0.0
# End input

import h5py
import numpy as np
from scipy.signal import tukey
from scipy.interpolate import interp1d
from sincbasis import SincBasis


def run_s3(noise_source_output_file, grid_file, surf_file,
           src_model_file, greens_function_file, n_basis, fmin):
    microseism_model = np.load(src_model_file)
    # frequency axis
    with h5py.File(greens_function_file) as wf:
        df = wf["stats"].attrs["Fs"]
        n = wf["stats"].attrs["npad"]

    freqs_new = np.fft.rfftfreq(n, d=1. / df)
    source_grid = np.load(grid_file)
    latmin = source_grid[1].min()
    latmax = source_grid[1].max()
    lonmin = source_grid[0].min()
    lonmax = source_grid[0].max()

    # pressure model
    print('min frequency', microseism_model['freqs'].min())
    print('max frequency', microseism_model['freqs'].max())
    print('min latitude', microseism_model['lats'].min())
    print('max latitude', microseism_model['lats'].max())
    print('min longitude', microseism_model['lons'].min())
    print('max longitude', microseism_model['lons'].max())
    print("Pressure source on ", microseism_model["date"], "\nUnit ",
          microseism_model["pressure_PSD_unit"])
    freqs_old = microseism_model["freqs"]
    sampling_frequency = np.diff(freqs_old).mean()
    lats = microseism_model["lats"]
    ixlats = np.intersect1d(np.where(lats >= latmin), np.where(lats <= latmax))
    lats = lats[ixlats]
    lons = microseism_model["lons"]
    ixlons = np.intersect1d(np.where(lons >= lonmin), np.where(lons <= lonmax))
    lons = lons[ixlons]
    p = microseism_model["pressure_PSD"]
    p = p[:, ixlats[0]: ixlats[-1]+1, ixlons[0]: ixlons[-1]+1]
    p = np.nan_to_num(p)

    # these data are gridded so we need to un-wrap them for plotting
    lats_c = np.zeros(len(lats) * len(lons))
    lons_c = np.zeros(len(lats) * len(lons))
    p_c = np.zeros(len(lats) * len(lons) * len(freqs_old))

    for k in range(len(freqs_old)):
        for i, lat in enumerate(lats):
            for j, lon in enumerate(lons):
                lats_c[j + len(lons) * i] = lat
                lons_c[j + len(lons) * i] = lon
                p_c[j + len(lons) * i + len(lats) * len(lons) * k] = p[k, i, j]

    # now, we first have to create a spectral basis for the new model
    sb = SincBasis(K=n_basis, N=len(freqs_new), freq=freqs_new,
                  fmin=fmin, fmax=freqs_new.max())

    spect_basis =np.zeros((n_basis, len(freqs_new)), dtype=np.float)

    # add vectors
    for i in range(n_basis):
        spect_basis[i, :] = sb.basis_func(k=i, n=len(freqs_new))

    # coefficient matrix
    n_loc = source_grid.shape[1]
    mod = np.zeros((n_loc, n_basis), dtype=np.float)

    # for each location, we have to fit the spectrum of the microseism model
    # this is a bit slow
    for i in range(n_loc):
        
        # nearest point, approximately (this is not exact, do we need more exact?)
        ix_lon = np.argmin((lons - source_grid[0, i]) ** 2)
        ix_lat = np.argmin((lats - source_grid[1, i]) ** 2)
        coeffs = np.zeros(n_basis)
        spec = p[:, ix_lat, ix_lon]
        
        # PSD-->amplitude
        spec *= sampling_frequency  # divide by time interval
        spec /= (50000. ** 2)  # divide by spatial sampling rate in 2 D APPROXIMATE
        spec = np.sqrt(spec)
        spec[spec <= 1.e-10] = 0.0

        if i % 1000 == 0:
            print(i, end=",")
        
        f = interp1d(freqs_old, spec, bounds_error=False, fill_value=0, kind="cubic")
        spec_n = f(freqs_new)

        for j in range(n_basis):
            coeffs[j] = np.dot(spect_basis[j, :], spec_n)
            
        # and finally enter the coefficients in the new model
        mod[i, :] = coeffs

    # Save to an hdf5 file
    with h5py.File(noise_source_output_file, 'w') as fh:
        fh.create_dataset('coordinates', data=np.load(grid_file))
        fh.create_dataset('frequencies', data=freqs_new)
        fh.create_dataset('surface_areas',
                          data=np.load(surf_file))
        fh.create_dataset("spectral_basis", data=spect_basis)
        fh.create_dataset("model", data=mod)

if __name__ == "__main__":

    run_s3(noise_source_output_file, grid_file, surf_file,
           src_model_file, greens_function_file, n_basis, fmin)