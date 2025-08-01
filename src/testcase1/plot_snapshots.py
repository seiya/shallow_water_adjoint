import glob
import os
import numpy as np
import matplotlib.pyplot as plt


def read_grid_params(fname="grid_params.txt"):
    """Return grid dimensions (nlon, nlat) read from a text file."""
    with open(fname) as f:
        parts = f.read().split()
    if len(parts) < 2:
        raise ValueError(f"{fname} must contain at least two integers")
    return int(parts[0]), int(parts[1])


nlon, nlat = read_grid_params()
dt = 600.0

# Construct longitude/latitude arrays in degrees
pi = np.pi
dlon = 2.0 * pi / nlon
dlat = pi / nlat
lon = np.arange(nlon) * dlon
lat = -pi / 2.0 + (np.arange(nlat) + 0.5) * dlat
lon_deg = np.degrees(lon)
lat_deg = np.degrees(lat)
lon2, lat2 = np.meshgrid(lon_deg, lat_deg, indexing='ij')

npts = nlon * nlat
for fname in sorted(glob.glob('snapshot_*.bin')):
    data = np.fromfile(fname, dtype=np.float32)
    if data.size != 3 * npts:
        raise ValueError(f'Unexpected data size in {fname}')
    h = data[0:npts].reshape((nlon, nlat), order='F')
    u = data[npts:2*npts].reshape((nlon, nlat), order='F')
    v = data[2*npts:].reshape((nlon, nlat), order='F')

    fig, ax = plt.subplots(figsize=(8, 4))
    cs = ax.pcolormesh(lon_deg, lat_deg, h.T, shading='auto')
    ax.quiver(lon2[::2, ::2], lat2[::2, ::2], u.T[::2, ::2], v.T[::2, ::2])

    nstep = int(os.path.splitext(fname)[0].split('_')[1])
    thours = nstep * dt / 3600.0
    ax.set_title(f't = {thours:.1f} h')
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    plt.colorbar(cs, ax=ax, label='Height (m)')
    pngname = os.path.splitext(fname)[0] + '.png'
    plt.savefig(pngname)
    plt.close(fig)
