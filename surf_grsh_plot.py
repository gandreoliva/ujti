import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

r"""
Plots the gravitational+Doppler shift as seen by the observer
"""
# make_ujti.sh -->  ztoa_make_run.sh  --> lchotspot_get.py --> lchotspot_plot.py
#                                     --> thspectrum_get_plot.py
#                                     --> ((surf_*_plot.py))

def main(filename):

    data_type = np.dtype([
        ('posix',float),
        ('posiy',float),
        ('posiz',float),
        ('z',float),
        ('posfx',float),
        ('posfy',float),
        ('posfz',float),
        ('x0',float)
        ])

    data = np.fromfile(filename,dtype=data_type)

    plot_x_grid = data['posiy'].reshape(nx,ny)
    plot_y_grid = data['posiz'].reshape(nx,ny)
    z = data['z'].reshape(nx,ny)
    zmax = np.amax(z[~np.isnan(z)])

    plt.pcolormesh(plot_x_grid,plot_y_grid,z, vmin=-zmax, vmax=zmax, cmap=cm.bwr)
    plt.colorbar(label=r"redshift")


if __name__ == "__main__":

    plt.figure(figsize=(6,5))
    plt.gca().set_aspect('equal')

    # (ztoa) Data ID
    dataid = r"SHFT-161-frutos-i0"
    main(f"data/zt/{dataid}.dat")
    metadata = np.loadtxt('data/zt/meta/'+dataid+'.txt')

    nx = ny = int(metadata[0])

    plt.show()
