import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import sys

r"""
Plots the light curve for a hotspot.
"""

# make_ujti.sh -->  ztoa_make_run.sh  --> lchotspot_get.py --> lchotspot_plot.py
#                                     --> ((thspectrum_get_plot.py))
#                                     --> surf_*_plot.py

def black_body_intensity(freq,temp):
    h = 6.626070e-34 # J s
    c = 299792458 # m/s
    kB = 1.380649e-23 # J/K
    B = 2*h*freq**3/c**2 * 1/( np.exp(h*freq/(kB*temp)) - 1 )
    return B


def main(filename,label='',ls='-',color='black'):
    h = 6.626070e-34 # J s
    kB = 1.380649e-23 # J/K

    data_type = np.dtype([
        ('posix',float),
        ('posiy',float),
        ('posiz',float),
        ('z',float),
        ('posfx',float),
        ('posfy',float),
        ('posfz',float),
        ('x0', float)
        ])

    data = np.fromfile(filename,dtype=data_type)

    plot_x_grid = data['posiy'].reshape(nx,ny)
    plot_y_grid = data['posiz'].reshape(nx,ny)
    z = data['z'].reshape(nx,ny)

    ##  z = om_surf/om_inf - 1 =>  fratio = om_inf/om_surf

    fratio = (z + 1)**-1

    fratio = np.nan_to_num(fratio)

    ## di, dj : differential of solid angle

    di = dj = 0.001
    D = 1
    temp = 5000

    freqs_inf = np.logspace(13,15,100)
    Fluxes = np.zeros(freqs_inf.shape)


    for k in np.arange(freqs_inf.shape[0]):
        Flux = np.sum(  1/D**2*fratio[fratio!=0]**3 * black_body_intensity(freqs_inf[k]/fratio[fratio!=0], temp) * di * dj )

        Fluxes[k] = Flux

    norm_energy = h*freqs_inf/(kB*temp)
    axp.plot(norm_energy,Fluxes,ls,label=label,color=color)



if __name__ == "__main__":


    plt.figure(figsize=(6,8))
    axp = plt.gca()

    axp.set_ylabel(r'Flux [arbitrary]')
    axp.set_yscale('log')
    axp.set_xscale('log')
    axp.set_xlabel(r'$h\nu_\infty/(k_B T)$')

    datadir="data/zt/"
    objid = "SHFT"
    resid = "-161"
    postfix = "-i0.dat"

    nx = ny = 161

    main(datadir+objid+resid+"-frutos"+postfix,'Fru16',color='dodgerblue'); print('.')

    axp.legend()

    plt.show()
