import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

r"""
Plots the light curve for a hotspot.
"""

# make_ujti.sh -->  ztoa_make_run.sh  --> lchotspot_get.py --> ((lchotspot_plot.py))
#                                     --> thspectrum_get_plot.py
#                                     --> surf_*_plot.py

default_norm_config={'t_offset':0, 't_norm':1, 'flux_norm':1}


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def plot_lc(dataid1,norm_config,color='black',label='',ls="-"):
    data_flux = np.loadtxt('data/lc/flux/{}.txt'.format(dataid1))

    every=1
    flux = np.convolve(data_flux[:,1], np.ones((every,))/every, mode='same')

    t = (data_flux[:,0]-norm_config['t_offset'])/norm_config['t_norm']
    flux_norm = flux/norm_config['flux_norm']

    axp.plot(t,flux_norm,ls,color=color,label=label)


def get_flux_normaliz(dataid1):
    data_flux = np.loadtxt('data/lc/flux/{}.txt'.format(dataid1))
    return np.amax(data_flux[:,1])

if __name__ == "__main__":

    plt.figure(figsize=(6,6))
    axp = plt.gca()

    axp.set_xlabel("Phase")

    hotspotid = "-circ_45_10-i0"
    objectid = "SHFT"
    resid = "-161"
    rot_freq = 716
    printing_press = False


    # normalization configuration
    norm_config = {
        't_offset': 0.00007,
        't_norm': 1/rot_freq,
        'flux_norm': get_flux_normaliz(objectid+resid+'-frutos'+hotspotid)
        }

    axp.set_xlim(0,1)


    plot_lc(objectid+resid+'-frutos'+hotspotid,norm_config,'dodgerblue',label='Fru16')


    axp.legend()
    plt.show()
