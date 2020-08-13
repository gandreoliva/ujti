import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle


r"""
Plots the shape of the neutron star by plotting a grid on the apparent surface
"""

# make_ujti.sh -->  ztoa_make_run.sh  --> lchotspot_get.py --> lchotspot_plot.py
#                                     --> thspectrum_get_plot.py
#                                     --> ((surf_*_plot.py))


def classify_surfpoints(cart,image):
    pattern = cycle((50,200))

    r = np.sqrt(cart[0,...]**2 + cart[1,...]**2 + cart[2,...]**2)
    th = np.arccos(cart[2,...]/r)%np.pi
    ph = np.arctan2(cart[1,...],cart[0,...])%(2*np.pi)

    # divisions
    div_ph = np.linspace(0,360,17)*np.pi/180
    div_th = np.linspace(0,180,9)*np.pi/180

    # regions = np.array(np.meshgrid(div_th,div_ph)).T.reshape(-1,2)

    for (di,), _ in np.ndenumerate(div_th[:-1]):
        for (dj,), _ in np.ndenumerate(div_ph[:-1]):
            with np.errstate(invalid='ignore'): # ignore nan
                region = (th >= div_th[di]) & (th < div_th[di+1]) &\
                (ph >= div_ph[dj]) & (ph < div_ph[dj+1])
                image[region] = next(pattern)
        next(pattern)




def main(filename,frame=0,gen_movie=False):
    r"""
    filename : name where the data is located
    frame: name of the frame, if gen_movie is True
    gen_movie: if True, the figure is saved as png for later movie generation
    """

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

    pxf = data['posfx'].reshape(nx,ny)
    pyf = data['posfy'].reshape(nx,ny)
    pzf = data['posfz'].reshape(nx,ny)

    ## Redshift. If it's not in the surface of the ns, it's a nan
    z = data['z'].reshape(nx,ny)
    image = np.zeros((nx,ny),dtype=int)+255

    pxf[np.isnan(z)] = np.nan
    pyf[np.isnan(z)] = np.nan
    pzf[np.isnan(z)] = np.nan
    cart = np.array([pxf,pyf,pzf])
    classify_surfpoints(cart,image)

    plt.cla()
    plt.clf()

    borders = 240-(np.diff(image,axis=1,append=255)+np.diff(image,axis=0,append=255))/2
    borders[borders != 240] = 0
    borders[borders == 240] = 255

    plt.imshow(borders, cmap='gray', vmin=0, vmax=255,origin="lower")




if __name__ == "__main__":

    plt.figure(figsize=(6,5))
    plt.gca().set_aspect('equal')

    import sys

    # (toa) data id
    dataid = f"SHFT-161-frutos-i0"
    metadata = np.loadtxt('data/zt/meta/'+dataid+'.txt')

    nx = ny = int(metadata[0])

    main("data/zt/"+dataid+".dat")
    plt.title(dataid)
    plt.show()
