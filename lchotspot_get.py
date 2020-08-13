import numpy as np
import os,sys

r"""
Calculates the light curve for a circular hotspot

Version optimized for memory management (it writes t_cap to the disk).
Warning: the size of 't_exposure' has to be equal than the size of 'azimuths'!
"""

# make_ujti.sh -->  ztoa_make_run.sh  --> ((lchotspot_get.py)) --> lchotspot_plot.py
#                                     --> thspectrum_get_plot.py
#                                     --> surf_*_plot.py


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def filter_cap(cart,cap_colat,cap_az,cap_aperture):
    r"""
    Determines which Cartesian points `cart` are in a general spherical cap

    cart: positions to examine in the form
            np.array([[x1,x2,...,xn],[y1,y2,...,yn],[z1,z2,...,zn]])
    cap_colat: angle between the z axis and the center of the cap
    cap_az: angle between the x axis (obs's. loc.) and the proj. of the cap
        onto the xy plane
    cap_aperture: angle between the center of the cap and the boundary

    returns: array of booleans (dim. n); True if the point is in the cap
    note: this function ignores np.nan (returns False)
    """

    rotmaty = np.zeros((3,3))
    rotmaty[0,0] = np.cos(-cap_colat)
    rotmaty[0,2] = np.sin(-cap_colat)
    rotmaty[1,1] = 1
    rotmaty[2,0] = -np.sin(-cap_colat)
    rotmaty[2,2] = np.cos(-cap_colat)

    rotmatz = np.zeros((3,3))
    rotmatz[0,0] = np.cos(cap_az)
    rotmatz[0,1] = -np.sin(cap_az)
    rotmatz[1,0] = np.sin(cap_az)
    rotmatz[1,1] = np.cos(cap_az)
    rotmatz[2,2] = 1

    cart = rotmatz @ cart
    cart = rotmaty @ cart

    r = np.sqrt(cart[0]**2 + cart[1]**2 + cart[2]**2)
    th = np.arccos(cart[2]/r)
    ph = np.arctan2(cart[1],cart[0])


    with np.errstate(invalid='ignore'):
        filtered = (th > 0) & (th < cap_aperture)
        return filtered


def get_photo_mask(t_cap,t,dt_exposure):
    r"""
    takes the time from the observer `t`, and the arrival times from all the
    snapshots of the caps, `t_cap`, and determines which geodesics have
    arrived to the observer.

    returns: np.array shape nx,ny with True if the geodesic at that point
    has arrived to the observer, and False if it has not.
    """
    photo_mask = np.full((nx,ny), False)
    for each_t_cap in t_cap:
        with np.errstate(invalid='ignore'):
            photo_mask = photo_mask | ((each_t_cap > t) & (each_t_cap < (t + dt_exposure)))
    return photo_mask



def load_necessary_t_cap(i,padding,imin,imax,t_cap_slice):
    r"""
    loads only the necessary t_cap, and puts the data into the dictionary
    t_cap_slice
    """
    frames = np.arange(i-padding,i+padding//2)
    to_delete = set(t_cap_slice.keys()).difference(frames)
    for mark in to_delete:
        del t_cap_slice[mark]
    for j in frames:
        if j >= imin and j <= imax and j not in t_cap_slice:
            try:
                t_cap_slice[j] = np.load('data/lc/t_cap/{}'.format(dataid)+'/{:03d}.npz'.format(j))['arr_0']
            except:
                print("  (i) Tried to access non existent data")



def main(nx,ny,dataid, ns_freq, ns_Re, toa_filename, cap_colat, cap_aperture,padding):

    print("> Info: dataid={}, ns_freq={}, ns_Re={}, toa_filename={}, cap_colat={}, cap_aperture={}".format(dataid,ns_freq,ns_Re,toa_filename, cap_colat, cap_aperture))

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

    data = np.fromfile("data/zt/{}.dat".format(toa_filename),dtype=data_type)

    # < get delay in time of arrival >

    az_max = 395
    daz_cap = 0.5
    azimuths = np.arange(0,az_max,daz_cap)

    # i: index of both the azimuths and the t_exposure (for naming files)
    imin = 0
    imax = azimuths.shape[0]-1


    print("Getting times...")

    # dtg : delta t gravitational, i.e., the gravitational time delay

    pxf = data['posfx']
    pyf = data['posfy']
    pzf = data['posfz']

    ## Time of arrival. If it's not in the surface of the ns, it's a nan
    x0 = data['x0'].reshape(nx,ny)

    x0_center = x0[nx//2+padding, ny//2]
    x0 = x0 - x0_center # Delay compared with the center

    ## x0 = t[geom]/R => t[SI] = x0*R/c
    dt_grav = x0*ns_Re*1e3/299792458
    max_delay = np.amax(dt_grav[~np.isnan(dt_grav)])



    for i,k in np.ndenumerate(azimuths):
        i = i[0]
        print(" {:.1f}%".format(100*(i-imin)/(imax-imin)),end="\r")
        cap_az = k*np.pi/180

        # reset dt_grav
        dt_grav = x0*ns_Re*1e3/299792458
        cart = np.array([pxf, pyf, pzf])
        on_cap = filter_cap(cart=cart,cap_colat=cap_colat,cap_az=cap_az,cap_aperture=cap_aperture).reshape(nx,ny)
        dt_grav[ on_cap == False ] = np.nan
        t_rotcap = cap_az/(2*np.pi*ns_freq)
        np.savez_compressed(ensure_dir('data/lc/t_cap/{}'.format(dataid))+'/{:03d}.npz'.format(i), t_rotcap + dt_grav)

    # </ get delay in time of arrival>



    # < get photo maks >
    # correct the hotspot shape with the time delays

    # (!) same size as azimuths in current implementation
    dt_exposure = (daz_cap*np.pi/180)/(2*np.pi*ns_freq)
    t_exposure_0 = 0
    t_exposure_f = (az_max*np.pi/180)/(2*np.pi*ns_freq)
    t_exposure = np.arange(t_exposure_0,t_exposure_f,dt_exposure)

    padding = int((2*np.pi*ns_freq*max_delay*180/np.pi)/daz_cap) + 2
    print("(i) Using padding of {} frames".format(padding))


    print("Getting photo masks...")
    t_cap_slice = {}

    for i,t in np.ndenumerate(t_exposure):
        i = i[0]
        # RESTARTING
        #if os.path.exists('data/lc/photo_masks/{}'.format(dataid)+'/{:03d}.npy'.format(i)):
        #    continue
        print(" {:.1f}%".format(100*(i-imin)/(imax-imin)),end="\r")
        load_necessary_t_cap(i,padding,imin,imax,t_cap_slice)
        photo_mask = get_photo_mask(t_cap_slice.values(),t,dt_exposure)
        np.savez_compressed(ensure_dir('data/lc/photo_masks/{}'.format(dataid))+'/{:03d}.npz'.format(i), photo_mask)

    # < / get photo masks >


    # < flux (t) >

    Fluxes = []
    z = data['z'].reshape(nx,ny)
    fratio = (z + 1)**-1
    fratio = np.nan_to_num(fratio)
    di = dj = 0.001
    D = 1

    print("Getting fluxes...")
    for i in range(t_exposure.shape[0]):
        print(" {:.1f}%".format(100*(i-imin)/(imax-imin)),end="\r")
        photo_mask = np.load('data/lc/photo_masks/{}'.format(dataid)+'/{:03d}.npz'.format(i))['arr_0']
        Flux = np.sum( 1/D**2*fratio[photo_mask==True]**3 * 1 * di*dj )
        Fluxes.append(Flux)


    # </ flux (t)>

    np.savetxt(ensure_dir('data/lc/flux')+'/{}.txt'.format(dataid), np.array([t_exposure,Fluxes]).T)

    print("")


if __name__ == "__main__":

    help=r"""
    Computes the light curve of a circular hotspot
    ----
    Usage: python lchotspot_get.py hotspot_dataid ztoa_dataid colatitude angular_radius central_geodesic
        hotspot_dataid: id for saving the generated data for the specific hotspot
        ztoa_dataid: dataid of the GR shift + toa data
        colatitude: colatitude of the hotspot in degrees (from positive z axis)
        angular_radius: angular_radius of the hotspot in degrees
        central_geodesic: (normally = 0) Reference geodesic for the time of arrival
    """

    if len(sys.argv) < 6:
        print(help)
        quit()

    dataid = sys.argv[1]
    toa_filename = sys.argv[2]
    cap_colat = float(sys.argv[3])*np.pi/180
    cap_aperture = float(sys.argv[4])*np.pi/180
    padding = int(sys.argv[5]) # 1 for Mink, 0 for the rest

    metadata = np.loadtxt('data/zt/meta/{}.txt'.format(toa_filename))
    nx = ny = int(metadata[0])
    ns_Re = metadata[5]
    ns_freq = metadata[6]

    print("(i) using nx=ny={}, Re = {}, freq = {}".format(nx,ns_Re,ns_freq))

    main(nx,ny,dataid, ns_freq, ns_Re, toa_filename, cap_colat, cap_aperture, padding)
