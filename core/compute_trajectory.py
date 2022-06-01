#==================================================================================================
# file compute_trajectory.py
#
# Computing the behavior of the reentry of a compact object for the AllBert EinStein mission
# Based on Kennewell - A flight of a meteor, Sky and Telescope (1987).
# A simple meteoroid/asteroid atmospheric entry computation
#
# V1.0, 2017-10-28, dvk: based on metpath_mars - plot angle versus mass
# V1.1, 2017-10-29, dvk: update to run 3-D parameter space
# V1.2, 2017-10-31, dvk: updated to plot as function of final radius
# V2.0, 2021-11-27, dvk: apply to AllBert EinStein - added doc strings
#==================================================================================================

#--------------------------------------------------------------------------------------------------
# import libraries
#--------------------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt


def compute_trajectory(atmo_density, r_planet, g0_planet, a, drag_coeff, heat_transfer,
                       heat_of_ablation, lum_eff, v_init, z_init, mass_init, rho_object,
                       h_obs, dt, array_length):
    """
    Computing the trajectory of an individual meteoroid.

    Parameters
    ----------
    atmo_density: nparray of float
        Density profile of the atmosphere in heights of 0 km, 10 km, ..., given in tbd.

    r_planet: float
        Radius of the planet

    g0_planet: float
        Gravitational constant of the planet, use 9.81 m/s2 for the Earth.

    a: float
        Shape factor of meteoroid - 1.2 for a sphere.

    drag_coeff: float
        Drag coefficient c_d: 0.05 (drop), 0.47 (sphere), 1.15 (short cylinder).

    heat_transfer: float
        Heat transfer coefficient

    heat_of_ablation: float
        Heat of ablation in J/kg. Depends on material, between 1e5 and 1e7 J/kg.

    lum_eff: float
        Luminous efficiency of the meteoroid (constant over flight path)

    v_init: float
        Initial velocity in km/s.

    z_init: float
        Inital angle w.r.t. horizontal, in degrees.

    mass_init: float
        Initial mass of meteoroid in kg.

    rho_object: float
        Density of the object in kg/m3.
        
    h_obs: float
    	The distance of the observer in m.
    	
    dt: float
    	Time step in s.
    	
    array_length: integer
    	Length of the output arrays. We initialise everything with np.nan, so that
    	they are all the same length.

    Returns
    -------
    ta: nparray of float
        An array with all time stamps of the computed points, in s.

    sa: nparray of float

    ha: nparray of float

    ra: nparray of float

    va: nparray of float
        An array with the veloicty of the object as a function of time

    thetaa: nparray of float
        An array with the flight path angle (w.r.t. horizontal) as a function of time

    a2a:
        tbd
        
    ma:
        tbd
        
    mva:
        tbd
        
    zda:
        tbd
        

    To do:
    ------
    - Make the lum_eff a function of e.g. velocity.
    - No fragmentation included.
        
    """
    # Helper function
    def init_array():
        x = np.empty(array_length)
        x[:] = np.nan
        return x


    sg = heat_transfer / (2 * drag_coeff * heat_of_ablation) #????

    h = 120000.0                         # initial altitude in m
    s = 0.0                              # initial ground track
    t = 0.0                              # initial time
    fm = 1                               # output will be printed every fm time steps
    mvmax = 9.0                          # faintest magnitude to be printed

    f1 = fm
    m = mass_init                        # intial mass is set
    v = v_init * 1000.                   # convert velocity to m/s
    zr = z_init * np.pi / 180.           # convert zenith angle to radians
    v_vertical = -v * np.cos(zr)         # initial velocity in vertical direction
    v_horizontal = v * np.sin(zr)        # initial velocity in horizontal direction

    first_time = 1              # flag to indicate that we run the first time

    # Initalize arrays for output
    ta = init_array()
    sa = init_array()
    ha = init_array()
    ra = init_array()
    va = init_array()
    thetaa = init_array()
    a2a = init_array()
    ma = init_array()
    mva = init_array()
    zda = init_array()

    #----------------------------------------------------------------------------------------------
    # do the computation
    #----------------------------------------------------------------------------------------------
    for step in np.arange(0, array_length):
        # find density of atmosphere from lookup table. The table is in steps of 10 km.
        i = int(np.floor(h/10000.))
        if i > 15:                                            # too high for noticeable atmosphere
            i = 15
        fr = (h / 10000.) - i

        # interpolate current density, 1 is the number one (not al 'el')
        rho_air =  10**(atmo_density[i] * (1 - fr) + fr * atmo_density[i + 1])

        # calculate rates of change
        a2 = drag_coeff * a * rho_air * v * v / (m * rho_object * rho_object)**0.33333
        gv = g0_planet / (1 + h / r_planet)**2
        av = -gv - a2 * v_vertical / v + v_horizontal * v_horizontal / (r_planet + h)
        ah = -a2 * v_horizontal / v - v_vertical * v_horizontal / (r_planet + h)
        m0 = sg * m * v * a2                                            # mass loss
        zd = np.arctan2(v_vertical, v_horizontal) * 180. / np.pi        # compute new zenith angle
        if (f1 == fm):
            # calculate visual magnitude via empirical formula
            luminous_energy = 0.5 * v * v * m0 * lum_eff * 1e10 / (h_obs * h_obs)
            if luminous_energy < 0.1:
                mv = np.nan
            else:
                mv = 6.8 - 2.086 * np.log10(luminous_energy)

            if mv < mvmax:
                if first_time == 1:
                    first_time = 0
                    t = 0.
                
                # print output
                # print('%7.1f  %8.1f   %6.1f  %6.1f  %7.1f   %6.1f %5.1f  %6.1f' %
                # (t, s/1000., h/1000., v/1000., a2, 100. * m / mass_init, mv, zd))
                
            # assign to output arrays
            ta[step] = t
            sa[step] = s
            ha[step] = h
            ra[step] = h + r_planet
            va[step] = v
            a2a[step] = a2
            ma[step] = m
            mva[step] = mv
            zda[step] = zd
            # end if mv < mvmax
            f1 = 0
        # endif f1 == fm
    
        # increment values
        t = t + dt
        s = s + v_horizontal * dt * r_planet / (r_planet + h)  # ground track
        h = h + v_vertical * dt                                # height
        m = m - m0 * dt                                        # new mass
        v_vertical = v_vertical + av * dt
        v_horizontal = v_horizontal + ah * dt
        v = np.sqrt(v_horizontal * v_horizontal + v_vertical * v_vertical)
        f1 = f1 + 1
        
        # exit criteria
        if h < 0:
        	# height in m
        	print(f"Initial mass: {mass_init} kg - height below 0. Final mass: {m} kg")
        	break
        
        if np.sqrt(v_vertical**2 + v_horizontal**2) < 55:
        	# velocity in m/s
        	print("velocity below 200 km/h")
        	break
        	
        if mv > mvmax:
            # magnitude below limit
            print("magnitude below limit.")
            break

    # Compute polar angle from ground flight path; '%' is the modulus, i.e. if value larger than
    # 360 deg it starts again at 0 deg.
    thetaa = [(s / r_planet) % (2 * np.pi) for s in sa]

    return ta, sa, ha, ra, va, thetaa, a2a, ma, mva, zda
