# ==================================================================================================
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
# V2.1, 2022-05-21, MAX: added variable observer position
# V2.2, 2022-07-25, MAX: added radius data and variable Cd with Kn actived for Cd<0
# V2.3, 2022-08-27, MAX: added longitude latitude data
# V2.4, 2022-09-20, MAX: added Elevation Azimuth data & Change equation of motions add Earth Coriolis and Eötvös effect
# ==================================================================================================

# --------------------------------------------------------------------------------------------------
# import libraries
# --------------------------------------------------------------------------------------------------
import numpy as np

def compute_trajectory(atmo_density, r_planet, g0_planet, a, drag_coeff, heat_transfer,
                       heat_of_ablation, lum_eff, v_init, z_init, mass_init, rho_object,
                       obs_altitude, dt, array_length, head_ang=90 ,ast_init_lat=0, 
                       ast_init_lon=0, OBS_lat_north=0, OBS_lon_west=0):
    
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
        
    obs_altitude: float
    	The altitude of the observer in m.
    	
    dt: float
    	Time step in s.
    	
    array_length: integer
    	Length of the output arrays. We initialise everything with np.nan, so that
    	they are all the same length.
    
    head_ang: float
        heading angle direction on earth horizontal plane
        
    ast_init_lat: float
        initial latitude of the asteroid in deg
        
    ast_init_lon: float
        initial longitude of the asteroid in deg
        
    OBS_lat_north: float
        latitude of the observer
    
    OBS_lon_west: float
        longitude of the observer

    Returns
    -------
    ta: nparray of float
        An array with all time stamps of the computed points, in s.

    sa: nparray of float
        An array with all the horizontal movement simulating the ground track

    ha: nparray of float
        An array with all the altitude of the object

    ra: nparray of float
        An array with all the altitude of the object plus the earth radius

    va: nparray of float
        An array with the veloicty of the object as a function of time

    thetaa: nparray of float
        An array with the flight path angle (w.r.t. horizontal) as a function of time

    a2a: nparray of float
        Drag deceleration of the object
        
    ma: nparray of float
        Current mass of the object
        
    mva: nparray of float
        visual magnitude of the object
        
    zda: nparray of float
        Zenith angle of the velocity direction
        
    h_dista: nparray of float
        distance of the oserver from the asteroid
        
    lona: nparray of float
        longitude of the asteroid in deg from 180 to -180
    
    lata: nparray of float
        latitude of the asteroid in deg from 90 to -90
    
    azimuth: nparray of float
        azimuth angle in deg from the observer point of view
        
    elevation: nparray of float
        elevation angle in deg from the observer point of view
        
    """

    # Helper function
    def init_array():
        x = np.empty(array_length)
        x[:] = np.nan
        return x
    
    # Cd coef. base on DRAMA-3.0.0-Final-Report.pdf
    def drag_coef_Kn(h,r,A=1.2): 
        # Mean free path of the 1976 Standard Atmosphere sectioned valid until 200 km
        if h > 0 and h < 11000:
            mean_free_path = 6.213e-08 * np.exp(0.0001027*h)
        elif h > 11000 and h < 20000:
            mean_free_path = 3.485e-08 * np.exp(0.0001569*h)
        elif h > 20000 and h < 47000:
            mean_free_path = 3.892e-08 * np.exp(0.0001533*h)
        elif h > 47000 and h < 51000:
            mean_free_path = 1.491e-07 * np.exp(0.0001243*h)
        elif h > 51000 and h < 61000:
            mean_free_path = 2.022e-07 * np.exp(0.0001181*h)
        elif h > 61000 and h < 71000:
            mean_free_path = 8.344e-08 * np.exp(0.0001322*h)
        elif h > 71000 and h < 81000:
            mean_free_path = 1.992e-08 * np.exp(0.0001519*h)
        elif h > 81000 and h < 120000:
            mean_free_path = 2.626e-09 * np.exp(0.0001762*h)
        elif h > 120000 and h < 160000:
            mean_free_path = 3.582e-09 * np.exp(0.0001739*h)
        else:
            mean_free_path = 5.357e-09 * np.exp(0.0001717*h)
            
        Kn = mean_free_path / (2 * r) #Kn for r*2 diameter of a sphere
        
        # Cd coef. from DRAMA-3.0.0-Final-Report.pdf
        # pag.33 Figure 7-5: Drag coefficient of a sphere 
        # and correction coeficient (A/1.2)**(1/5) by adapting Figure 7-6
        Cd = ((0.925 + np.arctan(Kn*35)**5 * 0.113) * (A/1.2)**(1/5))/2
        return Cd

    h = 120000.0  # initial altitude in m
    s = 0.0  # starting position of the asteroid ground track
    side = 0.0 # acrostrack starting position
    t = 0.0  # initial time
    fm = 1  # output will be printed every fm time steps
    mvmax = 9.0  # faintest magnitude to be printed

    f1 = fm
    m = mass_init  # intial mass is set
    v = v_init * 1000.  # convert velocity to m/s
    zr = z_init * np.pi / 180.  # convert zenith angle to radians
    v_vertical = -v * np.cos(zr)  # initial velocity in vertical direction
    v_horizontal = v * np.sin(zr)  # initial velocity in horizontal direction
    v_side = 0; # initial velocity in acrostrack
    radius=((m/rho_object)*3/(4*np.pi))**(1/3) # radius of the object
    om_earth=(2*np.pi)/(23.934472*3600) # angular speed of the planet
    head_ang_init=head_ang # initial headding angle
    h_side_ang=0           # initial angle between horizontal and side
    
    # latitude longitude data processing
    ast_lon = ast_init_lon
    ast_lat = ast_init_lat
    
    #in case OBS_lon_west OBS_lat_north are in deg
    OBS_lon = OBS_lon_west # observer longitude
    OBS_lat = OBS_lat_north # observer latitude
    
    Up_dist = (h - obs_altitude)
    North_dist =  (h+r_planet)*(ast_init_lon-OBS_lon) * np.pi / 180 
    West_dist =  (h+r_planet)*(ast_init_lat-OBS_lat) * np.pi / 180 
    h_dist = np.sqrt( (West_dist)**2 + (North_dist)**2 + (Up_dist)**2 ) # distance respec to the observer
    theta = 0
    azim = ( np.arctan2(North_dist , West_dist) ) * 180 / np.pi
    elev = ( np.arctan2(Up_dist , np.sqrt(West_dist**2 + North_dist**2)) ) * 180 / np.pi
    
    first_time = 1  # flag to indicate that we run the first time

    # Initalize arrays for output
    ta = init_array()
    sa = init_array()
    ha = init_array()
    ra = init_array()
    va = init_array()
    a2a = init_array()
    ma = init_array()
    mva = init_array()
    zda = init_array()
    h_dista = init_array()
    thetaa = init_array()
    lona = init_array()
    lata = init_array()
    UpDista = init_array()
    NorthDista = init_array()
    WestDista = init_array()
    azimuth = init_array()
    elevation = init_array()
    
    
    active_Kn_drag=0
    if drag_coeff<0: #flag to use a variable drag with Kn
        active_Kn_drag=1

    # ----------------------------------------------------------------------------------------------
    # do the computation
    # ----------------------------------------------------------------------------------------------
    for step in np.arange(0, array_length):
        # find density of atmosphere from lookup table. The table is in steps of 10 km.
        i = int(np.floor(h / 10000.))
        if i > 15:  # too high for noticeable atmosphere
            i = 15
        fr = (h / 10000.) - i

        # interpolate current density, 1 is the number one (not al 'el')
        rho_air = 10 ** (atmo_density[i] * (1 - fr) + fr * atmo_density[i + 1])
        
        if active_Kn_drag==1 and a==1.2:
            drag_coeff=drag_coef_Kn(h,radius)
        elif active_Kn_drag==1:
            drag_coeff=drag_coef_Kn(h,radius,a)


        # calculate rates of change
        a2 = drag_coeff * a * rho_air * v * v / (m * rho_object * rho_object) ** 0.33333  # deceleration same formula
        gv = g0_planet / (1 + h / r_planet) ** 2  # gravity at that height       
        # Gravity + Drag vertical + centrifugal + Coriolis Earth + Coriolis Earth
        av = -gv - a2 * v_vertical / v * np.cos(h_side_ang * np.pi / 180) + v_horizontal * v_horizontal / (r_planet + h) + 2*om_earth*np.cos( ast_lat * np.pi / 180)*v_horizontal*np.cos(head_ang * np.pi / 180) + 2* om_earth *np.cos( ast_lat * np.pi / 180 ) * v_side * np.sin( head_ang * np.pi / 180 )
        # Drag horizontal + Coriolis + Eötvös effect + Coriolis Earth
        ah = -a2 * v_horizontal / v - v_vertical * v_horizontal / (r_planet + h) - 2*om_earth*np.cos( ast_lat * np.pi / 180)*v_vertical*np.cos(head_ang * np.pi / 180) - 2* om_earth*v_side*np.sin( ast_lat * np.pi / 180 ) 
        # Drag horizontal + Coriolis Earth + Coriolis Earth
        aside= -a2 * v_horizontal / v * np.sin(h_side_ang * np.pi / 180) + 2*om_earth * v_horizontal * np.sin(ast_lat * np.pi / 180) - 2*om_earth * v_vertical * np.cos(ast_lat * np.pi / 180) * np.sin(head_ang * np.pi / 180) 
        m0 = heat_transfer * a * rho_air * v * v * v * ((m / rho_object)**0.66667) / (2 * heat_of_ablation)
        zd = np.arctan2(v_vertical, v_horizontal) * 180. / np.pi  # compute new zenith angle
        if (f1 == fm):
            # calculate visual magnitude via empirical formula
            luminous_energy = 0.5 * v * v * m0 * lum_eff * 1e10 / (h_dist * h_dist)
            if luminous_energy < 0.1:
                mv = np.nan
            else:
                mv = 6.8 - 1.086 * np.log(luminous_energy)  # magnitude before 2.086 but in -1.086 kampbel-brown & koshny

            if mv < mvmax:
                if first_time == 1:
                    first_time = 0
                    t = 0.

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
            thetaa[step] = theta
            h_dista[step] = h_dist
            lona[step] = ast_lon
            lata[step] = ast_lat
            UpDista[step] = Up_dist
            NorthDista[step] = North_dist
            WestDista[step] = West_dist
            azimuth [step] = azim
            elevation [step] = elev
            # end if mv < mvmax
            f1 = 0
        # endif f1 == fm

        # increment values
        t = t + dt
        s = s + v_horizontal * dt * r_planet / (r_planet + h)  # ground track
        h = h + v_vertical * dt                                # height
        side = side + v_side * dt                              # acrosstrack position
        ast_lat=180/np.pi*( (s * np.sin( head_ang * np.pi / 180)) / (h+r_planet) ) + ast_init_lat
        ast_lon=180/np.pi*( (s * np.cos( head_ang * np.pi / 180)) / (h+r_planet) ) + ast_init_lon
        m = m - m0 * dt  # new mass
        radius=((m/rho_object)*3/(4*np.pi))**(1/3) # new radius
        v_vertical = v_vertical + av * dt
        v_horizontal = v_horizontal + ah * dt
        v_side = v_side + aside * dt
        v = np.sqrt(v_horizontal * v_horizontal + v_vertical * v_vertical + v_side * v_side) 
        Up_dist = (h - obs_altitude)
        North_dist =  (h+r_planet)*(ast_lon-OBS_lon) * np.pi / 180 
        West_dist =  (h+r_planet)*(ast_lat-OBS_lat) * np.pi / 180 
        h_dist = np.sqrt( (West_dist)**2 + (North_dist)**2 + (Up_dist)**2 ) # distance respec to the observer
        # Compute polar angle from ground flight path; '%' is the modulus, i.e. if value larger than
        # 360 deg it starts again at 0 deg.
        theta = (s / r_planet) % (2 * np.pi)
        h_side_ang = np.arctan2(side,s) * 180 / np.pi
        head_ang = head_ang_init - h_side_ang
        azim = ( np.arctan2(North_dist , West_dist) ) * 180 / np.pi
        elev = ( np.arctan2(Up_dist , np.sqrt(West_dist**2 + North_dist**2)) ) * 180 / np.pi
        f1 = f1 + 1

        # exit criteria
        if h < 0: # height in m
            # touch the ground
            print(f"Initial mass: {mass_init} kg - height below 0. Final mass: {m} kg")
            break

        if np.sqrt(v_vertical ** 2 + v_horizontal ** 2) < 55: # 200 km/h in m/s
            # reach terminal velocity
            print("velocity below 200 km/h")
            break
        
        if first_time == 0:  # wait to reach the maximum
            if mv > mvmax:
                # magnitude below limit
                print("magnitude below limit.")
                break
        
        if first_time == 0:  # need two entries for distance
            if elev < 1 and h_dista[step]-h_dista[step-1]>0:
                # if it is too low in the horizon and is getting away from us
                print("too close to the horizon.")
                break
        
    return ta, sa, ha, ra, va, thetaa, a2a, ma, mva, zda, h_dista, lona, lata, azimuth, elevation
    # UpDista, NorthDista, WestDista
