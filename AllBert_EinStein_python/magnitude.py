"""
File: magnitude.py

Version
-------
01 Dec 2021, dvk

"""

import numpy as np
import matplotlib.pyplot as plt
from core.compute_trajectory import compute_trajectory

# -------------------------------------------------------------------------------------------------
# set variables
# -------------------------------------------------------------------------------------------------
# The following reads the file config.py as if it were part of this code.
exec(compile(source=open("config.py").read(), filename="config.py", mode="exec"))

h_obs = 140000.0      # observers distance in m to calculate apparent magnitude this should really
                      # be computed at every step! The original paper just used 'h', i.e. the
                      # observer just below the object.
                      

mass_objects = []             # This is computed from the radius
ballistic_coeffs = []
final_angle = []              # Impact angle (not really needed for AllBert EinStein)
final_r = []                  # Radius of the final object.

final_angles = []
final_mass = []
final_rs = []

for r_object in r_objects:
        
    mass_object = 4. / 3 * np.pi * r_object ** 3 * rho_object  # mass of projectile in kg
    ballistic_coeff = mass_object / (drag_coeff * np.pi * r_object**2)

    ta, sa, ha, ra, va, thetaa, a2a, ma, mva, zda = compute_trajectory(atmo_density, r_planet,
                    g0_planet, a, drag_coeff, heat_transfer, heat_of_ablation, lum_eff, v_init,
                    z_init, mass_object, rho_object, h_obs, dt, array_length)

    mass_objects.append(mass_object)
    ballistic_coeffs.append(ballistic_coeff)    
    final_angle.append(zda[-1])
    final_mass.append(ma[-1])
    final_r.append((3. / 4. * ma[-1] / rho_object / np.pi)**0.333)

    plt.plot(ta, mva)

final_angles.append(final_angle)
final_rs.append(final_r)


# -------------------------------------------------------------------------------------------------
# plot something
# -------------------------------------------------------------------------------------------------
plt.ylim([20, -20])
plt.grid()
plt.xlabel('time in s')
plt.ylabel('magnitude')
plt.show()

# ================================================ eof ============================================
