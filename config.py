"""
File: config.py
Configuration values for meteor computation

"""

### Geometry
v_init = 8.0         # Speed of projectile in km/s - 7.8112 for h = 155 km @ Earth.
                      # 3.54 km/s for r_Mars; around 9.0 for AllBert Einstein
z_init = 85.0         # Zenith angle of projectile path in degrees.
h_obs = 100000.0      # observers distance in m to calculate apparent magnitude this should really
                      # be computed at every step! The original paper just used 'h', i.e. the
                      # observer just below the object.

fig_size = (8, 10)

dt = 0.025             # time step in s.
t_max = 300            # Maximum time to consider in seconds.
array_length = int(t_max / dt)

# Define Earth atmosphere - numbers are log of density in kg/m3 for 0 km, 10 km, 20 km, ...
atmo_density = [0.09, -0.38, -1.05, -1.74, -2.39, -2.98, -3.50, -4.07, -4.72, -5.45, \
                -6.30, -7.00, -7.62, -7.97, -8.32, -8.67, -8.81]
g0_planet = 9.81      # Earth acceleration in m/s2.
r_planet = 6371000.0  # Earth radius in m.

### Meteoroid parameters
rho_object = 2800.0   # density of projectile in kg/m3 - stone: 2800.0; iron: 7800.0.
a = 1.2               # shape factor (1.2 for sphere).
drag_coeff = 0.5      # drag coefficient c_d: 0.05 (drop), 0.47 (sphere), 1.15 (short cylinder).
heat_transfer = 0.02  # heat transfer coefficient, between 0.05 and 0.6.
heat_of_ablation = 500000        # heat of ablation in J/kg, depends on material, 1e5..1e7 J/kg.

lum_eff = 0.1        # luminous efficency, 0.0002 (dust) to 0.002 (stone, iron). According to
                      # Ott et al. (2021), values between 0.1 and 0.01 should be used.

r_objects = np.array([0.025, 0.05, 0.10, 0.20])       # Radius of object in m
