"""
File: plot_everything.py

Version
-------
V0.1, 01 Dec 2021, dvk - implement just one plot
V0.2, 04 Dec 2021, dvk - continue

"""
import numpy as np
import matplotlib.pyplot as plt
from core.compute_trajectory import compute_trajectory

# -------------------------------------------------------------------------------------------------
# set variables
# -------------------------------------------------------------------------------------------------
# The following reads the file config.py as if it were part of this code.
exec(compile(source=open("config.py").read(), filename="config.py", mode="exec"))

linestyles = ["solid", "dashed", "dashdot", "dotted"]
legend_text = [f"r = {radius * 100} cm" for radius in r_objects]


initial_mass_objects = []           # This is computed from the radius
initial_ballistic_coeffs = []
final_angle = []              # Impact angle (not really needed for AllBert EinStein)
final_r = []                  # Radius of the final object.

final_angles = []
final_mass = []
final_rs = []

# Initialize the output variables to be a list.
ta_r = []
sa_r = []
ha_r = []
ra_r = []
va_r = []
thetaa_r = []
a2a_r = []
ma_r = []
mva_r = []
zda_r = []

for i, r_object in enumerate(r_objects):

    initial_mass_object = 4. / 3 * np.pi * r_object**3 * rho_object  # mass of projectile sphere in kg
    initial_ballistic_coeff = initial_mass_object / (drag_coeff * np.pi * r_object**2)  # normal formula balistic coef

    ta, sa, ha, ra, va, thetaa, a2a, ma, mva, zda = compute_trajectory(atmo_density, r_planet,
                    g0_planet, a, drag_coeff, heat_transfer, heat_of_ablation, lum_eff, v_init,
                    z_init, initial_mass_object, rho_object, h_obs, dt, array_length)

    initial_mass_objects.append(initial_mass_object)
    initial_ballistic_coeffs.append(initial_ballistic_coeff)    
    final_angle.append(zda[-1])
    final_mass.append(min(ma))
    final_r.append((3. / 4. * ma[-1] / rho_object / np.pi)**0.333)

    ta_r.append(ta)
    sa_r.append(sa)
    ha_r.append(ha)
    ra_r.append(ra)
    va_r.append(va)
    thetaa_r.append(thetaa)
    a2a_r.append(a2a)
    ma_r.append(ma)
    mva_r.append(mva)
    zda_r.append(zda)

final_angles.append(final_angle)
final_rs.append(final_r)


# -------------------------------------------------------------------------------------------------
def single_plot(col_no, row_no, x_label, y_label, x_data, y_data, \
                x_lim=None, y_lim=None, x_log=False, y_log=False):
    """Generic routine to plot one single plot

    Parameters
    ----------
    col_no, row_no : integer
        Column and row number, counting from 0 and from the upper left.

    x_label, y_label : string
        Text to write as label for x- and y axis.

    x_data : numpy array or list
        Array or list of x data to be plotted

    y_data : numpy array or list, or list of numpy array or list
        The y data to be plotted. Can be a list of list (or array of array). Note that the
        number of elements must match the x data.

    x_lim, y_lim : Two-element tuple
        The minimum and maximum values to be plotted on the respective axis.

    Returns
    -------
    A plot. Note that plt.show() needs to be called to show everything.
    
    """
    for i, y in enumerate(y_data):
        axs[col_no, row_no].plot(x_data, y, linestyle=linestyles[i], label=legend_text[i])
        
    axs[col_no, row_no].set_xlabel(x_label)
    axs[col_no, row_no].set_ylabel(y_label)
    axs[col_no, row_no].set_xlim(x_lim)
    axs[col_no, row_no].set_ylim(y_lim)
    axs[col_no, row_no].legend(fontsize=8)
    axs[col_no, row_no].grid()
    if y_log == True:
        axs[col_no, row_no].set_yscale("log")


# -------------------------------------------------------------------------------------------------
# plot everything
# -------------------------------------------------------------------------------------------------
fig, axs = plt.subplots(3, 2, figsize=fig_size, constrained_layout=True)
fig.set_facecolor("whitesmoke")
fig.suptitle(f"AllBert EinStein ablation simulation\n$\\rho$ = {rho_object} kg/m$^3$, " + \
             f"v = {v_init} km/s, $\\alpha$ = {z_init} deg, \n" + \
             f"tau = {lum_eff}, $c_d$ = {drag_coeff} " + \
             f"shape factor = {a}, \nheat transfer coeff = {heat_transfer}, heat of ablation =" + \
             f" {heat_of_ablation} J/kg\n Final masses: {[round(fm, 3) for fm in final_mass]} kg", \
             fontsize=10)

single_plot(0, 0, "time in s", "magnitude", np.array(ta), np.array(mva_r), \
            x_lim=(0, t_max), y_lim=(10, -15))
single_plot(0, 1, "height in km", "magnitude", np.array(ha) / 1000, np.array(mva_r), \
            x_lim=(100, 0), y_lim=(10, -15))
single_plot(1, 0, "time in s", "mass in kg", np.array(ta), np.array(ma_r), \
            x_lim=(0, t_max), y_lim=(1e-5, 1e3), y_log=True)
single_plot(1, 1, "height in km", "mass in kg", np.array(ha) / 1000, np.array(ma_r), \
            x_lim=(100, 0), y_lim=(1e-5, 1e3), y_log=True)
single_plot(2, 0, "time in s", "velocity in km/s", np.array(ta), np.array(va_r) / 1000, \
            x_lim=(0, t_max))
single_plot(2, 1, "height in km", "velocity in km/s", np.array(ha) / 1000, np.array(va_r) / 1000, \
            x_lim=(100, 0))

"""
plt.grid()
plt.title(f"$\\rho$ = {rho_object} kg/m$^3$, shape factor = {a}, c$_d$ = {drag_coeff}, " + \
          f"tau = {lum_eff}\n heat transfer coefficent = {heat_transfer}\n" +
          f"heat of ablation = {heat_of_ablation} J/kg",
          fontsize=8)
plt.xlim([100, 0])
for i, r_object in enumerate(r_objects):
    plt.semilogy([x / 1000 for x in ha_r[i]], ma_r[i], label=f"r = {r_object * 100} cm")
plt.xlabel('height in km')
plt.ylabel('mass in kg')
plt.legend()
plt.show()
"""

"""
plt.grid()
plt.title(f"tau = {lum_eff}")
plt.ylim([10, -15])
for i, r_object in enumerate(r_objects):
    plt.plot(ta_r[i], mva_r[i], label=f"r = {r_object * 100} cm", linestyle=linestyles[i])
plt.xlabel("time in s")
plt.ylabel("magnitude")
plt.legend()
plt.show()


plt.grid()
plt.title(f"rho = {rho_object} kg/m$^3$, shape factor = {a}, c$_d$ = {drag_coeff}, " + \
          f"tau = {lum_eff}\n heat transfer coefficent = {heat_transfer}\n" +
          f"heat of ablation = {heat_of_ablation} J/kg",
          fontsize=8)

for i, r_object in enumerate(r_objects):
    plt.semilogy(ta_r[i], ma_r[i], label=f"r = {r_object * 100} cm")
plt.xlabel('time in s')
plt.ylabel('mass in kg')
plt.legend()
plt.show()
"""

plt.show()

# ================================================ eof ============================================
