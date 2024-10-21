# AllBertEinStein_Thesis_code

## Introduction

This repository serves as a collection of all the calculations, graphs, and design considerations made during the thesis work, developed using MATLAB and Python scripts. All the code is based on freely available libraries.

## Structure

The repository is organized into two main folders:

1. **DRAMA-SARA_Reader**: MATLAB scripts to read the DRAMA SARA output and put them in plots
2. **FRIPON_Reader**: Contains all the MATLAB functions use to fit the meteors we tested our code with
3. **LongMarch3B_Reader**: Contains all the MATLAB functions use to fit the LongMarch3B reentry data
4. **Meteor_code-AllBert_EinStein**: Contains all the Python code reated by Detlef and modified by Max during the thesis on the AllBert EinStein mission
5. **SCARAB-Data_Reader**: MATLAB scripts to read the SCARAB output and put them in plots

### Meteor_code-AllBert_EinStein/code: compute_trajectory.py

- Capable of simulating reentry for objects ranging from cm to m in size.
- Initial altitudes that can reach up to 150 km.
- Velocity ranging from 7.5 km/s to 72 km/s.
- The code will terminate based on several conditions:
  - Reaching the maximum computation time.
  - The object reaches the ground.
  - The object reaches terminal velocity.
  - The object's final magnitude is below a set minimum observable limit.
  
- Provides output in terms of longitude, latitude, azimuth, and elevation based on the local observer.
- Incorporates the Earth’s Coriolis effect into the equations of motion.

The Python code does not propagate the initial orbit, meaning the initial position and heading angle must be selected carefully; otherwise, the final position may be invalid. Propagation based on the initial conditions in the geodetic reference frame will allow simulation at higher altitudes, but it was deemed superfluous for the scope of this thesis.

### Meteor_code-AllBert_EinStein: plot_everything.py

This script is the main driver for plotting various results computed by compute_trajectory.py. It generates a variety of subplots that provide detailed insights into the simulated reentry trajectory base also on the visibility map.

## Contributors

- Detlef Koschny: Original developer of the Python Meteor code.
- Maximilian Vovk: Developed the MATLAB scripts and modified the Python code for the purposes of this thesis.

Special thanks to all the team members and contributors who made this thesis possible.
