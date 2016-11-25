A directory to hold input files for Codor's Modelling:

NCPA_canonical_profile_zuvwtdp.dat -> the NCPA 'toy' atmosphere. See the NCPAPROP package for details

Note: I have provided atmospheric velocity profiles up to 200km in altitude. This is probably too much,
the profiles could be truncated accordingly. The profiles are really intended for the long-range propagation
problems - more specific profiles for  local (10's km propagation) may be more appropriate.

NCPA_canonical_profile_zuvwtdp.dat -> the full NCPA toy atmospheric model (see ncpaprop codes)

atmos_2d_adiabatic_sound_sp.dat -> adiabatic sound speed profile (altitude (m), sound speed (m/s))
atmos_2d_eff_sound_sp.dat -> effective sound speed from toy model (dirn chosen in generate_atmos_files.py) (altitude (m), sound speed (m/s))
atmos_2d_isothermal_sound_sp.dat -> isothermal sound speed (temperature chosen in generate_atmos_files.py) (altitude (m). sound speed (m/s))
example_atmospheres.png -> An example figure showing the altitude dependence of the different models
generate_atmos_files.py -> A python script to generate the atmospheric profile files

