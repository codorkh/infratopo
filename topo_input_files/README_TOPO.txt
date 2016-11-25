22nd Nov 2016

Aim is to generate a set of test case topography input files that Codor can use
to test his program.

Below, I list the files in each directory. New files can be created using the python (and GMT) scripts provided. However
they will need modification by the user.


2D Cases
--------

1. Flat Case
2. Gaussian Hill - regular x-spacing
3. Gaussian Hill - irregular x-spacing (testing interpolation function in Codor's program)
4. Ascension - 2D profiles between wind farm and the sensors (L configuration)

The first three should be relatively trivial, but we should take care in providing cases that
can test the code.

I have split the files into two directories: 'synthetics' and 'ascension'

Synthetics
----------
flat_and_gaussian_create.py -> python file to create inputs (note you may need to change path definitions if you want to run it)
flat_topo.dat -> An example flat topography file
gauss_200m_hill.dat -> A 2d profile with a 200m high gaussian hill, with a half-width of 600m (see python code for definitions)
gauss_200m_hill_irreg.dat -> Same as gauss_200m_hill.dat, but not regularly spaced x-points

flat_topo_3d.dat -> An example 3D flat topography
gauss3d_200m_hill.dat -> An example 3D file with a Gaussian hill that has half-widths the same in x- and y- directions (of 600m)


Note that I've parameterised the Gaussians in terms of maximum altitude, rather than
Codor's maximum angle

Ascension
---------
topo_notes_4codor.pdf -> Notes (including a map of Ascension) about the topographic profiles
plot_profiles_Lx.py -> A python script that plots profiles, and outputs ASCII 2D profiles
topo_profiles_asc.gmt -> the GMT file to extract the profiles from the SRTM grd file. Note I have not included the grid file, this
	was created in the usual way. The data was downloaded, then the esri-asc2xyz.pl script was used to generate an .xyz file. This
	was then converted to a .grd file using the GMT command:
	gmt xyz2grd srtm_34_14_plout.xyz' -Gsrtm_34_14_pl.grd -R-14.45/-14.27/-8.00/-7.88 -I3c
asc_L1_2d.dat -> 2d profile (created from plot_profiles_Lx.py) for wind turbine to L1 path
asc_L2_2d.dat -> 2d profile for wind turbine to L2 path
asc_L3_2d.dat -> 2d profile for wind turbine to L3 path
asc_L4_2d.dat -> 2d profile for wind turbine to L4 path

srtmgrd_to_metres_grid.csh -> A script that takes the SRTM grid file (not included here) and forms a UTM x,y,z file (with
coordinates in metres, of a small subset of these values). Note that the xyz file has the origin (0,0) at an arbitrary 
latitude/longitude point (14.39W,7.97S) and to use with Codor's programs we will need to transform so that the origin is at
the source (wind turbine).

ascension.xyz -> the x,y,z output from srtmgrd_to_metres_grid.csh
temp_L1_L2_L3_L4_windfarm.dat -> the coordinates (in the same frame of reference as ascension.xyz) for the Ascension Island
infrasound array (I50GB) elements L1,L2,L3,L4 and the windfarm (in that order)
asc_xyz_plot.py -> A python routine to plot a contour plot of the map, and the stations.
ascension_xyz_example_crude.png -> A crude plot (output from asc_xyz_plot.py). Stations are black dots, windfarm is red star
