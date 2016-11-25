#!/bin/tcsh -f


#----------------------------------------------------
# Script to take the grd file made from SRTM
# data, and convert it to a geographic grid
#
# Note: SRTM data is formed on the WGS-84 ellipsoid
# This is datum 219 in the GMT database.
#---------------------------------------------------

# Form an xyz file in the restricted area around the windfarm and the
# infrasound array
gmt grd2xyz srtm_34_14_pl.grd -R-14.385/-14.353/-7.97/-7.925 > test.xyz

# Just get the xy coordinates (in lon / lat, hence .ll)
awk '{print $1, $2}' test.xyz > test.ll
awk '{print $3}' test.xyz > test.z

gmt mapproject -R -Ju28M/1:1 -F test.ll > test.xy

# Add the heights back onto the xy values.

paste test.xy test.z > ascension.xyz


