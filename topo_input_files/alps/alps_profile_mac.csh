#!/bin/tcsh -f

# A GMT script to pull out an Alpine profile from the
# ETOPO grid.

# Uncomment next line to generate the xydz line 
# gmt grdtrack -E11.56/44.57/8.88/49.17+d -f0x,1y,2f,3f -G/sharedprograms/data/srtm30/w020n90.nc > alp.xydz

awk '{print $3,$4}' alp.xydz > alp.dz

./alp_profile.py
