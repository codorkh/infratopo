#!/bin/tcsh -f

# A GMT script to pull out an Alpine profile from the
# ETOPO grid.

gmt grdtrack -E11.56/44.57/8.88/49.17+d -f0x,1y,2f,3f -G/sharedprograms/data/srtm30/w020n90.nc > alp.xydz


