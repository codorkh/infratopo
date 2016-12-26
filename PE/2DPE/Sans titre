gmt gmtset FONT_ANNOT_PRIMARY 18p,Helvetica,black
gmt gmtset FONT_LABEL 18p,Helvetica,black
gmt gmtset FONT_TITLE 18p,Helvetica,black
gmt gmtset PS_CHAR_ENCODING Standard+
gmt gmtset FORMAT_GEO_MAP ddd.xxx
#./esri-asc2xyz.pl < ./srtm_23_17/srtm_23_17.asc > srtm_23_17.xyz

#gmt xyz2grd srtm_23_17.xyz -Gsrtm_23_17.grd -R-69.90/-65.1/-24.90/-21.00 -I3s
gmt grdtrack -E-24.5/-69.8/-24.5/-65.2+d -f0x,1y,2f,3f -Getopo2.grd > testing.xydz

awk '{print $3, $4}' testing.xydz > testing.dz

gmt psxy testing.dz -JX15/8 -R0/500/0/7000 -Ba100f50:"Dist (km)":/a1000b500 > tmp.ps
