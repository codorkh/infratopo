#set terminal x11
#set terminal post eps rounded color dashed dl 4 "Times-Bold" 35
set terminal png enhanced font "Cambria:Bold, 16" size 1024, 1024
set output "LPg_SnT_2017.png"
#set output "LPg_L4.eps"

file11="LPg_L1_f1.dat"
file12="LPg_L1_f2.dat"
file13="LPg_L1_f3.dat"
file14="LPg_L1_f4.dat"
topo11="TOPO_L1_f1.dat"
topo12="TOPO_L1_f2.dat"
topo13="TOPO_L1_f3.dat"
topo14="TOPO_L1_f4.dat"

file21="LPg_L2_f1.dat"
file22="LPg_L2_f2.dat"
file23="LPg_L2_f3.dat"
file24="LPg_L2_f4.dat"
topo21="TOPO_L2_f1.dat"
topo22="TOPO_L2_f2.dat"
topo23="TOPO_L2_f3.dat"
topo24="TOPO_L2_f4.dat"

file31="LPg_L3_f1.dat"
file32="LPg_L3_f2.dat"
file33="LPg_L3_f3.dat"
file34="LPg_L3_f4.dat"
topo31="TOPO_L3_f1.dat"
topo32="TOPO_L3_f2.dat"
topo33="TOPO_L3_f3.dat"
topo34="TOPO_L3_f4.dat"

file41="LPg_L4_f1.dat"
file42="LPg_L4_f2.dat"
file43="LPg_L4_f3.dat"
file44="LPg_L4_f4.dat"
topo41="TOPO_L4_f1.dat"
topo42="TOPO_L4_f2.dat"
topo43="TOPO_L4_f3.dat"
topo44="TOPO_L4_f4.dat"

###############

set size 1,1
set multiplot layout 1,4
set style func linespoints
set border lw 2
#set label "Transmission loss at z = 0 along paths L1-4" font "Times,22" at screen 0.50,0.82 center

###############
#PLOT1
###############

set tmargin at screen 0.80; set bmargin at screen 0.55
set lmargin at screen 0.10; set rmargin at screen 0.50
set origin 0.0,0.5
set size 0.5,0.5
set key out horiz
set key samplen 1
set key spacing -10
set key width -20
set key at screen 1.46,0.83
set xlabel ''
set format x "%g"
set xrange [0:4500]
set xtics format ''
set xtics scale 0
set xtics ("" 0, "" 500, "" 1000, "" 1500, "" 2000, "" 2500, "" 3000, "" 3500, "" 4000, "" 4500, "" 5000) nomirror

set ylabel "TL (dB)" offset 1,0
set format y "%g"
set yrange [-70:0]
set ytics format '%g'
set ytics ("" -120, "" -110, "" -100, "" -90, "" -80, "" -70, "-60" -60, "-50" -50, "-40" -40, "-30" -30, "-20" -20, "-10" -10, "0" 0) nomirror

set y2label ''
set y2range [100:450]

plot topo11 using 1:2 axes x1y2 with lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "Terrain", \
file11 using 1:2 with line lt -1 lw 4 lc 7 title "1.90 Hz", \
file12 using 1:2 with line lt -1 lw 4 lc 1 title "3.79 Hz", \
file13 using 1:2 with line lt -1 lw 4 lc 3 title "5.69 Hz", \
file14 using 1:2 with line lt -1 lw 4 lc 4 title "7.59 Hz"

set label "{/:Bold L1}" at screen 0.45,0.77 font "Cambria, 16"

###############
###PLOT 2
###############
set tmargin at screen 0.80; set bmargin at screen 0.55
set lmargin at screen 0.50; set rmargin at screen 0.90
set origin 0.5,0.5
set size 0.5,0.5
unset key
unset ylabel

set xtics format ''
set ytics format '' nomirror
set ytics scale 0
set ytics ("" -120, "" -110, "" -100, "" -90, "" -80, "" -70, "" -60, "" -50, "" -40, "" -30, "" -20, "" -10, "" 0) nomirror

set y2label "Elevation (m)"
set y2tics ("" 100, "150" 150, "200" 200, "250" 250, "300" 300, "350" 350, "400" 400, "" 450, "" 500) nomirror
set y2range [100:450]
set xtics ("" 0, "" 500, "" 1000, "" 1500, "" 2000, "" 2500, "" 3000, "" 3500, "" 4000, "" 4500, "" 5000) nomirror
set xrange [0:4500]
plot topo21 using 1:2 axes x1y2 with lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "terr.", \
file21 using 1:2 with line lt -1 lw 4 lc 7 title "TL(z=0)" , \
file22 using 1:2 with line lt -1 lw 4 lc 1 title "TL(z=0)", \
file23 using 1:2 with line lt -1 lw 4 lc 3 title "TL(z=0)", \
file24 using 1:2 with line lt -1 lw 4 lc 4 title "TL(z=0)"

set label "{/:Bold L2}" at screen 0.85,0.77 font "Cambria,16"

###############
#PLOT 3
###############

set tmargin at screen 0.55; set bmargin at screen 0.30
set lmargin at screen 0.10; set rmargin at screen 0.50
set origin 0.0,0.0
set size 0.5,0.5
unset y2label
unset y2tics

set xlabel 'range (km)'
set format x "%g"
set xrange [0:4500]
set xtics format ''
set xtics scale 0
set xtics ("0" 0, "0.5" 500, "1" 1000, "1.5" 1500, "2" 2000, "2.5" 2500, "3" 3000, "3.5" 3500, "4" 4000, "" 4500, "" 5000) nomirror

set ylabel "TL (dB)" offset 1,0
set format y "%g"
set yrange [-70:0]
set ytics format '%g'
set ytics ("" -120, "" -110, "" -100, "" -90, "" -80, "" -70, "-60" -60, "-50" -50, "-40" -40, "-30" -30, "-20" -20, "-10" -10, "0" 0) nomirror

set y2label ''
set y2range [100:450]

plot topo31 using 1:2 axes x1y2 with lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "terrain", \
file31 using 1:2 with line lt -1 lw 4 lc 7 title "TL(z=0)", \
file32 using 1:2 with line lt -1 lw 4 lc 1 title "TL(z=0)", \
file33 using 1:2 with line lt -1 lw 4 lc 3 title "TL(z=0)", \
file34 using 1:2 with line lt -1 lw 4 lc 4 title "TL(z=0)"

set label "{/:Bold L3}" at screen 0.45,0.52 font "Cambria,16"
set label "{/:Bold L4}" at screen 0.85,0.52 font "Cambria,16"

###############
#PLOT 4
###############

set tmargin at screen 0.55; set bmargin at screen 0.30
set lmargin at screen 0.50; set rmargin at screen 0.90
set origin 0.5,0.0
set size 0.5,0.5
#set key box
#set key at screen 0.87,0.45
set xlabel "range (km)"
unset ylabel
set ytics format '' nomirror
set ytics scale 0
set ytics ("" -120, "" -110, "" -100, "" -90, "" -80, "" -70, "" -60, "" -50, "" -40, "" -30, "" -20, "" -10, "" 0) nomirror

set y2label "Elevation (m)"
set y2tics ("" 100, "150" 150, "200" 200, "250" 250, "300" 300, "350" 350, "400" 400, "" 450, "" 500) nomirror
set y2range [100:450]
set xrange [0:4500]
plot topo41 using 1:2 axes x1y2 with lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "terrain", \
file41 using 1:2 with line lt -1 lw 4 lc 7 title "1.90 Hz" , \
file42 using 1:2 with line lt -1 lw 4 lc 1 title "3.79 Hz", \
file43 using 1:2 with line lt -1 lw 4 lc 3 title "5.69 Hz", \
file44 using 1:2 with line lt -1 lw 4 lc 4 title "7.59 Hz"

unset multiplot
