#set terminal x11
#set terminal post eps rounded color dashed dl 4 "Times-Bold" 35
set terminal png enhanced font "Cambria:Bold, 16" size 1024, 1024
set output "LPg_Validation.png"
#set output "LPg_L4.eps"

file11="Atten_Case1_10.dat"
file12="Atten_Case1_100.dat"
file21="Atten_Case2_10.dat"
file22="Atten_Case2_100.dat"

result11="LPg_W0_A1_f1.dat"
result12="LPg_W0_A1_f2.dat"
result21="LPg_W0_A2_f1.dat"
result22="LPg_W0_A2_f2.dat"

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
#set xlabel ''
#set format x "%g"
#set xrange [0:10000]
#set xtics format ''
#set xtics scale 0
#set xtics ("" 0, "" 500, "" 1000, "" 1500, "" 2000, "" 2500, "" 3000, "" 3500, "" 4000, "" 4500, "" 5000) nomirror

#set ylabel "TL (dB)" offset 1,0
#set format y "%g"
#set yrange [-70:0]
#set ytics format '%g'
#set ytics ("" -120, "" -110, "" -100, "" -90, "" -80, "" -70, "-60" -60, "-50" -50, "-40" -40, "-30" -30, "-20" -20, "-10" -10, "0" 0) nomirror

plot file11 using 1:2 with lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "Benchmark", \
result11 using 1:3 with line lt -1 lw 4 lc 7 title "PE"

set label "{/:Bold L1}" at screen 0.45,0.77 font "Cambria, 16"

###############
###PLOT 2
###############
set tmargin at screen 0.80; set bmargin at screen 0.55
set lmargin at screen 0.50; set rmargin at screen 0.90
set origin 0.5,0.5
set size 0.5,0.5
unset key
#unset ylabel

set xtics format ''
set ytics format '' nomirror
set ytics scale 0
set ytics ("" -120, "" -110, "" -100, "" -90, "" -80, "" -70, "" -60, "" -50, "" -40, "" -30, "" -20, "" -10, "" 0) nomirror

set xtics ("" 0, "" 500, "" 1000, "" 1500, "" 2000, "" 2500, "" 3000, "" 3500, "" 4000, "" 4500, "" 5000) nomirror

plot file12 using 1:2 w lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "Benchmark", \
result12 using 1:3 w l lt -1 lw 4 lc 7 title "PE"

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
set xtics format ''
set xtics scale 0
set xtics ("0" 0, "0.5" 500, "1" 1000, "1.5" 1500, "2" 2000, "2.5" 2500, "3" 3000, "3.5" 3500, "4" 4000, "" 4500, "" 5000) nomirror

set ylabel "TL (dB)" offset 1,0
set format y "%g"
set yrange [-70:0]
set ytics format '%g'
set ytics ("" -120, "" -110, "" -100, "" -90, "" -80, "" -70, "-60" -60, "-50" -50, "-40" -40, "-30" -30, "-20" -20, "-10" -10, "0" 0) nomirror

plot file21 using 1:2 with lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "Benchmark", \
result21 using 1:3 w l lt -1 lw 4 lc 7 title "PE"

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

plot file22 using 1:2 with lp lt -1 lw 4 pi -20 pt 7 ps 0.2 title "Benchmark", \
result22 using 1:3 with l lt -1 lw 4 lc 7 title "PE"

unset multiplot
