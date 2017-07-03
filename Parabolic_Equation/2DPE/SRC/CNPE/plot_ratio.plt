set terminal png enhanced font "Cambria:Bold, 11"
set output "Amp_ratio.png"

set multiplot layout 1,3 scale 1,0.5
file2="SnT_AWE_L2.dat"
file3="SnT_AWE_L3.dat"
file4="SnT_AWE_L4.dat"
result2="SnT_P_L2.dat"
result3="SnT_P_L3.dat"
result4="SnT_P_L4.dat"

#set label "Amplitute ratios at microbarographs L1-4" font "Times,14" at screen 0.50,0.77 center

###############

set tmargin at screen 0.60
set bmargin at screen 0.20
set lmargin at screen 0.15
set rmargin at screen 0.40
set border lw 2
set border 7
unset key
#set xlabel "f (Hz)" offset 0.0,-1.0
set ylabel "P/P1" offset 0.8,0.0
set xrange [0:9.59]
set yrange [0:2.5]
set ytics 0.5 nomirror
set xtics ("1.9" 1.90, "3.79" 3.79, "5.69" 5.69, "7.59" 7.59) rotate by 45 offset -0.8,-1.5
set title "L2/L1" offset 0.0,-3
set bars 3.0
plot 0.89 lw 2 lc 1 t "flat", file2 u 1:2:3 lc 4 lw 2 pt 1 w e t "Data", result2 u 1:3 lt -1 pt 5 lc 7 t "PE"

###############

set border 5
set arrow 5 from 0,0 to 0,2.5 nohead
set key out horiz
set key width -20 spacing 1
set key at screen 1.10,0.66
unset ylabel
unset ytics
set xlabel "f (Hz)" offset 0.0,-1.0
set format y ''
set tmargin at screen 0.60
set bmargin at screen 0.20
set lmargin at screen 0.40
set rmargin at screen 0.65
set title "L3/L1" offset 0.0,-3
unset ylabel
plot 0.63 lw 2 lc 1 t "Flat", file3 u 1:2:3 lc 4 lw 2 pt 1 w e t "Data", result3 u 1:3 lt -1 pt 5 lc 7 t "PE"
###############

set border 13
unset key
set tmargin at screen 0.60
set bmargin at screen 0.20
set lmargin at screen 0.65
set rmargin at screen 0.90
unset xlabel
set title "L4/L1" offset 0.0,-3
plot 1.58 lw 2 lc 1 t "flat", file4 u 1:2:3 lc 4 lw 2 pt 1 w e t "Data", result4 u 1:3 lt -1 pt 5 lc 7 t "PE"
