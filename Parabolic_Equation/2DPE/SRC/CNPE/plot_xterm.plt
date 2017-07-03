set terminal x11 enhanced font "Times,12"
file1="CNPE2_TL.dat"
file2="CNPE2_FIELD.dat"
set size 1,1
set grid
set title "Transmission Loss along (z = 0)"
set xlabel "r (m)" offset 0,0
set ylabel "LP (dB)" offset -3,0
set y2label "H (m)" offset 3,0
set xrange [0:10000]
set yrange [-120:0]
set y2range [0:1000]
set format x "%g"
set format y "%g"
set format y2 "%g"
set xtics 1000 nomirror
set ytics 20 nomirror
set y2tics 200 nomirror
plot file1 using 1:3 with lines lt -1 lc 1 lw 2 title "PE", \
file1 using 1:2 with lp lt -1 lw 3 pi -30 pt 7 ps 0.2 axes x1y2 title "Terrain"
#file3 using 1:3 axes x1y2 with lines title "Slope", \
#file3 using 1:4 axes x1y2 with lines title "Curvature"
