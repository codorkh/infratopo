set terminal png enhanced font "Cambria, 12"
set output "topo.png"
file="../results/Topo_LP.dat"
file1="../results/Topo_LPg.dat"
file2="../results/Topo.dat"
set size 1,1
set tmargin 3
set grid
set title "Transmission Loss at ground LP(x,0)"
set xlabel "r (km)" offset 0,1
set ylabel "LP (dB)"
set y2label "H (km)"
set xrange [0:5000]
set yrange [-120:10]
set y2range [0:2000]
set format x "%g"
set format y "%g"
set xtics ("0" 0, "1" 1000, "2" 2000, "3" 3000, "4" 4000, "5" 5000) nomirror
set ytics 20 nomirror
set y2tics ("0" 0, "0.4" 400, "0.8" 800, "1.2" 1200, "1.6" 1600, "2" 2000)
set tics font "Cambria,12"
plot file1 using 1:2 with lines title "LP(x,0)", \
file1 using 1:3 with lines title "LPrms(x,0)", \
file2 using 1:2 axes x1y2 with lines title "Terrain"
#file3 using 1:3 axes x1y2 with lines title "Slope", \
#file3 using 1:4 axes x1y2 with lines title "Curvature"
