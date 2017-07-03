set terminal x11
#set terminal png enhanced font "Cambria, 12"
#set output "topo.png"
file1="Gauss_LPg.dat"
file2="Gauss.dat"
set style func linespoints
set size 1,1
set tmargin 3
set grid
set title "Transmission Loss at ground - L4 -3.79Hz"
set xlabel "r (km)" offset 0,1
set ylabel "LP (dB)"
set y2label "Elevation (m)"
set autoscale x
set autoscale y
#set yrange [-120:10]
set y2range [0:1000]
set format x "%g"
set format y "%g"
set xtics ("0" 0, "1" 1000, "2" 2000, "3" 3000, "4" 4000, "5" 5000) nomirror
set ytics 10 nomirror
set y2tics ("0" 0, "50" 50, "100" 100, "150" 150, "200" 200, "250" 250, "300" 300)
set tics font "Cambria,12"
plot file1 using 1:2 with lines title "LP(x,0)" , \
file2 using 1:2 axes x1y2 with lines title "Terrain"
