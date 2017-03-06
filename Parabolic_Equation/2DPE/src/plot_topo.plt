#set terminal x11 font "Cambria,12" enhanced
set terminal png enhanced font "Cambria, 12"
set output "topo.png"
file="../results/Topo_LP.dat"
file1="../results/Topo_LPg.dat"
file2="../results/Topo_Cond.dat"
file3="../results/Topo.dat"
set size 1,1
#set multiplot
#set size 1,0.5
#set origin 0.0,0.5
set tmargin 3
#########
set title "Transmission Loss at ground LP(x,0)"
set xlabel "r (km)" offset 0,1
set ylabel "LP (dB)"
set y2label "H (km)"
set xrange [0:5000]
#set yrange [-120:10]
set y2range [0:2000]
set format x "%g"
set format y "%g"
set xtics ("0" 0, "1" 1000, "2" 2000, "3" 3000, "4" 4000, "5" 5000) nomirror
set ytics 20 nomirror
set y2tics ("0" 0, "0.4" 400, "0.8" 800, "1.2" 1200, "1.6" 1600, "2" 2000)
set tics font "Cambria,12"
plot file1 using 1:2 with lines title "LP(x,0)", \
file1 using 1:3 with lines title "LPrms(x,0)", \
file3 using 1:2 axes x1y2 with lines title "Terrain"
#file3 using 1:3 axes x1y2 with lines title "Slope", \
#file3 using 1:4 axes x1y2 with lines title "Curvature"
#########
#set size 1,0.5
#set origin 0.0,0.0
#set title "Transmission Loss LP(x,z)"
#set palette color positive
#set view map
#set cblabel "LP (dB)" offset 2,0
#set xlabel "r (m)" offset 0,1
#set ylabel "z (m)" offset -2,0
#set xrange [0:5000]
#set yrange [0:3000]
#set zrange [-120:20]
#set format x "%g"
#set format y "%g"
#set xtics ("0" 0, "1" 1000, "2" 2000, "3" 3000, "4" 4000, "5" 5000) nomirror
#set ytics ("0" 0, "1" 1000, "2" 2000, "3" 3000) nomirror
#set tics font "Cambria,12"
#set pm3d map
#set samples 100
#set isosamples 100
#splot file using 1:2:3 with pm3d 
#set pm3d interpolate 0,0
#########
#unset multiplot
#set output
