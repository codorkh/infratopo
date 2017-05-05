#set terminal x11 font "Cambria,12" enhanced
set terminal png enhanced font "Cambria, 12"
set output "test2.png"
set rmargin 5
file="../results/Flat_LP.dat"
file1="../results/Flat_LPg.dat"
file2="../results/Flat_Cond.dat"
file3="../results/Flat.dat"
set size 1,1
set origin 0.0,0.0
set title "Transmission Loss LP(x,z)"
set palette color positive
set view map
#set cblabel "LP (dB)" offset 2,0
set xlabel "r (km)" offset 0,1
set ylabel "z (km)" offset -2,0
set xrange [0:5000]
set yrange [0:3000]
set zrange [-120:20]
set format x "%g"
set format y "%g"
set xtics ("0" 0, "1" 1000, "2" 2000, "3" 3000, "4" 4000, "5" 5000) nomirror
set ytics ("0" 0, "1" 1000, "2" 2000, "3" 3000) nomirror
set tics font "Cambria,12"
set pm3d map
set samples 100
set isosamples 100
splot file using 1:2:3 with pm3d 
set pm3d interpolate 0,0
