set terminal x11 font "Cambria,12" enhanced
#set terminal latex
#set output "test.tex"
file="./Test_LP.dat"
file1="./Test_LPg.dat"
file2="./Test_Cond.dat"
file3="./Test_Terr.dat"
set multiplot layout 2,1
set tmargin 3
#########
set title "Sound Pressure Level LP(x,0)"
set xlabel "r (m)"
set ylabel "LP (dB)"
set y2label "H (m)"
set xrange [0:5000]
set yrange [-80:10]
set y2range [0:2000]
set format x "%g"
set format y "%g"
set xtics 500
set ytics 10 nomirror
set y2tics 400 nomirror
set tics font "Cambria,12"
plot file1 using 1:2 with lines title "LP(x,0)", \
file1 using 1:3 with lines title "LPrms(x,0)", \
file3 using 1:2 axes x1y2 with lines title "Terrain"
#########
set title "Sound Pressure Level field LP(x,z)"
set palette color positive
set view map
set cblabel "LP (dB)"
set xlabel "r (m)"
set ylabel "z (m)"
set xrange [0:5000]
set yrange [0:3000]
set format x "%g"
set format y "%g"
set xtics 0,500,5000 nomirror
set ytics 0,500,3000 nomirror
set tics font "Cambria,12"
set pm3d interpolate 0,0
splot file
#########
unset multiplot
