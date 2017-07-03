set terminal x11 enhanced font "Times,12"
file="CNPE2_FIELD.dat"
set size 1,1
set origin 0.0,0.0
set title "Transmission Loss LP(x,z)"
set palette defined ( 0 "#000090",\
                      1 "#000fff",\
                      2 "#0090ff",\
                      3 "#0fffee",\
                      4 "#90ff70",\
                      5 "#ffee00",\
                      6 "#ff7000",\
                      7 "#ee0000",\
                      8 "#7f0000")
set view map
set xlabel "r (m)" offset 0,0
set ylabel "z (m)" offset -3,0
set xrange [0:10000]
set yrange [0:5000]
set zrange [-120:0]
set format x "%g"
set format y "%g"
set xtics 1000 nomirror
set ytics 200 nomirror
set tics font "Cambria,12"
set pm3d map
set samples 100
set isosamples 100
splot file using 1:2:3 title "Transmission Loss" with pm3d
set pm3d interpolate 0,0


#set terminal png enhanced font "Cambria, 12"
#set output "topo2.png"
#set rmargin 5
#file="LP_L1_f1.dat"
#set size 1,1
#set origin 0.0,0.0
#set title "Transmission Loss LP(x,z)"
#set palette color positive
#set view map
#set xlabel "r (km)" offset 0,0.5
#set ylabel "z (m)" offset -2,0
#set xrange [0:5000]
#set yrange [100:500]
#set zrange [-120:20]
#set format x "%g"
#set format y "%g"
#set xtics ("0" 0, "1" 1000, "2" 2000, "3" 3000, "4" 4000, "5" 5000) nomirror
#set ytics 100 nomirror
#set tics font "Cambria,12"
#set pm3d map
#set samples 100
#set isosamples 100
#splot file using 1:2:3 with pm3d
#set pm3d interpolate 0,0
