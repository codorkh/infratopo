set terminal png enhanced font "Cambria:Bold, 16" size 1024, 1024
set output "LP_SnT_2017.png"

file="temp.dat"
file1="LP_L3_f4.dat"

set multiplot
unset key
set lmargin at screen 0.10
set rmargin at screen 0.50
set tmargin at screen 0.60
set bmargin at screen 0.20
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
set xlabel "r (km)" offset 0,0.5
set ylabel "z (m)" offset -2,0
set xrange [0:4500]
set yrange [100:800]
set zrange [-100:20]
set format x "%g"
set format y "%g"
set xtics ("0" 0, "1" 1000, "2" 2000, "3" 3000, "4" 4000, "5" 5000) nomirror
set ytics 100 nomirror
set tics font "Cambria,12"
set pm3d map
set samples 100
set isosamples 100
splot file using 1:2:3 with pm3d
set pm3d interpolate 10,10

set key
set lmargin at screen 0.50
set rmargin at screen 0.90
set tmargin at screen 0.60
set bmargin at screen 0.20
unset ylabel
unset ytics
splot file1 using 1:2:3 with pm3d
set pm3d interpolate 10,10
unset multiplot
