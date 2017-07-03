set terminal png enhanced font "Cambria:Bold, 14" size 1024, 1024
set output "Alps.png"
file1="LP_ALPS.dat"
file2="LPg_ALPS.dat"

set multiplot
set rmargin at screen 0.85; set lmargin at screen 0.15
set bmargin at screen 0.50; set tmargin at screen 0.80
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
set ylabel "z (km)" offset -0.7,0
unset xlabel
set xrange [0:550000]
set yrange [-500:24000]
set ytics ("" 0, "5.0" 5000, "10.0" 10000, "15.0" 15000, "20.0" 20000)
set zrange [-120:0]
set format x "%g"
set format y "%g"
set xtics format ''
set ytics  nomirror
set pm3d map
set samples 100
set isosamples 100
splot file1 using 1:2:3 with pm3d
set pm3d interpolate 0,0

set rmargin at screen 0.85; set lmargin at screen 0.15
set bmargin at screen 0.30; set tmargin at screen 0.50
set grid xtics
set xlabel "r (km)" offset 0,0
set ylabel "TL (dB)" offset 0.1,0
set y2label "H (km)" offset 1,0
set xrange [0:550000]
set yrange [-120:0]
set y2range [0:6000]
set format x "%g"
set format y "%g"
set format y2 "%g"
set xtics ("0" 0, "50" 50000, "100" 100000, "150" 150000, "200" 200000, "250" 250000, "300" 300000, "350" 350000, "400" 400000, "450" 450000, "500" 500000, "550" 550000) nomirror
set ytics 20 nomirror
set y2tics ("" 0, "2"  2000, "4" 4000) nomirror
plot file2 using 1:3 w l lt -1 lc 7 lw 4 title "PE", \
file2 using 1:2 w lp lt -1 lw 3 lc -1 axes x1y2 title "Terrain"

unset multiplot
