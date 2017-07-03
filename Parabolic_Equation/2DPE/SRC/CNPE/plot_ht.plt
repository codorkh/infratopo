set terminal png enhanced font "Cambria:Bold, 14" size 1024, 1024
set output "HT.png"
file1="CNPE2_FIELD.dat"
file2="PE.dat"
file3="FEM.dat"
file4="CNPE2_TL.dat"

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
set ylabel "z (km)" offset -1,0
unset xlabel
set xrange [0:6680]
set yrange [-20:3000]
set zrange [-120:0]
set format x "%g"
set format y "%g"
set xtics format ''
set ytics ("" 0, "0.8" 800, "1.6" 1600, "2.4" 2400, "3.2" 3200)
set pm3d map
set samples 100
set isosamples 100
splot file1 using 1:2:3 with pm3d
set pm3d interpolate 0,0

set rmargin at screen 0.85; set lmargin at screen 0.15
set bmargin at screen 0.30; set tmargin at screen 0.50
set grid xtics
set xlabel "r (m)" offset 0,0
set ylabel "TL (dB)" offset 0.7,0
set y2label "H (km)" offset 1,0
set xrange [0:6680]
set yrange [-120:0]
set y2range [0:2000]
set format x "%g"
set format y "%g"
set format y2 "%g"
set xtics 1000 nomirror
set ytics 20 nomirror
set y2tics ("" 0, "0.5"  500, "1.0" 1000, "1.5" 1500) nomirror
plot file2 w l lt -1 lc 7 lw 4 title "PE", \
file4 w lp lt -1 lw 3 lc -1 axes x1y2 title "Terrain", \
file3 u 1:($2-10) w l lt -1 lc 6 lw 4 title "FEM"
unset multiplot
