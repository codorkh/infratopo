set xlabel "r (m)"
set ylabel "LP (dB)"
set format x "%g"
set format y "%g"
set xtic 0,500,5000
m="./PE2D_HFR_LPg.dat"
set terminal x11 font "Helvetica,16"
set tics font "Helvetica,16"
set output "test.tex"
set nokey
set grid
set title "LP at ground level (z=0)"
plot m using 1:2 with lines
