set key autotitle columnheader
set format y "10^{%L}"
set log y

set style line 80 lt 0 lc rgb "#808080"

set style line 1 lt 1 lc rgb "#A00000" lw 2 pt 7 ps 1.0
set style line 2 lt 1 lc rgb "#00A000" lw 2 pt 11 ps 1.0
set style line 3 lt 1 lc rgb "#5060D0" lw 2 pt 9 ps 1.0
set style line 4 lt 1 lc rgb "#0000A0" lw 2 pt 8 ps 1.0
set style line 5 lt 1 lc rgb "#D0D0FA" lw 2 pt 13 ps 1.0
set style line 6 lt 1 lc rgb "#0066FE" lw 2 pt 12 ps 1.0
set style line 7 lt 1 lc rgb "#B200B2" lw 2 pt 5 ps 1.0
set style line 8 lt 1 lc rgb "#C200B2" lw 2 pt 6 ps 1.0
set style line 9 lt 1 lc rgb "#FF0000" lw 2 pt 4 ps 1.0

set grid xtics
set grid ytics

# nomirror means do not put tics on the opposite side of the plot
set xtics nomirror
set ytics nomirror

set title 'Species number in simulation cell'
set xlabel 'time (days)'

f='output/1.dat'
plot for [i=3:9] f  u ($1*1.15741e-5):i w li ls (i-2), for[i=10:11] f u ($1*1.15741e-5):i w lp lt i
