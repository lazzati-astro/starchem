set key autotitle columnheader
set format y "10^{%L}"
set log y

plot for [i=2:7] "<(tail -n +2 output1.dat)" u 1:i w l, for [i=8:9] "<(tail -n +2 output1.dat)" u 1:i w lp
