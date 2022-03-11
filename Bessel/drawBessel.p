
#set datafile separator ','
set terminal pngcairo size 1000, 800 enhanced dashed font "Helvetica,15"
set output 'plots/bessel.png'

set title "Bessel functions"
#set xlabel "x"
#set ylabel "j(x)"
set xrange [0:15]
set yrange [-0.4:1.2]
set grid

set label 1 'j_0' at 1.9, 0.6 font "Helvetica,19"
set label 2 'j_1' at 2.8, 0.45 font "Helvetica,19"
set label 3 'j_2' at 3.7, 0.35 font "Helvetica,19"
set label 4 'j_3' at 4.9, 0.29 font "Helvetica,19"
set label 5 'j_4' at 6.2, 0.24 font "Helvetica,19"

plot "bessel.csv" using 1:2 with lines linewidth 4 lc "orange" title "recursion formula" ,\
     for [l=3:6] "bessel.csv" using 1:l with lines linewidth 4 lc "orange" notitle ,\
     "besselteo.csv" using 1:2 with lines dashtype 2 lc "black" title "theoretical formula" ,\
     for [l=3:6] "besselteo.csv" using 1:l with lines dashtype 2 lc "black" notitle

pause -1
