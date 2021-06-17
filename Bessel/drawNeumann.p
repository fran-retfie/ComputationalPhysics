
#set datafile separator ','
set terminal pngcairo size 1000, 800 enhanced dashed
set output 'plots/neumann.png'

set title "Neumann functions"
#set xlabel "x"
#set ylabel "j(x)"
set xrange [0:15]
set yrange [-0.4:0.4]
set grid

set label 1 'n_0' at 2.6, 0.36 font "Helvetica,19"
set label 2 'n_1' at 4.1, 0.26 font "Helvetica,19"
set label 3 'n_2' at 5.4, 0.21 font "Helvetica,19"
set label 4 'n_3' at 6.6, 0.18 font "Helvetica,19"
set label 5 'n_4' at 7.8, 0.16 font "Helvetica,19"


plot "neumann.csv" using 1:2 with lines linewidth 4 lc "0xAAAAFF" title "recursion formula" ,\
     for [l=3:6] "neumann.csv" using 1:l with lines linewidth 4 lc "0xAAAAFF" notitle ,\
     "neumannteo.csv" using 1:2 with lines dashtype 2 lc "black" title "theoretical formula" ,\
     for [l=3:6] "neumannteo.csv" using 1:l with lines dashtype 2 lc "black" notitle

pause -1
