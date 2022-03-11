
set datafile separator ','

set terminal pngcairo size 700, 550 enhanced
set output 'plots/HO1DerrorN.png'

set title "Error in determination of Energy varing mesh size"
set xlabel "N"
set ylabel "Error"
set grid xtics ytics mytics
set mytics 5

set format y '%g'
set logscale y

set xrange [8:400000]
set format x '%g'
set logscale x

plot "dati/energiesN.csv" using 2:3 with points pointtype 2 notitle

pause -1
