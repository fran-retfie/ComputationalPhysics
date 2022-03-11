
set datafile separator ','

set xrange [-0.5:2.5]
set yrange [1:8]
set title "3D Harmonic oscillator states"
set ylabel "E"
set grid ytics mytics
set mytics 2
unset xtics

plot "energies.csv" using 2:3:(sprintf("k=%d, l=%d", $1, $2)) with labels point pt 7 ps 1.5 lc rgb "blue" offset char 1,1 notitle

pause -1
