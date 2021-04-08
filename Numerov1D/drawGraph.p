
set datafile separator ','

set title "1D Harmonic oscillator states"
set xlabel "r"
set ylabel "{/Symbol y}(r)"

filename(n) = sprintf("dati/plot%02i",n)
plot for [n=0:4] filename(n) using 2:3 with lines

pause -1
