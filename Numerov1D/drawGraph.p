
set datafile separator ','

set terminal pngcairo size 700, 550 enhanced
set output 'plots/HO1D.png'

set title "1D Harmonic oscillator states"
set xlabel "r"
set ylabel "{/Symbol y}(r)"

filename(n) = sprintf("dati/plot%02i",n)
plot for [n=0:4] filename(n) using 2:3 with lines title sprintf("n = %d",n)

pause -1
