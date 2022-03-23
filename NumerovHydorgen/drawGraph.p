
set datafile separator ','
set terminal qt enhanced
set termoption enhanced

set title "3D Harmonic oscillator states"
set xlabel "r"
set ylabel "{/Symbol y}(r)"

filename(n,l) = sprintf("dati/plot%02i%02i",n,l)
plot for [l=0:3] for [n=1:4] filename(n,l) using 2:3 with lines title (sprintf("n = %i, l = %i", n,l))

pause -1
