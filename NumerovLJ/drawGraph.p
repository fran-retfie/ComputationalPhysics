
set datafile separator ','

set title "Lennard-Jones potential bounded states"
set xlabel "r"
set ylabel "{/Symbol y}(r)"

filename(n,l) = sprintf("dati/plot%02i%02i",n,l)
plot for [n=0:2] for [l=0:2] filename(n,l) using 2:3 with lines title (sprintf("n = %i, l = %i", n,l))

pause -1
