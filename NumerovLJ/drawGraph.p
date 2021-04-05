
set datafile separator ','

filename(n,l) = sprintf("dati/plot%02i%02i",n,l)
#plot for [n=0:2] for [l=0:2] filename(n,l) using 2:3 with lines title sprintf("n = %i, l = %i", n,l)

plot "dati/plot0000" using 2:3 with lines

pause -1
