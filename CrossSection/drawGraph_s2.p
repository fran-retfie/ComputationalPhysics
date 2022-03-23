
set datafile separator ','

filename(l,E) = sprintf("dati/plot%02i_%02i",l,E)
l=4;
plot for [E=0:10] filename(l,E) using 2:3 with lines title sprintf("l = %i",l)

pause -1
