
set datafile separator ','

filename(l,f) = sprintf("dati/plot%02i_%02i",l,f)
plot for [l=0:3] for [f=0:6] filename(l,f) using 2:3 with lines title sprintf("l = %i, f = %i",l,f)

pause -1
