
set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

set title "Independent electron model"
set xlabel "r"
set ylabel "{/Symbol y}(r)"

filename(n,l) = sprintf("dati/plot%02i%02i%02i", 6, n, l)
#plot for [n=1:2] for [l=0:3] filename(n,l) using 2:3 with lines title (sprintf("n = %i, l = %i", n, l))
plot for [n=1:2] for [l=0:3] filename(n,l) using 2:4 with lines title (sprintf("n = %i, l = %i", n, l))

pause -1