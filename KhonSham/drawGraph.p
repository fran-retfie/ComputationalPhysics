
set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

set title "Independent electron model, Ground state"
set xlabel "r"
set ylabel "{/Symbol y}(r)"

filename(n) = sprintf("dati/plot%02i",n)
plot for [n=1:4] filename(n) using 2:3 with lines title sprintf("n = %d",n)
#plot for [n=1:4] filename(n) using 2:4 with lines title sprintf("n = %d",n)

pause -1
