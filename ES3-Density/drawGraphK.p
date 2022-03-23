
set datafile separator ','

set terminal pngcairo size 1000, 1000 font 'Verdana,12' enhanced
set output 'plots/K-wavefunctions.png'


#set xlabel "r [a.u]"
#set ylabel "{/Symbol y}(r)"
set grid

set multiplot layout 2,2 title "Independent electron model - Electron wavefuction - K cluster" font 'Verdana,17'

set title "N_e = 2"
filename(n,l) = sprintf("dati/plot%02i%02i%02i", 1, n, l)
plot for [n=1:1] for [l=0:0] filename(n,l) using 2:3 with lines title (sprintf("n = %i, l = %i", n, l)) lw 2
#plot for [n=1:4] filename(n) using 2:4 with lines title sprintf("n = %d",n)

set title "N_e = 8"
filename(n,l) = sprintf("dati/plot%02i%02i%02i", 2, n, l)
plot for [n=1:1] for [l=0:1] filename(n,l) using 2:3 with lines title (sprintf("n = %i, l = %i", n, l)) lw 2

set title "N_e = 18"
filename(n,l) = sprintf("dati/plot%02i%02i%02i", 3, n, l)
plot for [n=1:1] for [l=0:2] filename(n,l) using 2:3 with lines title (sprintf("n = %i, l = %i", n, l)) lw 2

set title "N_e = 20"
filename(n,l) = sprintf("dati/plot%02i%02i%02i", 3, n, l)
plot for [n=1:2] for [l=0:(2-n)*2] filename(n,l) using 2:3 with lines title (sprintf("n = %i, l = %i", n, l)) lw 2
