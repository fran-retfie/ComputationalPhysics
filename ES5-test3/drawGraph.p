
#set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

#set title "Electron density"
#set xlabel "r"
#set ylabel "{/Symbol r}(r)"
set key autotitle columnheader


set grid

FILES = system("ls -1 dati/*.csv")
#plot for [data in FILES] data using 1:2 title data w line
plot "dati/E1.01.csv" using 1:2 w line, "dati/E1.1.csv" using 1:3 w line

pause -1
