
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
plot "dati/E1.15.csv" using 1:6 w line title "kin", "dati/E1.15.csv" using 1:7 w line title "kin2", "dati/E1.15.csv" using 1:8 w line title "pot" ,\
     "dati/E1.15.csv" using 1:($6 + $8) w line title "E", "dati/E1.15.csv" using 1:($7 + $8) w line title "E2"

pause -1
