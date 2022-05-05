
#set datafile separator ','

#set terminal pngcairo size 700, 550 enhanced
#set output 'plots/HO1D.png'

#set title "Electron density"
#set xlabel "r"
#set ylabel "{/Symbol r}(r)"
set key autotitle columnheader


set grid

data = "dati/E1.12.csv"

FILES = system("ls -1 dati/*.csv")
#plot for [data in FILES] data using 1:2 title data w line
plot data using 1:6 w line title "kin", data using 1:7 w line title "kin2", data using 1:8 w line title "pot", data using 1:9 w line title "pot2" ,\
     data using 1:($6 + $8) w line title "E", data using 1:($7 + $8) w line title "E2"

pause -1
