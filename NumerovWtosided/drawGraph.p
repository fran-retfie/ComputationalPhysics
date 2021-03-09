
set datafile separator ','

plot "dati.csv" using 2:3 with lines, "dati.csv" using 2:4 with lines, "dati.csv" using 2:5 with lines

pause -1
