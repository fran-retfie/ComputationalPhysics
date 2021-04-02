
set datafile separator ','

plot "energies.csv" using 2:3 with points pointtype 5

# Draw 5 vertical lines
n = 30

# ... evenly spaced between x0 and x1
x0 = 2.5
x1 = 3.5
dx = (x1-x0)

# ... each line going from y0 to y1
y0 = 0
y1 = 10

do for [i = 0:n-1] {
    x = x0 + i*dx
    set arrow from y0,x to y1,x #nohead linecolor "blue" # add other styling options if needed
}

pause -1
