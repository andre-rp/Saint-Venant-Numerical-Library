set lmargin 8
set tics font "Helvetica,12"

set terminal qt 1
set xlabel "x (m)" font "Helvetica,14"
set ylabel "h+zb (m)" font "Helvetica,14"
set xrange [-30:30]
set yrange [0:3.5]
plot "db_Lax.out" using 1:2 with line lt rgb "blue" title "M Lax D"

set terminal qt 2
set xlabel "x (m)" font "Helvetica,14"
set ylabel "q (m^2/s)" font "Helvetica,14"
set xrange [-30:30]
unset yrange
plot "db_Lax.out" using 1:3 with line lt rgb "blue" title ""

set terminal qt 3
set xlabel "x (m)" font "Helvetica,14"
set ylabel "v (m/s)" font "Helvetica,14"
set xrange [-30:30]
plot "db_Lax.out" using 1:4 with line lt rgb "blue" title ""
