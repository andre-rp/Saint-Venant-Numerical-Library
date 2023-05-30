#set terminal qt 0
#set xlabel "x (m)" font "Helvetica,14"
#set ylabel "zb+h (m)" font "Helvetica,14"
#set xrange [-30:30]
#set yrange [0:3.5]
#set title "t = 2s"
#set tics font "Helvetica,10"
#plot "db_ERiemann.out" using 1:2 with line lt rgb "blue" title "Método de Hartree"

set terminal qt 1
set lmargin 8
set tics font "Helvetica,12"
set multiplot layout 3,1
unset xlabel
set ylabel "h+zb (m)" font "Helvetica,14"
set xrange [-30:30]
set yrange [0:3.5]
plot "db_ERiemann.out" using 1:2 with line lt rgb "blue" title "S Exata"
unset title
set ylabel "q (m^2/s)" font "Helvetica,14"
set xrange [-30:30]
set yrange [0:5]
plot "db_ERiemann.out" using 1:3 with line lt rgb "blue" title ""
set xlabel "x (m)" font "Helvetica,14"
set ylabel "v (m/s)" font "Helvetica,14"
set xrange [-30:30]
set yrange [0:11]
plot "db_ERiemann.out" using 1:4 with line lt rgb "blue" title ""
