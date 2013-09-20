reset
dx=5.
n=2
total_box_width_relative=0.70
gap_width_relative=0.3
d_width=(gap_width_relative+total_box_width_relative)*dx/2.
reset

#set terminal postscript eps
#set output 'Mixed_DataRate.eps'

set terminal png
set output 'Mixed_DataRate.png'

set xlabel "Number of Users"
set ylabel "Aggreagated Throughput for Seed 1 [Mbps]"
set grid
set boxwidth total_box_width_relative/n relative

set xtics ("1" 2.25,"8" 9.25,"15" 16.25,"22" 23.25,"29" 30.25 ,"36" 37.25,"43" 44.25,"50" 51.25,"57" 58.25)
set yrange [0:26]

set style fill transparent solid 0.8 noborder
plot "Mixed_dataRate.dat" u 1:2 w boxes lc rgb"black" title "Simulation result", "Mixed_dataRate.dat" u ($1+d_width):3 w boxes lc rgb"gray" title "Analytical result"
