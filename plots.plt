reset
set terminal pngcairo enhanced font "Arial, 16" size 1280, 960
set encoding utf8
set output "firststates.png"
set size 1,1
set grid
set key outside;
set key right top;
set title "First states of the 1-D system"
set xlabel "Position"
max(a,b) = ( a > b ) ? a : b
m = "./plottingdata.dat"
stats m u 1:9 nooutput
set xrange [STATS_min_x:STATS_max_x]
set yrange [*:STATS_max_y*1.3]
plot \
m us 1:(max($3,$2)) notitle w lines, m us 1:($4*3.0+$3) title "{/Symbol Y}_0" w lines lc "dark-blue",\
m us 1:(max($5,$2)) notitle w lines, m us 1:($6*3.0+$5) title "{/Symbol Y}_1" w lines lc "dark-green",\
m us 1:(max($7,$2)) notitle w lines, m us 1:($8*3.0+$7) title "{/Symbol Y}_2" w lines lc "dark-red",\
m us 1:(max($9,$2)) notitle w lines, m us 1:($10*3.0+$9) title "{/Symbol Y}_3" w lines lc "dark-yellow",\
m us 1:3 notitle w lines lc "blue" dt 2, m us 1:5 notitle w lines lc "green" dt 2, m us 1:7 notitle  w lines lc "red" dt 2, m us 1:9 notitle w lines lc "yellow" dt 2,\
m us 1:2  lc "black" w lines title "V(x)",\
