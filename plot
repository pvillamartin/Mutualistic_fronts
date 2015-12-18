
set style line 1 lc rgb '#D53E4F' lt 1 lw 10 # red
set style line 2 lc rgb '#F46D43' lt 1 lw 10 # orange
set style line 3 lc rgb '#FDAE61' lt 1 lw 10 # pale orange
set style line 4 lc rgb '#FEE08B' lt 1 lw 10 # pale yellow-orange
set style line 5 lc rgb '#E6F598' lt 1 lw 10 # pale yellow-green
set style line 6 lc rgb '#ABDDA4' lt 1 lw 10 # pale green
set style line 7 lc rgb '#66C2A5' lt 1 lw 10 # green
set style line 8 lc rgb '#3288BD' lt 1 lw 10 # blue

set term postscript enhanced color
set output 'Eq_Pairs.eps'
set title 'Equally distributed initially'
set xlabel 'Time'
set ylabel 'Pair in collisions'
plot 'Collision_detailed_100000_3_1_5000_5.txt' u 1:($3/$11) w l ls 1 title 'AB+BA',\
'Collision_detailed_100000_3_1_5000_5.txt' u 1:($4/$11) w l ls 2 title 'AC+CA',\
'Collision_detailed_100000_3_1_5000_5.txt' u 1:($9/$11) w l ls 3 title 'BC+CB'

set output 'Uneq_Pairs.eps'
set title 'Unequally distributed initially'
set xlabel 'Time'
set ylabel 'Pair in collisions'
plot 'Collision_detailed_100000_3_2_5000_5.txt' u 1:($3/$11) w l ls 1 title 'AB+BA',\
'Collision_detailed_100000_3_2_5000_5.txt' u 1:($4/$11) w l ls 2 title 'AC+CA',\
'Collision_detailed_100000_3_2_5000_5.txt' u 1:($9/$11) w l ls 3 title 'BC+CB'


set style line 1 lc rgb '#D53E4F' lt 1 lw 10 # red
set style line 2 lc rgb '#F46D43' lt 1 lw 10 # orange
set style line 3 lc rgb '#FDAE61' lt 1 lw 10 # pale orange
set style line 4 lc rgb '#FEE08B' lt 1 lw 10 # pale yellow-orange
set style line 5 lc rgb '#E6F598' lt 1 lw 10 # pale yellow-green
set style line 6 lc rgb '#ABDDA4' lt 1 lw 10 # pale green
set style line 7 lc rgb '#66C2A5' lt 1 lw 10 # green
set style line 8 lc rgb '#3288BD' lt 1 lw 10 # blue

set log
set term postscript enhanced color
set output 'Eq_Pairs_total.eps'
set title 'Equally distributed initially'
set xlabel 'Time'
set ylabel 'Pair in collisions'
plot 'Collision_detailed_total_pairs_100000_3_1_5000_6.txt' u 1:($3/$11) w l ls 1 title 'AB+BA',\
'Collision_detailed_total_pairs_100000_3_1_5000_6.txt' u 1:($4/$11) w l ls 2 title 'AC+CA',\
'Collision_detailed_total_pairs_100000_3_1_5000_6.txt' u 1:($9/$11) w l ls 3 title 'BC+CB'

set output 'Uneq_Pairs_total.eps'
set title 'Unequally distributed initially'
set xlabel 'Time'
set ylabel 'Pair in collisions'
plot 'Collision_detailed_total_pairs_100000_3_2_5000_6.txt' u 1:($3/$11) w l ls 1 title 'AB+BA',\
'Collision_detailed_total_pairs_100000_3_2_5000_6.txt' u 1:($4/$11) w l ls 2 title 'AC+CA',\
'Collision_detailed_total_pairs_100000_3_2_5000_6.txt' u 1:($9/$11) w l ls 3 title 'BC+CB'





set term png
set output 'Attempt1_L1000_V48_S3_M10_D01_T100.png'
set xlabel 'Sites'
set ylabel 'Time'
plot 'Mutual_front_1000_4_3_10_0.1_100_1.txt' u 2:1:3 lc palette

set term png
set output 'Attempt1_L1000_V48_S3_M0_D01_T100.png'
set xlabel 'Sites'
set ylabel 'Time'
plot 'Mutual_front_1000_4_3_0_0.1_100_1.txt' u 2:1:3 lc palette

set term png
set output 'Attempt1_L1000_V48_S3_M0_D0_T100.png'
set xlabel 'Sites'
set ylabel 'Time'
plot 'Mutual_front_1000_4_3_0_0_100_1.txt' u 2:1:3 lc palette


set term png
set output 'Attempt1_L1000_V8_S3_M10_D01_T100.png'
set xlabel 'Sites'
set ylabel 'Time'
plot 'Mutual_front_1000_1_3_10_0.1_100_1.txt' u 2:1:3 lc palette

set term png
set output 'Attempt1_L1000_V8_S3_M10_D0_T100.png'
set xlabel 'Sites'
set ylabel 'Time'
plot 'Mutual_front_1000_1_3_10_0_100_1.txt' u 2:1:3 lc palette

set term png
set output 'Attempt1_L1000_V48_S3_M10_D0_T100.png'
set xlabel 'Sites'
set ylabel 'Time'
plot 'Mutual_front_1000_4_3_10_0_100_1.txt' u 2:1:3 lc palette









