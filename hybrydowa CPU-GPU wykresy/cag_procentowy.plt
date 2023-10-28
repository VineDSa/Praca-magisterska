reset

set encoding utf8
set term epscairo enh color dl 3 font "Helvetica, 10" size 4.5, 4.3

set output 'cag_2_procs_procent.eps'

set nokey 
set xlabel 'rozmiar macierzy - N'
set ylabel 'procent pełnego czasu(%)'
set label 'rozwiązywanie' at 7000,60 front
set label 'wypełnianie' at 1500,15 front
set yrange [0:100]
set xrange [512:14336]
set style data lines
plot 'cag_2_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'cag_2_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"

set output 'cag_3_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 6000,60 front
set label 'wypełnianie' at 1500,15 front
set xrange [512:12288]
plot 'cag_3_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'cag_3_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"    
     
set output 'cag_4_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 5000,60 front
set label 'wypełnianie' at 1500,15 front
set xrange [512:10240]
plot 'cag_4_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'cag_4_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"     
     
set output 'cag_6_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 4000,60 front
set label 'wypełnianie' at 1500,15 front
set xrange [512:8192]
plot 'cag_6_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'cag_6_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"

set output 'cag_12_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 2500,60 front
set label 'wypełnianie' at 1500,15 front
set xrange [512:5120]
plot 'cag_12_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'cag_12_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE" 
     
     
