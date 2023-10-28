reset

set encoding utf8
set term epscairo enh color dl 3 font "Helvetica, 10" size 4.5, 4.3

set output 'gpu_2_procs_procent.eps'

set nokey 
set xlabel 'rozmiar macierzy - N'
set ylabel 'procent pełnego czasu(%)'
set label 'rozwiązywanie' at 7000,60 front
set label 'wypełnianie' at 550,2 front
set yrange [0:100]
set xrange [512:14336]
set style data lines
plot 'gpu_2_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'gpu_2_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"

set output 'gpu_3_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 6000,60 front
set label 'wypełnianie' at 550,2 front
set xrange [512:12288]
plot 'gpu_3_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'gpu_3_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"    
     
set output 'gpu_4_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 5000,60 front
set label 'wypełnianie' at 550,2 front
set xrange [512:10240]
plot 'gpu_4_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'gpu_4_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"     
     
set output 'gpu_6_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 4000,60 front
set label 'wypełnianie' at 550,2 front
set xrange [512:8192]
plot 'gpu_6_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'gpu_6_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"

set output 'gpu_12_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 2500,60 front
set label 'wypełnianie' at 550,2 front
set xrange [512:5120]
plot 'gpu_12_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'gpu_12_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE" 
     
     
