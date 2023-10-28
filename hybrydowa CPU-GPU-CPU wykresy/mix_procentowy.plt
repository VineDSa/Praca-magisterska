reset

set encoding utf8
set term epscairo enh color dl 3 font "Helvetica, 10" size 4.5, 4.3

set output 'mix_2_procs_procent.eps'

set nokey 
set xlabel 'rozmiar macierzy - N'
set ylabel 'procent pełnego czasu(%)'
set label 'rozwiązywanie' at 10000,85 front
set label 'wypełnianie' at 6000,45 front
set yrange [0:100]
set xrange [512:14336]
set style data lines
plot 'mix_2_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'mix_2_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"

set output 'mix_3_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 9000,90 front
set label 'wypełnianie' at 6000,50 front
set xrange [512:12288]
plot 'mix_3_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'mix_3_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"
     
set output 'mix_4_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 7000,90 front
set label 'wypełnianie' at 5000,40 front
set xrange [512:10240]
plot 'mix_4_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'mix_4_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"     
     
set output 'mix_6_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 6000,90 front
set label 'wypełnianie' at 4000,40 front
set xrange [512:8192]
plot 'mix_6_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'mix_6_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"

set output 'mix_12_procs_procent.eps'
unset label
set label 'rozwiązywanie' at 4000,90 front
set label 'wypełnianie' at 2500,40 front
set xrange [512:5120]
plot 'mix_12_procs.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'mix_12_procs.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE" 
     
     
