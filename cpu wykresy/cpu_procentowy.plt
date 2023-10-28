reset

set encoding utf8
set term epscairo enh color dl 3 font "Helvetica, 10" size 4.5, 4.3

set output 'CPU_procent.eps'

set nokey 
set xlabel 'rozmiar macierzy - N'
set ylabel 'procent pełnego czasu(%)'
set label 'rozwiązywanie' at 13000,85 front
set label 'wypełnianie' at 6000,40 front
set yrange [0:100]
set xrange [512:20480]
set style data lines
plot 'cpu.txt' u 1:($5*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#F5F5DC", \
     'cpu.txt' u 1:(($2+$3)*100/$4):($2*100/$4) w filledcu lt 1 lc rgb "#B0C4DE"

