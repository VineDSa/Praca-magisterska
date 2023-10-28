reset

set encoding utf8
set term epscairo enh dl 3 font "Helvetica, 10" size 4.5, 4.3

set output 'gpu_2_procs.eps'

set grid
set xlabel 'rozmiar macierzy - N'
set ylabel 'czas (s)'
set xrange [512:14848]
set xtics 512,2048,14848
set key left width -2 height 0.2 box

plot 'gpu_2_procs.txt' u 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'gpu_2_procs.txt' u 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'gpu_2_procs.txt' u 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3

