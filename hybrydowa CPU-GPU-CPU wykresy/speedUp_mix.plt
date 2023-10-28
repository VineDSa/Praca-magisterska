reset

set encoding utf8
set term epscairo enh dl 3 font "Helvetica, 11" size 4.5, 4.3

set output 'SpeedUp_mix_2_procs.eps'

set grid
set xlabel 'rozmiar macierzy - N'
set ylabel 'przyśpieszenie względem CPU'
unset yrange
set xrange [512:14848]
set yrange [0:4]
set xtics 512,2048,14848
set key left width -2 height 0.2 box

plot 'speedUp_mix_2_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_mix_2_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_mix_2_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_mix_2_procs.eps'
set yrange [0:5]
plot 'speedUp_mix_2_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości' w l lt -1 dt 1 lw 3
     
set output 'SpeedUp_mix_3_procs.eps'
set key left width -2 height 0.2 box
unset yrange  
set xrange [512:12800]
set yrange [0:3]
set xtics 512,1024,12800    
plot 'speedUp_mix_3_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_mix_3_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_mix_3_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_mix_3_procs.eps'
set key left width -2 height 0.2 box
set yrange [0:7]
plot 'speedUp_mix_3_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości' w l lt -1 dt 1 lw 3
  
set output 'SpeedUp_mix_4_procs.eps' 
unset yrange  
set xrange [512:10752]
set yrange [0:4]
set xtics 512,1024,10752  
plot 'speedUp_mix_4_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_mix_4_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_mix_4_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_mix_4_procs.eps'
set yrange [0:12]
plot 'speedUp_mix_4_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 8 częstotliwości' w l lt -1 dt 1 lw 3
     
set output 'SpeedUp_mix_6_procs.eps'  
unset yrange  
set xrange [512:8704]
set yrange [0:3]
set xtics 512,1024,8704
plot 'speedUp_mix_6_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_mix_6_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_mix_6_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_mix_6_procs.eps'
set yrange [0:13]
plot 'speedUp_mix_6_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości' w l lt -1 dt 1 lw 3

set output 'SpeedUp_mix_12_procs.eps' 
unset yrange   
set xrange [512:5120]
set yrange [0:3]
set xtics 512,512,5120
plot 'speedUp_mix_12_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_mix_12_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_mix_12_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_mix_12_procs.eps'
set yrange [0:14]
plot 'speedUp_mix_12_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 12 częstotliwości' w l lt -1 dt 1 lw 3      

set output 'gain_mix_all_procs.eps'    
set yrange [0:15]
set xrange [512:14848]
set xtics 512,2048,14848
plot 'speedUp_mix_2_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości - 2 procs' w l lt -1 dt 1 lw 3,\
     'speedUp_mix_3_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości - 3 procs' w l lt 1 dt 1 lw 3,\
     'speedUp_mix_4_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 8 częstotliwości - 4 procs' w l lt 2 dt 1 lw 3,\
     'speedUp_mix_6_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości - 6 procs' w l lt 3 dt 1 lw 3,\
     'speedUp_mix_12_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 12 częstotliwości - 12 procs' w l lt 4 dt 1 lw 3       
