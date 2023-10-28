reset

set encoding utf8
set term epscairo enh dl 3 font "Helvetica, 11" size 4.5, 4.3

set grid
set xlabel 'rozmiar macierzy - N'
set ylabel 'przyśpieszenie względem CPU'
unset yrange
set xrange [512:14848]
set xtics 512,2048,14848
set key center top width -2 height 0.2 box
set output 'SpeedUp_gpu_2_procs.eps'
plot 'speedUp_gpu_2_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_gpu_2_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_gpu_2_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_gpu_2_procs.eps'
set key left top width -2 height 0.2 box
set yrange [0.4:6]
plot 'speedUp_gpu_2_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości' w l lt -1 dt 1 lw 3
     
set output 'SpeedUp_gpu_3_procs.eps'
unset yrange  
set xrange [512:12800]
set xtics 512,1024,12800    
plot 'speedUp_gpu_3_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_gpu_3_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_gpu_3_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_gpu_3_procs.eps'
set yrange [0.4:6]
plot 'speedUp_gpu_3_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości' w l lt -1 dt 1 lw 3
set output 'SpeedUp_gpu_4_procs.eps' 
set key center top width -2 height 0.2 box
unset yrange  
set xrange [512:10752]
set xtics 512,1024,10752  
plot 'speedUp_gpu_4_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_gpu_4_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_gpu_4_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_gpu_4_procs.eps'
set key left top width -2 height 0.2 box
set yrange [0.4:6]
plot 'speedUp_gpu_4_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 8 częstotliwości' w l lt -1 dt 1 lw 3
     
set output 'SpeedUp_gpu_6_procs.eps'  
unset yrange  
set xrange [512:8704]
set xtics 512,1024,8704
plot 'speedUp_gpu_6_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_gpu_6_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_gpu_6_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_gpu_6_procs.eps'
set yrange [0.4:6]
plot 'speedUp_gpu_6_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości' w l lt -1 dt 1 lw 3

set output 'SpeedUp_gpu_12_procs.eps' 
unset yrange   
set xrange [512:5120]
set xtics 512,512,5120
plot 'speedUp_gpu_12_procs.txt' using 1:2 smooth bezier t 'Z matrix - Wypełnianie' w l lt -1 dt 4 lw 3, \
     'speedUp_gpu_12_procs.txt' using 1:3 smooth bezier t 'Rozwiązywanie' w l lt -1 dt 2 lw 3, \
     'speedUp_gpu_12_procs.txt' using 1:4 smooth bezier t 'Suma' w l lt -1 dt 1 lw 3
     
set output 'gain_gpu_12_procs.eps'
set yrange [0.4:6]
plot 'speedUp_gpu_12_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 12 częstotliwości' w l lt -1 dt 1 lw 3    

set output 'gain_gpu_all_procs.eps'    
set yrange [0.4:6]
set xrange [512:14848]
set xtics 512,2048,14848
plot 'speedUp_gpu_2_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości - 2 procs' w l lt -1 dt 1 lw 3,\
     'speedUp_gpu_3_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości - 3 procs' w l lt 1 dt 1 lw 3,\
     'speedUp_gpu_4_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 8 częstotliwości - 4 procs' w l lt 2 dt 1 lw 3,\
     'speedUp_gpu_6_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 6 częstotliwości - 6 procs' w l lt 3 dt 1 lw 3,\
     'speedUp_gpu_12_procs.txt' using 1:5 smooth bezier t 'Zysk czasu dla 12 częstotliwości - 12 procs' w l lt 4 dt 1 lw 3        
