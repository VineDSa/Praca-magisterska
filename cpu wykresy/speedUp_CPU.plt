reset
set term post enh dl 3 'Helvetica' 22
#set encoding iso_8859_2
#set encoding cp1250
set output 'TimeCPU.eps'
set grid
set xlabel 'matrix size - N'
set ylabel 'time (s)'
set xrange [256:15285]
set xtics 256,2048,15285
set key left width -2 height 0.2 box
plot 'speedUP_CPU680.txt' u 2:3 smooth bezier t 'Z matrix fill' w l lt 4 lw 3, \
     'speedUP_CPU680.txt' u 2:4 smooth bezier t 'solution' w l lt 6 lw 3, \
     'speedUP_CPU680.txt' u 2:5 smooth bezier t 'total' w l lt 1 lw 3
set output 'TimeGPU680.eps'
set key left width 1. height 1 box
plot 'speedUP_CPU680.txt' u 2:6 smooth bezier t 'Solution' w l lt 4 lw 3, \
     'speedUP_CPU680.txt' u 2:7 smooth bezier t 'Z matrix fill' w l lt 6 lw 3, \
     'speedUP_CPU680.txt' u 2:8 smooth bezier t 'Total' w l lt 1 lw 3
set output 'SpeedUP680.eps'
set ylabel 'speedup ratio'
set key right bottom width -2 height 0.2 box box
plot 'speedUP_CPU680.txt' u 2:($3/$6) smooth bezier t 'Z matrix fill' w l lt 4 lw 3, \
     'speedUP_CPU680.txt' u 2:($4/$7) smooth bezier t 'solution' w l lt 6 lw 3, \
     'speedUP_CPU680.txt' u 2:($5/$8) smooth bezier t 'total' w l lt 1 lw 3
#
set xrange [256:9072]
set xtics 256,2048,9072
set ylabel 'Time (s)'
set output 'TimeGPU580.eps'
set key left width 1. height 1 box
plot 'speedUP_CPU580.txt' u 2:6 smooth bezier t 'Solution' w l lt 4 lw 3, \
     'speedUP_CPU580.txt' u 2:7 smooth bezier t 'Z matrix fill' w l lt 6 lw 3, \
     'speedUP_CPU580.txt' u 2:8 smooth bezier t 'Total' w l lt 1 lw 3
set output 'SpeedUP580.eps'
set ylabel 'speedup ratio'
set key right bottom width -2 height 0.2 box box
plot 'speedUP_CPU580.txt' u 2:($3/$6) smooth bezier t 'Z matrix fill' w l lt 4 lw 3, \
     'speedUP_CPU580.txt' u 2:($4/$7) smooth bezier t 'solution' w l lt 6 lw 3, \
     'speedUP_CPU580.txt' u 2:($5/$8) smooth bezier t 'total' w l lt 1 lw 3
#
set xrange [256:6571]
set xtics 256,2048,6571
set output 'TimeGPU295.eps'
set key left top width 1. height 1 box
set ylabel 'Time (s)'
plot 'speedUP_CPU295.txt' u 2:6 smooth bezier t 'Solution' w l lt 4 lw 3, \
     'speedUP_CPU295.txt' u 2:7 smooth bezier t 'Z matrix fill' w l lt 6 lw 3, \
     'speedUP_CPU295.txt' u 2:8 smooth bezier t 'Total' w l lt 1 lw 3
set output 'SpeedUP295.eps'
set ylabel 'speedup ratio'
set key right bottom width -2 height 0.2 box box
plot 'speedUP_CPU295.txt' u 2:($3/$6) smooth bezier t 'Z matrix fill' w l lt 4 lw 3, \
     'speedUP_CPU295.txt' u 2:($4/$7) smooth bezier t 'solution' w l lt 6 lw 3, \
     'speedUP_CPU295.txt' u 2:($5/$8) smooth bezier t 'total' w l lt 1 lw 3
