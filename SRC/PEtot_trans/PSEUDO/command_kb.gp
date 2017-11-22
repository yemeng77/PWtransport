set term x11
set nolabel
set title 'Kleinman & Bylander Operators'
set nokey
set notime
set noparametric
set size
set autoscale y
set xrange [0: 3.0]
set format xy '%.1f'
set nogrid
set xlabel 'r [a.u.]'
set xtics
set key
set ylabel 'V(r) [(Ry/a^3)^0.5]' 7

'kbop.gp' using  1: 2 ti 'p - Kl & By '
'kbop.gp' using  3: 4 ti 'd - Kl & By '
0 ti ' ' w li 1
pause -1
#
set nolabel
set title 'Kleinman & Bylander Transforms'
set format y '%.2f'
set xtics
set key
set autoscale y
set autoscale x
set xlabel 'q [1/a.u.]'
set ylabel 'Bessel V(q) [Ry^0.5]' 7
set nolabel

'kbtr.gp' using  1: 2 ti 'trans p K&B '
'kbtr.gp' using  3: 4 ti 'trans d K&B '
0 ti ' ' w li 1
pause -1
#
set title 'Kleinman & Bylander Area'
set nolabel
set key
set autoscale y
set xrange [0: 3.0]
set xlabel 'r [a.u.]'
set ylabel 'V(r) [Ry]' 5
set format xy '%.1f'
set xtics

'kbar.gp' using  1: 2 ti 'ppt  '
'kbar.gp' using  3: 4 ti 'pdt  '
0 ti ' ' w li 1
