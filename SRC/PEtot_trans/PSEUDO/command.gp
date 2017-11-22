set term x11
set nolabel
set title 'Wave Functions'
set nokey
set notime
set noparametric
set size
set yrange [  -.5:  1.0]
set xrange [0:5.5]
set format xy '%.1f'
set nogrid
set xlabel 'r [a.u.]'
set xtics 0,1
set ylabel 'rR(r)' 7
set label ' true wave functions ' at 3.7,5.6
set label '_____________________' at 3.7,5.4
set label 'pseudo wave functions' at 3.7,4.6
set label '.....................' at 3.7,4.3
plot'wfct.gp' using 1:2 with lines  ,)
'wfct.gp' using 1: 4 with lines 1,)
'wfct.gp' using 1: 6 with lines 1,)
'wfct.gp' using 1: 8 with lines 1,)
'wfct.gp' using 1:10 with lines 1,)
'wfct.gp' using 1: 3 with lines 5,)
'wfct.gp' using 1: 5 with lines 5,)
'wfct.gp' using 1: 7 with lines 5,)
'wfct.gp' using 1: 9 with lines 5,)
'wfct.gp' using 1:11 with lines 5,)
0 ti ' ' w li 5
pause -1
#
set title 'Wave Function Transforms'
set format y '%.2f'
set xtics 0,5
set yrange [-1.5:1.5]
set xrange [0:10]
set xlabel 'q [1/a.u.]'
set ylabel 'R(q)' 7
set nolabel
plot'fwfct.gp' using 1:2 with lines 4,)
'fwfct.gp' using 1: 3 with lines 4,)
'fwfct.gp' using 1: 4 with lines 4,)
'fwfct.gp' using 1: 5 with lines 4,)
'fwfct.gp' using 1: 6 with lines 4,)
0 ti ' ' w li 1
pause -1
#
set title 'Pseudopotentials'
set nolabel
set key 3.5,   -5.0
set yrange [ -10: 5]
set xrange [0:4]
set xlabel 'r [a.u.]'
set ylabel 'V(r) [Ry]' 5
set format xy '%.1f'
set xtics 0,1.0
plot )
'ppot.gp' using 1: 2 ti 's - nonlocal' w li 1,)
'ppot.gp' using 1: 3 ti 'p - nonlocal' w li 2,)
'ppot.gp' using 1: 4 ti 'p - nonlocal' w li 3,)
'ppot.gp' using 1: 5 ti 'd - nonlocal' w li 4,)
'ppot.gp' using 1: 6 ti 'd - nonlocal' w li 5,)
-12.00/x ti 'Z - ion' w li 5, 0 ti ' ' w li 1
pause -1
#
set nolabel
set key 10,-0.4
set size 2.8/5.,3/3.
set title 'Fourier Transforms'
set xrange [0:12]
set yrange [-1.5:1.5]
set xlabel 'q [1/a.u.]'
set ylabel 'q^2/4piZion V(q) [Ry]' 4
set xtics 0,5.0
plot )
'fppot.gp' using 1: 2 ti 's - nonlocal' w li 1,)
'fppot.gp' using 1: 3 ti 'p - nonlocal' w li 2,)
'fppot.gp' using 1: 4 ti 'p - nonlocal' w li 3,)
'fppot.gp' using 1: 5 ti 'd - nonlocal' w li 4,)
'fppot.gp' using 1: 6 ti 'd - nonlocal' w li 5,)
0 ti ' ' w li 1
