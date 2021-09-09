set title 'K18'
set ylabel 'Iq/I0'
set xlabel 'q(angstrom^-1)'
p 'K18.csv' u ($1)/10:2 w p pt 7 t 'exp', 'saxs_sim.dat' u 1:2 w l lw 2 t 'sim'
pause -1