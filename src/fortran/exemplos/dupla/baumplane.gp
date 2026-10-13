# gnuplot baumplane.gp  ->  baumplane_perfil.png e baumplane_convergencia.png
set terminal pngcairo size 900,500 font "Sans,12"
set output "baumplane_perfil.png"
set xlabel "x (mm)"; set ylabel "T (C)"; set grid; set key top left
set arrow from 28.2000,graph 0 to 28.2000,graph 1 nohead dt 2
set arrow from 28.6000,graph 0 to 28.6000,graph 1 nohead dt 2
plot "baumplane_perfil.dat" w l lw 2 lc rgb "#2a78d6" t "Analitica (6x6)", \
     "baumplane_nos.dat" w p pt 7 ps 1.2 lc rgb "#d95926" t "MEF, n = 1"
unset arrow; set output "baumplane_convergencia.png"
set logscale xy; set format y "10^{%L}"; set xlabel "h_{II} (mm)"; set ylabel "E_r (L2)"
f(x) = A*x**2; A = 1; fit [0.1:0.4] f(x) "baumplane_convergencia.dat" u 3:5 via A
plot "baumplane_convergencia.dat" u 3:5 w lp pt 7 lw 2 lc rgb "#d95926" t "MEF", \
     f(x) w l dt 2 lc rgb "#666666" t "inclinacao 2 (h^2)"
