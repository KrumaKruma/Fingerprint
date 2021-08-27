set term png
set output "Path_Integraten.png"

plot "Path_Integration.dat" u 1:2 w l title "Eigenvalue",\
"Path_Integration.dat" u 1:3 w l title "Path Integral"

replot
