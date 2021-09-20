set term png
set output "Path_Integration_Simplex.png"

plot "Path_Integration_Simplex.dat" u 1:2 w l title "Eigenvalue",\
"Path_Integration_Simplex.dat" u 1:3 w l title "Path Integral"

replot
