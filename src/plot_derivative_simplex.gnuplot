set term png
set output "Path_Integration_Simplex.png"

set xlabel "Integration Steps" offset 0.0,-0.0,0 font "Times,10"  #linewidth 2.5
set ylabel "Fingerprint/Path Integral" offset -0.0,-0.0,0 font "Times,10" #linewidth 2.5

plot "Path_Integration_Simplex.dat" u 1:2 w l title "Eigenvalue",\
"Path_Integration_Simplex.dat" u 1:3 w l title "Path Integral"

replot
