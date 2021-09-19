# plot the DOS as a sum of Gaussians
# Comment out to generate a PNG file
set term png font 'Times,30' size 2500,2000
#set pm3d map corners2color c1 #interpolate 10,10
#set dgrid3d
#set contour 
set output 'Simplex_Correlation.png' 
set border linewidth 6.0
#set xlabel "Fingerprint Distance (arb. u.)" offset 0.0,-0.2,0 font "Times,30"  #linewidth 2.5
#set ylabel "Energy (arb. u.)" offset -1.5,-0.0,0 font "Times,30" #linewidth 2.5
set xlabel "Simplex FP distance" offset 0.0,-0.2,0 font "Times,30"  #linewidth 2.5
set ylabel "OM FP distance" offset -1.5,-0.0,0 font "Times,30" #linewidth 2.5
#set samples 1600; set isosamples 1600
set logscale cb
#set logscale x
#set logscale y

set xtics font "Times,30" nomirror offset -0.0,-0.2 
set ytics font "Times,30" nomirror 
set cbtics font ",30"
set key top left
set key spacing 1.0 font "Times,30.0"
#spl "FPSC10_vs_energy.dat" u 1:2:3  notitle  
set palette rgb 33,13,10
plot 'Correlation.dat' u 1:2:3 w p ps 0.5 pt 7 lc palette z notitle



replot
pause -10 "Hit any key to continue"

#pause -10 "Hit return to continue"
