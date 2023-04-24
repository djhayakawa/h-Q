# log_sepa_mass
reset
set terminal qt font "Helvetica"
set xl "X" 
set yl "Y"
set zeroaxis ls -1;
set xl font "Arial,15"
set yl font "Arial,15"
set tics font "Arial,15"
#set log
#set format x "10^{%L}"
#set format y "10^{%L}" 
set key font"Arial,15" 
set xr [-20:20]
set yr [-20:20]
set title "hayashi-disk-10AU" font"Arial,15"

set parametric
set size square 

plot "10_1.d" using 2:3 w l  title "0.001μm"
replot "10_2.d" using 2:3 w l  title "0.01μm"
replot "10_3.d" using 2:3 w l  title "0.1μm"
replot "10_4.d" using 2:3 w l  title "1μm"
replot "10_5.d" using 2:3 w l  title "10μm"
replot "10_6.d" using 2:3 w l  title "100μm"
replot "10_7.d" using 2:3 w l  title "0.1cm"
replot "10_8.d" using 2:3 w l  title "1cm"
