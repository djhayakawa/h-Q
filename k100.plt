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
set xr [-200:200]
set yr [-200:200]
set title "hayashi-disk-100AU" font"Arial,15"

set parametric
set size square 

plot "100_1.d" using 2:3 w l  title "0.001μm"
replot "100_2.d" using 2:3 w l  title "0.01μm"
replot "100_3.d" using 2:3 w l  title "0.1μm"
replot "100_4.d" using 2:3 w l  title "1μm"
replot "100_5.d" using 2:3 w l  title "10μm"
replot "100_6.d" using 2:3 w l  title "100μm"
replot "100_7.d" using 2:3 w l  title "0.1cm"
replot "100_8.d" using 2:3 w l  title "1cm"