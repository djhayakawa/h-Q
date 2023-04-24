reset
set terminal qt font "Helvetica"
set logscale x
set logscale y
set xlabel font "Arial,25"
set ylabel font "Arial,25"
#ticsはメモリ文字
set tics font "Arial,15"
#keyは凡例
set key right outside font"Arial,15"
# ラベルの位置
set xlabel offset 0,0
set ylabel offset 0,0
#Y軸の余白
set lmargin 10
#X軸の余白
set bmargin 4
set xlabel "yr"
set ylabel "ad" offset -3,0
set title "hayashi-disk-100AU" font"Arial,15"
set format x "10^{%L}"
set format y "10^{%L}"
set yr [1e-8:10]
#set size square

plot "100_1.d" using 50:7 w l  title "0.001μm"
replot "100_2.d" using 50:7 w l  title "0.01μm"
replot "100_3.d" using 50:7 w l  title "0.1μm"
replot "100_4.d" using 50:7 w l  title "1μm"
replot "100_5.d" using 50:7 w l  title "10μm"
replot "100_6.d" using 50:7 w l  title "100μm"
replot "100_7.d" using 50:7 w l  title "0.1cm"
replot "100_8.d" using 50:7 w l  title "1cm"