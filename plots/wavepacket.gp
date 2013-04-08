set terminal postscript enhanced color solid "Helvetica" 22
set output '| ps2pdf - wavepacket.pdf'
set samples 1000
unset xtics
unset ytics
unset key
set xlabel "x"
set ylabel "Re {/Symbol y}"
plot [-0.4:0.4] [-1.25:1.25] exp(-50*x**2)*cos(100*x)
