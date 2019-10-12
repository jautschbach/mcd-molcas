#!/bin/bash

# prepare gnuplot graph from "plotspec" output,

# adjust this for the location of your gnuplot binary:
PL=/usr/bin/gnuplot
#PL=/user/fvgendro/bin/bin/gnuplot

ylabel="MCD terms"  # default y2 label
scint=1    # default scaling of stick intensities
y2scale=1000 # default inverse scaling of y2 axis versys y
shiftspec=0  # default shift of excitation energies
lscale="set nologscale y"
mcdunit="{/Symbol De} / M^{-1} cm^{-1} T^{-1}"

# ----------------------
# parse the command line
# ----------------------
 
#set -- ` getopt lr:y:s:m: $*`
#if [ $? != 0 ]
#then
#        exit 2
#fi
#while [ $1 != -- ]
while getopts "ltr:y:s:m:" flag
do    
    case $flag in           
        r)
            scint="${OPTARG}"
            ;;
        l)
            lscale="set logscale y"
            ;;
	t)
            mcdunit="[{/Symbol q}] / deg cm^2 dmol^{-1} G^{-1}"
            ;;
	y)
            ylabel="${OPTARG}"
            ;;
	s)
            y2scale="${OPTARG}"
            ;;     
	m)
            shiftspec="${OPTARG}"
            ;;
#	*)
#	    echo 'unknown option'; exit
     esac
done

#echo $ylabel

shift $((OPTIND-1))

if test $scint != "1"
then
  yarg="$ylabel (x $scint) / Debye^2"
else
  yarg="$ylabel / Debye^2"
fi


# ---------------------------
# check for filename argument
# ---------------------------

if test -z "$1"
then
echo " "
echo " output filename specification is missing ... aborting"
echo " "
echo " usage: graphit filespec"
echo " "
exit 1
else

# --------------------------------
# check if input files are present 
# --------------------------------

if test -r "graph.dat"
then

if test -r "impulses.dat"
then

#echo $scosc
#echo $scrot
#echo $shiftev

# -------------------------
# plot on energy scale (wavenumbers)
# -------------------------

$PL << eor
set loadpath '/usr/bin/gnuplot'
c2f(y) = $y2scale*y
f2c(y) = y/$y2scale
set link y2 via f2c(y) inverse c2f(y)

#set xrange[10:60]
#set yrange[-4:4]

set xzeroaxis
set term postscript eps enhanced "Helvetica" 20
set xlabel "{/Symbol D}E / 10^3 cm^{-1}"
set ylabel "$mcdunit"
set y2label "$yarg"
set nokey 
set xzeroaxis
set ytics nomirror
set y2tics nomirror
set title "$1"



set out "$1-cm.eps"
plot "graph.dat" using ((\$1 + ($shiftspec))/1000):2   w l lt 1 lw 4 , \
"impulses.dat" using ((\$1 + ($shiftspec))/1000):(\$3 * $scint) axes x1y2  w imp lt 1 lw 2
set out
#
eor


else # test -r impulses.dat
echo " "
echo " input file impulses.dat not found ... aborting"
echo " "
exit 2
fi # test -r impulses.dat

else # test -r graph.dat
echo " "
echo " input file graph.dat not found ... aborting"
echo " "
exit 3
fi # test -r graph.dat

fi # test -z $1
rm -f tmp.$$
exit 0


