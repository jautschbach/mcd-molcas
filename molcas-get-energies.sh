#!/usr/bin/env bash

if [ $# -ne 1 ]; then
    echo "Usage: molcas-get-energies.sh rassi.out"
    exit 1
fi


nstates=$(grep -c "SO-RASSI State" $1)
echo "# $nstates (atomic units)" > energies.txt
awk '{if(/SO-RASSI/) {print $7}}' $1 >> energies.txt

nstates_SF=$(grep "RASSI State" $1 |grep -v SO -c)
echo "# $nstates_SF (atomic units)" > energies-SF.txt
awk '{if($0 ~ /RASSI State/ && $0 !~ /SO/) {print $7}}' $1 >> energies-SF.txt


#nstates=$(grep -c "SO-RASSI" $1)
#grepa=$(echo "5+$nstates" | bc)
#energies=($(grep -A$grepa "Eigenvalues of complex Hamiltonian" $1 | tail -$nstates | awk '{print $2}'))
#echo "# $nstates (atomic units)" > energies.txt
#for i in ${energies[@]}; do
#    echo $i >> energies.txt
#done

#nstates_SF=$(grep -c "RASSI State" $1)
#nstates_SF=$(expr $nstates_SF-$nstates|bc)
#grepa=$(echo "4+$nstates_SF" | bc)
#energies=($(grep -A$grepa "SPIN-FREE ENERGIES:" $1 | tail -$nstates_SF | awk '{print $2}'))
#echo "# $nstates_SF (atomic units)" > energies-SF.txt
#for i in ${energies[@]}; do
#    echo $i >> energies-SF.txt
#done

