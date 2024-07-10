#!/bin/bash

####################### export ######################
export Project=cas22scf_twisted-ethene
export HomeDir=$PWD
export WorkDir=/tmp
export MOLCAS_WORKDIR=$WorkDir
export MOLCAS_MEM=12000
export MOLCAS_PRINT=3
######################################################

echo "start up"
echo "HOME=$HOME"
echo "HOSTNAME=$HOSTNAME"

export | grep MOLCAS

cd $WorkDir

echo "Working Directory=`pwd`"

###################### run molcas ######################

$MOLCAS/pymolcas  $HomeDir/${Project}.inp >  $HomeDir/${Project}.log

##################### copy #############################
mv *.RasOrb $HomeDir
mv *.JobIph $HomeDir
mv spin*.txt $HomeDir
mv angmom*.txt $HomeDir
mv dipole*.txt $HomeDir
mv velocity*.txt $HomeDir
mv quadrupole*.txt $HomeDir
mv *.rasscf.h5 $HomeDir
mv *.lus $HomeDir

exit
#######################################################

