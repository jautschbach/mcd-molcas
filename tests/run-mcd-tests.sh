#!/bin/bash

mkdir -p mcd-tests
cd mcd-tests
home=`pwd`
test_dirs=('ucl6-f-d' 'ucl6-lmct')

for d in ${test_dirs[*]}; do
    cp -r ${MCDHOME}/tests/$d .
    cd $d
    ./run-mcd.sh
    cd $home
done

pytest -v ${MCDHOME}/tests

