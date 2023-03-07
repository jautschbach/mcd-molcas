#!/bin/bash

pushd . > /dev/null
curdir="$PWD"
cd "`dirname "${BASH_SOURCE[0]}"`" > /dev/null
MCDHOME="$PWD"
MCDBIN=$MCDHOME
popd > /dev/null

export MCDHOME MCDBIN

