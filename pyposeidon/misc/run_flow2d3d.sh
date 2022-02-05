#!/bin/bash

export argfile=$1

export ncores=$2

if [ $# -gt 2 ]; then
       D3D="$3"
       LD3D="$4"
fi

#if [ -z $D3D ];then
#   echo 'no executable'
#   exit 1
#fi

    #
    # Set the directory containing delftflow.exe here
    #
exedir=$D3D/bin
libdir=$LD3D/lib

    #
    # No adaptions needed below
    #

    # Set some (environment) parameters
export LD_LIBRARY_PATH=$exedir:$libdir:$LD_LIBRARY_PATH
export PATH=$exedir:$PATH

    # Run
mpiexec -np $ncores d_hydro $argfile
