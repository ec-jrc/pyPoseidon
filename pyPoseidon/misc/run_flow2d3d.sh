#!/bin/bash
argfile=config_d_hydro.xml

export ncores=$1

source activate pyPoseidon
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

d3d=$D3D

if [ $# -eq 2 ]; then
#    if [[ -z "${D3D}" ]]; then
       d3d="$2"
#    fi
fi

if [ -z $d3d ];then
   echo 'no executable'
   exit 1 
fi

    #
    # Set the directory containing delftflow.exe here
    #
exedir=$d3d/flow2d3d/bin
libdir=$d3d/flow2d3d/bin
scrlib=$d3d/../../src/lib
 
    #
    # No adaptions needed below
    #

    # Set some (environment) parameters
export LD_LIBRARY_PATH=$exedir:$libdir:$scrlib:$LD_LIBRARY_PATH 

    # Run
mpiexec -np $ncores $exedir/d_hydro.exe $argfile
