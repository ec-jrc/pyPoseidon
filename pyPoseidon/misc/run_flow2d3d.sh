#!/bin/bash
argfile=config_d_hydro.xml

export ncores=$1

if [ $# -gt 2 ]; then
	source activate $3 #relevant conda env
	export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
fi

if [ $# -gt 1 ]; then
       d3d="$2"
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
