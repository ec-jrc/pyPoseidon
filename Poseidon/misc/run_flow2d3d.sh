#!/bin/bash
    #
    # This script is an example for running Delft3D-FLOW
    # Adapt and use it for your own purpose
    #
    # adri.mourits@deltares.nl
    # 27 Dec 2010
    # 
    #
    # This script starts a single-domain Delft3D-FLOW computation on Linux
    #


    #
    # Set the config file here
    # 
#argfile=config_flow2d3d.xml
argfile=config_d_hydro.xml


export D3D=$1
export ncores=$2

    #
    # Set the directory containing delftflow.exe here
    #
#export D3D_HOME=$D3D/bin/$ARCH
exedir=$D3D/flow2d3d/bin
libdir=$D3D/flow2d3d/bin
scrlib=$D3D/../../src/lib
 
    #
    # No adaptions needed below
    #

    # Set some (environment) parameters
export LD_LIBRARY_PATH=$exedir:$libdir:$scrlib:$LD_LIBRARY_PATH 

    # Run
mpiexec -np $ncores $exedir/d_hydro.exe $argfile
