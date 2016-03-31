#!/bin/bash

export LD_LIBRARY_PATH=/MPATHB/resource/NCBI/CDD/usr/lib/x86_64-linux-gnu/:${LD_LIBRARY_PATH}

ln -fs /MPATHB/resource/NCBI/CDD/rpsbproc.ini .

if test $# -lt 2; then 
	/MPATHB/resource/NCBI/CDD/rpsbproc -h
else
	/MPATHB/resource/NCBI/CDD/rpsbproc $@
fi

