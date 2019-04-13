#!/bin/bash
#-----------------------------------------------------------------------
# bootstrap utility for proj02 autotools
# source from terminal: source bootstrap.sh
#-----------------------------------------------------------------------
if [ $# -ge 2 ] ; then
	echo "Usage: source bootstrap.sh [Option: --petsc, none]"
	exit 1
fi 
echo "purging and reloading TACC"
module purge
module load TACC

echo "loading module autotools"
module load autotools

echo "Setting PKGPATH variable"
export PKGPATH=/work/00161/karl/stampede2/public

cmd="autoreconf -f -i"
echo "Bootstrapping using $cmd"

$cmd

if [ "$1" = "--petsc" ] ; then
	module load petsc
	echo "loading hdf5"
	module load hdf5
	cmd="./configure --with-masa=$PKGPATH/masa-gnu7-0.50 --with-grvy=$PKGPATH/grvy-gnu7-0.34 --with-hdf5=$TACC_HDF5_DIR --with-petsc=$TACC_PETSC_DIR"
else
	echo "adding to module path for gcc toolchain"
	export CLASSPATH=/work/00161/karl/stampede2/public
	export MODULEPATH=$CLASSPATH/ohpc/pub/modulefiles/:$MODULEPATH
	module swap intel gnu7
	echo "loading hdf5"
	module load hdf5
	cmd="./configure --with-masa=$PKGPATH/masa-gnu7-0.50 --with-grvy=$PKGPATH/grvy-gnu7-0.34 --with-hdf5=$TACC_HDF5_DIR --enable-coverage"
fi
echo "Runnig configure using $cmd"

$cmd
