#!/bin/bash
#SBATCH -N 2
#SBATCH -n 32
#SBATCH -o 1d4thbatch.log
#SBATCH -J heatFD
#SBATCH -p skx-normal
#SBATCH -A cse38018
#SBATCH -t 01:00:00

echo "Master Host = "`hostname`
echo "PWD = "`pwd`
echo "Time run = "`date`
echo
source bootstrap.sh
make
./src/heatFD ./inputs/input1d4th20.dat
./src/heatFD ./inputs/input1d4th40.dat
./src/heatFD ./inputs/input1d4th80.dat
./src/heatFD ./inputs/input1d4th160.dat
./src/heatFD ./inputs/input1d4th320.dat
