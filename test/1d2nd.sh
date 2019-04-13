#!/bin/bash

executable="../src/heatFD ./1d2ndinput.dat"

# Diff solution with ref output from 1d 2nd order problem
$executable
h5diff .temp1d2nd.h5 ref1d2nd.h5
status=$?

if [ $status -ne 0 ]; then
	echo "Error: h5diff returns non-zero exit code"
	h5diff -v .temp1d2nd.h5 ref1d2nd.h5
	exit 1;
fi
