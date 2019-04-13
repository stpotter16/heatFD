#!/bin/bash

executable="../src/heatFD"

# verify executible exists
if [ ! -x "$executable" ]; then
	echo "Error: expecting executable -> $executable"
	exit 1;
fi
