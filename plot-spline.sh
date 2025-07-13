#!/usr/bin/env bash

set -exu

# See also gnuplot wrapper scripts from fynth/scripts/

data=plot-spline-8.txt

gnuplot -e 'plot
	"'$data'" using 1:2 with lines,
	"'$data'" using 1:3 with lines;
	pause -1;
'

