#!/usr/bin/env bash

set -exu

# See also gnuplot wrapper scripts from fynth/scripts/

gnuplot -e 'plot
	"plot-spline.txt" using 1:2 with lines,
	"plot-spline.txt" using 1:3 with lines;
	pause -1;
'

