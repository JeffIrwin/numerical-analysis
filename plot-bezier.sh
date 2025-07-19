#!/usr/bin/env bash

set -exu

# See also gnuplot wrapper scripts from fynth/scripts/

data=plot-bezier-1.txt

gnuplot -e 'plot
	"'$data'" using 1:2 with lines;
	pause -1;
'

