#!/usr/bin/env bash

set -exu

# See also gnuplot wrapper scripts from fynth/scripts/

data=plot-poly-1.txt
datd=plot-poly-data-1.txt

gnuplot -e 'plot
	"'$data'" using 1:2 with lines,
	"'$datd'" using 1:2 with circles;
	pause -1;
'

