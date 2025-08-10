#!/usr/bin/env bash

set -exu

# Configuration is set in Doxyfile.  The output directory, set in Doxyfile, is
# "./doxy/", which can be cleaned with `rm -rf doxy`
#
# To view the docs after generating them with this script, just open the index
# html in a browser:
#
#     ./doxy/html/index.html

## Doxygen v1.9.1 is the version installed on Ubuntu as of 2025-08-10, but it lacks features like HTML_COLORSTYLE to set dark mode
#doxygen-1.9.1


# Built from Doxygen github commit f1ec8b74 using their Dockerfile.  This is
# doxygen v1.15.0.  I probably didn't need docker if I just apt installed bison
doxygen

