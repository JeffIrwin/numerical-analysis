#!/usr/bin/env bash

set -exu

flags=""
flags+=" -Wno-tabs "
flags+=" -cpp "  # fortran pre-processing

# TODO: add debug/release arg option to switch flags
#flags+=" -O3 "

flags+=" -fbounds-check "
flags+=" -Wall -Wextra "
flags+=" -Werror "
flags+=" -Wno-unused-variable "
flags+=" -Wno-maybe-uninitialized "
flags+=" -Wno-uninitialized "

fpm test --flag "$flags"

