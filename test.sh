#!/usr/bin/env bash

set -exu

flags=""
flags+=" -Wall -Wextra "
flags+=" -Werror "
flags+=" -Wno-tabs "
flags+=" -Wno-maybe-uninitialized "
flags+=" -Wno-uninitialized "

fpm test --flag "$flags"

