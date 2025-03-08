#!/usr/bin/env bash

set -exu

fpm test --flag "-Wall -Wextra -Werror -Wno-tabs -Wno-maybe-uninitialized -Wno-uninitialized"

