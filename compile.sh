#!/bin/bash

mkdir -p build
cd build || exit

# shellcheck disable=SC2199
if [[ "$@" == *"--optimize-build"* ]]; then
    echo "Configuring for an optimized build..."
    cmake -DOPTIMIZE_BUILD=ON ..
else
    echo "Configuring for a standard build..."
    cmake ..
fi

make


