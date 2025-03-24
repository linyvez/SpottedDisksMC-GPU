#!/bin/bash

set -e  # Stop script if any command fails

# Create build directory
mkdir -p build && cd build

# Run CMake and build
cmake ..
make

# Return to project root and execute
cd ..
./bin/main "$@"

# Cleanup
rm -rf build bin