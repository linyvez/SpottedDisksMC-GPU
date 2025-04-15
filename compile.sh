#!/bin/bash

set -e

BUILD_TYPE="Debug"
BUILD_DIR="build-debug"
PROFILE_FLAGS=""

print_help() {
  cat << EOF
Usage: ./compile.sh [OPTIONS]

Options:
  -d        --debug-build           Compile with debug options (default)
  -o        --optimize-build        Compile with optimization (-O3)
  -i        --relwithdebinfo-build  Compile with release optimizations and debug info
  -c        --clean                 Remove all build directories and executable
  -h        --help              Show this help message
EOF
}

clean_build() {
  echo "Cleaning all builds..."
  rm -rf build-* bin/
  echo "Clean complete!"
  exit 0
}

while [[ $# -gt 0 ]]; do
  case $1 in
    -d | --debug-build)
      BUILD_TYPE="Debug"
      BUILD_DIR="build-debug"
      shift ;;
    -o | --optimize-build)
      BUILD_TYPE="Release"
      BUILD_DIR="build-release"
      shift ;;
    -i | --relwithdebinfo-build)
      BUILD_TYPE="RelWithDebInfo"
      BUILD_DIR="build-relwithdebinfo"
      PROFILE_FLAGS="-pg"
      shift ;;
    -c | --clean)
      clean_build ;;
    -h|--help)
      print_help
      exit 0 ;;
    *)
      echo "Unknown option: $1"
      print_help
      exit 1 ;;
  esac
done


echo "Configuring $BUILD_TYPE build in '$BUILD_DIR' directory..."
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DCMAKE_C_FLAGS="$PROFILE_FLAGS" ..

make


echo "Build successful! Executable located at: ./bin/aks_project"

