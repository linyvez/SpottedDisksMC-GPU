mkdir -p build
cd build || exit
cmake ..
make
cd ..
rm -rf build
