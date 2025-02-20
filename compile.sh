mkdir -p build
cd build || exit
cmake ..
make
cd ..
./bin/main
rm -rf build bin

