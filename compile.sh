mkdir -p build
cd build || exit
cmake ..
make
cd ..
./bin/ProjectAKS
rm -rf build bin

