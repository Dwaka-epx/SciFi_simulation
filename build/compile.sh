#!/bin/bash 

#printenv | grep "GEANT4_DIR"
#GEANT4_DIR=/home/t2k/ogawat/myt2kwork/software_cos7/geant4.10.1.3_install/

rm -rf CMake* cmake_install.cmake
cmake -DGeant4_DIR=$GEANT4_DIR/lib64 ../
make -j 5
