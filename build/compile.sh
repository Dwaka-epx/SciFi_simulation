#!/bin/bash 

#printenv | grep "GEANT4_DIR"
#GEANT4_DIR=/home/t2k/ogawat/myt2kwork/software_cos7/geant4.10.1.3_install/
GEANT4_DIR=/opt/MyGeant4/geant4.10.05.p01-install/

rm -rf CMakeFiles CMakeCache.txt cmake_install.cmake
cmake -DGeant4_DIR=$GEANT4_DIR/lib64/Geant4-10.5.1/ ../
make -j 
