#!/bin/bash

#cmake -DCMAKE_CXX_COMPILER=/c/mingw64/bin/g++.exe -DCMAKE_MAKE_PROGRAM=/c/mingw64/bin/mingw32-make.exe -G "MinGW Makefiles" .

qmake 2>/dev/null

if [[ $? -eq 127 ]]
then
  echo "remove cmake files"
  rm -rf CMakeCache CMakeFiles/

  echo "cmake"
  cmake .

  echo "make"
  make
else
  echo "qmake"
  qmake ico.pro
  make
fi