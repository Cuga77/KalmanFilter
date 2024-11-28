#!/bin/bash

mkdir -p build
cd build
cmake ..
make
if [ $? -eq 0 ]; then
    echo -e "\nЗапускаем программу...\n"
    ./kalmanFilter

    rm -f CMakeCache.txt
    rm -rf CMakeFiles
    rm -f cmake_install.cmake
    rm -f Makefile
else
    echo -e "\nОшибка при компиляции\n"
    exit 1
fi