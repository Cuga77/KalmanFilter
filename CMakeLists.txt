cmake_minimum_required(VERSION 3.29)
project(kalmanFilter)

set(CMAKE_CXX_STANDARD 20)

add_executable(kalmanFilter main.cpp
        model/vector.cpp
        model/matrix.cpp
        model/Visualize.cpp
        model/KalmanFilter.cpp
        model/vector.h
        model/matrix.h
        model/Visualize.h
        model/KalmanFilter.h
)
