cmake_minimum_required(VERSION 3.10)
project(kalmanFilter)

set(CMAKE_CXX_STANDARD 17)

add_executable(kalmanFilter main.cpp
        model/vector.cpp
        model/matrix.cpp
        model/Visualize.cpp
        model/KalmanFilter.cpp
        model/vector.h
        model/matrix.h
        model/Visualize.h
        model/KalmanFilter.h
        model/MotionIntegrator.cpp
        model/MotionIntegrator.h
)

target_link_libraries(kalmanFilter)