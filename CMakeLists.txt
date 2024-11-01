cmake_minimum_required(VERSION 3.29)
project(kalmanFilter)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost REQUIRED COMPONENTS iostreams)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost libraries: ${Boost_LIBRARIES}")
else()
    message(FATAL_ERROR "Could not find Boost libraries")
endif()

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

target_link_libraries(kalmanFilter
        ${Boost_LIBRARIES}
)
