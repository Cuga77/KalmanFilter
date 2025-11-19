cmake_minimum_required(VERSION 3.10)
project(KalmanGyroProject)

set(CMAKE_CXX_STANDARD 17)

# Находим пакеты для визуализации (если они у вас используются через Boost)
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)

# Включаем директорию model, чтобы работали #include "model/..."
include_directories(.)

# Добавляем исполняемый файл
add_executable(kalmanFilter
    main.cpp
    model/matrix.cpp
    model/matrix.h
    model/vector.cpp
    model/vector.h
    model/KalmanFilter.cpp 
    model/KalmanFilter.h
    model/Visualize.cpp
    model/Visualize.h
    model/ParameterEstimator.cpp 
    model/ParameterEstimator.h
)

# Линкуем библиотеки (Boost нужен для gnuplot-iostream)
target_link_libraries(kalmanFilter ${Boost_LIBRARIES})