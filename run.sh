#!/bin/bash

# Создаем директорию build, если её нет
mkdir -p build

# Переходим в build
cd build

# Генерируем файлы сборки
cmake ..

# Компилируем проект
make

# Проверяем успешность компиляции
if [ $? -eq 0 ]; then
    echo -e "\nЗапускаем программу...\n"
    ./kalmanFilter
else
    echo -e "\nОшибка при компиляции\n"
    exit 1
fi