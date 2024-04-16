cmake_minimum_required(VERSION 3.27)
project(Code)

set(CMAKE_CXX_STANDARD 17)

# Добавляем путь к файлу example.txt
#set(TEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Output/")


# Сюда включаем файлы проекта
add_executable(Code main.cpp
        LTE.h
        LTE.cpp
        Output/test1/test1_LD2e.txt
        Output/test1/test1_LD2i.txt
        Output/test1/test1_LD3e.txt
        Output/test1/test1_LD3i.txt
        Output/test1/test1_Lax.txt
        Output/test1/test1_LaxWen.txt
)

