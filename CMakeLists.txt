cmake_minimum_required(VERSION 3.16)
project(ParallelCapstone C)

set(CMAKE_C_STANDARD 99)
SET(CFLAG -Wall -pedantic -std=c99 -I.. -O3 -ffast-math -Wsuggest-attribute=const)

include_directories(.)

add_library(go.o STATIC go.c)

add_executable(go_ai main_ai.c)
target_link_libraries(go_ai go.o)
target_link_options(go_ai PRIVATE -lm -fopenmp)

add_executable(go_omp main_omp.c)
target_link_libraries(go_omp go.o)
target_link_options(go_omp PRIVATE -lm -fopenmp)

add_executable(go_serial main_serial.c)
target_link_libraries(go_serial go.o)
target_link_options(go_serial PRIVATE -lm)
