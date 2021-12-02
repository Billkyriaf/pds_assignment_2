cmake_minimum_required(VERSION 3.20)
project(Distance_sort C)

set(CMAKE_C_STANDARD 99)

include_directories(src/data)
include_directories(src/logs)

include_directories(${MPI_INCLUDE_PATH})
find_package(MPI REQUIRED)

add_executable(Distance_sort src/main.c src/data/points.c src/data/points.h src/data/read_points.c src/data/read_points.h src/logs/logging.c src/logs/logging.h)
target_link_libraries(Distance_sort ${MPI_LIBRARIES})