cmake_minimum_required(VERSION 3.30)
project(ProjectAKS C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
file(MAKE_DIRECTORY ${BIN_DIR})

add_executable(main
        main.c
        src/periodic_boundary.c
        src/square_config.c
        src/square_move.c
        src/circle_config.c
        src/circle_move.c)

include_directories(${CMAKE_SOURCE_DIR}/include)

target_link_libraries(main m)

