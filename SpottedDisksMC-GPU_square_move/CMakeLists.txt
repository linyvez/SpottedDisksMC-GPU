cmake_minimum_required(VERSION 3.29.6)
project(ProjectAKS C)

set(CMAKE_C_STANDARD 11)

option(OPTIMIZE_BUILD "Build with optimization flags" OFF)

if(OPTIMIZE_BUILD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_executable(main
        src/periodic_boundary.c
        src/square_config.c
        src/square_move.c
        main.c
)

include_directories(${CMAKE_SOURCE_DIR}/include)

target_link_libraries(main m)

