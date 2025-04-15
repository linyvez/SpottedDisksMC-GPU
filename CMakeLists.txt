cmake_minimum_required(VERSION 3.30)
project(aks_project C)

set(CMAKE_C_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(
        ${CMAKE_PROJECT_NAME}
        src/periodic_boundary.c
        src/square_config.c
        src/square_move.c
        src/circle_config.c
        src/circle_move.c
        main.c)

include_directories(${CMAKE_SOURCE_DIR}/include)

target_link_libraries(${CMAKE_PROJECT_NAME} m)

