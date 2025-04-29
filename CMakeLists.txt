cmake_minimum_required(VERSION 3.30)
project(particle_simulation C)

set(CMAKE_C_STANDARD 11)

add_executable(
        ${CMAKE_PROJECT_NAME}
        src/periodic_boundary.c
        src/square_config.c
        src/square_move.c
        src/circle_config.c
        src/circle_move.c
        main.c
        src/patch.c)

include_directories(${CMAKE_SOURCE_DIR}/include)

target_link_libraries(${CMAKE_PROJECT_NAME} m)

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
)