cmake_minimum_required(VERSION 3.28.3)
project(particle_simulation C CUDA)

set(CMAKE_C_STANDARD 11)
set(CUDA_ROOT $ENV{CUDA_HOME})

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)


if (CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

find_path(CURAND_INCLUDE_DIR
        NAMES curand.h
        HINTS ${CUDA_ROOT}/include)

find_library(CURAND_LIBRARY
        NAMES curand
        HINTS ${CUDA_ROOT}/lib/x64
)

add_executable(
        ${CMAKE_PROJECT_NAME}
        src/periodic_boundary.c
        src/square_config.c
        src/square_move.cu
        include/square_move.cuh
        include/periodic_boundary.h
        include/patch.h
        include/general_config.h
        include/square_config.h
        main.c
        src/patch.c
        src/shared_utilities.c
        include/shared_utilities.h
        src/circle_config.c
        src/circle_move.cu
        include/circle_config.h
        include/circle_move.cuh
)



set_target_properties(${PROJECT_NAME} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
)


target_include_directories(${PROJECT_NAME} PRIVATE
        ${CUDA_ROOT}/include
        ${CURAND_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/include
)


target_link_libraries(${PROJECT_NAME}
        PRIVATE
        ${CURAND_LIBRARY}
)


install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
)