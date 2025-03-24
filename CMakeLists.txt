cmake_minimum_required(VERSION 3.10)
project(ParticleSimulation C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Ensure the bin directory exists
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# Add include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add executable
add_executable(main
        src/generate_config.c
        src/periodic_boundary.c
)

# Link math library
target_link_libraries(main m)