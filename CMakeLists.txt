# Minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(HP_model VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directory for header files
include_directories(src)

# Define the source files
set(SOURCE_FILES
    src/HP_model.cpp
    src/sequence_reader.cpp
    src/run.cpp
)

# Add the executable
add_executable(HP_model ${SOURCE_FILES})

# Copy input.txt to the build directory for runtime
add_custom_command(
    TARGET HP_model POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/input $<TARGET_FILE_DIR:HP_model>/input
)

# Create the output and conformations directories in the build folder
add_custom_command(
    TARGET HP_model POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:HP_model>/output
    $<TARGET_FILE_DIR:HP_model>/output/conformations
)

# Copy Python scripts to the build directory
add_custom_command(
    TARGET HP_model POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/2D_animator.py
    $<TARGET_FILE_DIR:HP_model>/2D_animator.py
)

add_custom_command(
    TARGET HP_model POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/3D_animator.py
    $<TARGET_FILE_DIR:HP_model>/3D_animator.py
)
