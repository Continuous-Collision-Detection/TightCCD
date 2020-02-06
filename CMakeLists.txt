cmake_minimum_required(VERSION 3.1)
project(TightCCD)

################################################################################
# Library
################################################################################

add_library(${PROJECT_NAME}
  src/bsc_tightbound.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC src)

################################################################################
# Compiler options
################################################################################

# Use C++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
