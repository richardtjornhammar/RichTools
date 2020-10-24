cmake_minimum_required(VERSION 3.7)
project(My_RT_Project)

include_directories(src)
file(GLOB SOURCES "src/*.cc")

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(richtool ${SOURCES})

find_package(GSL REQUIRED)
target_link_libraries(richtool GSL::gsl GSL::gslcblas)
