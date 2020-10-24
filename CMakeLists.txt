cmake_minimum_required(VERSION 2.8.5)
project(My_RT_Project)

include_directories(src)
file(GLOB SOURCES "src/*.cc")

#SET(GCC_COVERAGE_LINK_FLAGS "-lgsl -lblas")
#SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

add_executable(richtool ${SOURCES})
target_link_libraries(richtool -lgsl -lblas)
