cmake_minimum_required(VERSION 2.8.5)
project(My_MD_Project)

include_directories(src)
file(GLOB SOURCES "src/*.c")

add_executable(stochMD ${SOURCES})
target_link_libraries(stochMD -lgsl -lblas)
