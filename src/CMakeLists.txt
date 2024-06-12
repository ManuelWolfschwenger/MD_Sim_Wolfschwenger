cmake_minimum_required (VERSION 3.0.0)
project(test)

#use the c++ 14 standard
set(CMAKE_CXX_STANDARD 14)

#Bring the headers, such as game.h into the project
include_directories( lib )

#Find all cpp files
file(GLOB SOURCES "*.cpp")

#Building SeWorms
add_executable(test ${SOURCES})