cmake_minimum_required(VERSION 3.4)
project(homoligcpp LANGUAGES CXX)

add_subdirectory(pybind11)
pybind11_add_module(splicemutrcpp interval.cpp)