 # cmake 3.15 has problems with the matlab interface (known issue) that were resolved in 3.16
 # versions before 3.15 might work but were not tested
cmake_minimum_required(VERSION 3.16)

project(pushpull)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/third_party)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/third_party/common)
ADD_SUBDIRECTORY(third_party)

find_package(CGAL REQUIRED)
find_package(Matlab)
matlab_add_mex(NAME blueNoise SRC
    src/myInterface.cpp
    src/ppo.h
    src/ppo.cpp)
target_include_directories(blueNoise PRIVATE src)
target_link_libraries(blueNoise CGAL::CGAL)
