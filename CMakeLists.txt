 # cmake 3.15 has problems with the matlab interface (known issue) that were resolved in 3.16
 # versions before 3.15 might work but were not tested
cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)

project(pushpull)

find_package(CGAL REQUIRED)
find_package(Matlab)

find_package(Python COMPONENTS Development NumPy) # python is required for matplotlibcpp

matlab_add_mex(NAME blueNoiseTwoTiles SRC
    src/myInterface.cpp
    src/pointSet.h
    src/pointSet.cpp
    src/ppo.h
    src/ppo.cpp
    src/utility.h)
target_include_directories(blueNoiseTwoTiles PRIVATE src)
target_link_libraries(blueNoiseTwoTiles CGAL::CGAL)

add_executable(blueNoiseStandalone
    src/main.cpp
    src/pointSet.h
    src/pointSet.cpp
    src/ppo.h
    src/ppo.cpp
    src/utility.h
    thirdParty/matplotlibcpp.h)
target_include_directories(blueNoiseStandalone PRIVATE src thirdParty)
target_link_libraries(blueNoiseStandalone CGAL::CGAL)

target_include_directories(blueNoiseStandalone PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(blueNoiseStandalone Python::Python Python::NumPy)
