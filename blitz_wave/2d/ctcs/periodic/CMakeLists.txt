cmake_minimum_required(VERSION 3.0.0)
project(wave VERSION 0.1.0)

include(CTest)
enable_testing()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(wave wave_periodic.cc)


include_directories(/home/ratnesh/blitz/include)
link_directories(/home/ratnesh/blitz/lib)

target_link_libraries(wave blitz)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
