project(StepCounting)
cmake_minimum_required(VERSION 2.8.0)


#set(sc_src_dir ${CMAKE_CURRENT_SOURCE_DIR}/src)
#set(pugixml_src_dir ${CMAKE_CURRENT_SOURCE_DIR}/3rd/pugixml)
#set(cbuf_src_dir ${CMAKE_CURRENT_SOURCE_DIR}/3rd/circularbuffer)

#set(src_list sc_src_dir/*.cpp pugixml_src_dir/*.cpp cbuf_src_dir/*.cpp)
#include_directories(sc_src_dir pugixml_src_dir cbuf_src_dir)

find_package(Boost 1.54.0 REQUIRED filesystem)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(include)
set(sc_lib_name StepCounting_API)
set(sc_test_name StepCounting_TEST)
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(test/data)
