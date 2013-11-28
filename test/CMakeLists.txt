aux_source_directory(. test_src_list)

#find_package(Boost REQUIRED COMPONENTS filesystem )	#顶层文件夹已经有这句
#LINK_DIRECTORIES(${Boost_LIBRARY_DIR})
LINK_DIRECTORIES(${BOOST_LIBRARYDIR})

add_executable(${sc_test_name} ${test_src_list})

target_link_libraries(${sc_test_name} ${sc_lib_name})
target_link_libraries(${sc_test_name} ${Boost_FILESYSTEM_LIBRARY})
#target_link_libraries(${sc_test_name} ${BOOST_LIBRARYDIR})		#×

