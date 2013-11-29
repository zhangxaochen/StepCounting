aux_source_directory(. lib_src_list)
add_library(${sc_lib_name} STATIC ${lib_src_list})

link_directories(${BOOST_LIBRARYDIR})
add_library(${shared_sc_lib_name} SHARED ${lib_src_list})
target_link_libraries(${shared_sc_lib_name} ${boost_filesystem_library})
