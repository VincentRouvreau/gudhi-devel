project(Persistence_matrix_benchmark)

add_executable(persistence_matrix_benchmark persistence_matrix_benchmark.cpp)
if(TARGET TBB::tbb)
  target_link_libraries(persistence_matrix_benchmark TBB::tbb)
endif()
