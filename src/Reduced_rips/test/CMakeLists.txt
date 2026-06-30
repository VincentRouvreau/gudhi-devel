add_executable_with_targets(Reduced_rips_test_unit test_reduced_rips.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  include(GUDHI_boost_test)
  # Copy the data file used by the OFF-file integration test into the test working directory.
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
  gudhi_add_boost_test(Reduced_rips_test_unit)
endif()
