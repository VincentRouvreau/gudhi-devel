add_executable_with_targets(Reduced_rips_example_from_off example_reduced_rips_from_off.cpp CGAL::CGAL Eigen3::Eigen TBB::tbb)

if(TARGET CGAL::CGAL AND TARGET Eigen3::Eigen)
  # Copy a small data file next to the binary so the test below can find it.
  file(COPY "${CMAKE_SOURCE_DIR}/data/points/tore3D_300.off" DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

  add_test(NAME Reduced_rips_example_from_off COMMAND $<TARGET_FILE:Reduced_rips_example_from_off>
      "${CMAKE_CURRENT_BINARY_DIR}/tore3D_300.off")

  install(TARGETS Reduced_rips_example_from_off DESTINATION bin)
endif()
