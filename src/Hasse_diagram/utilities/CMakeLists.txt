project(Hasse_diagram_utilities)

add_executable ( Hasse_diagram_from_file Hasse_diagram_from_file.cpp )
if (TBB_FOUND)
  target_link_libraries(Hasse_diagram_from_file ${TBB_LIBRARIES})
endif()

add_test(NAME Hasse_diagram_from_torus_file COMMAND $<TARGET_FILE:Hasse_diagram_from_file>
  "${CMAKE_SOURCE_DIR}/data/hasse_diagram/cw_decomposition_of_torus.hasse")

install(TARGETS Hasse_diagram_from_file DESTINATION bin)
