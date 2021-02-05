project(Hasse_diagram_utils)

file(COPY cw_decomposition_of_torus.hasse DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)

add_executable ( Hasse_diagram_from_file Hasse_diagram_from_file.cpp )
if (TBB_FOUND)
  target_link_libraries(Hasse_diagram_from_file ${TBB_LIBRARIES})
endif()

add_test(NAME Hasse_diagram_utilities_from_file COMMAND $<TARGET_FILE:Hasse_diagram_from_file>
         "cw_decomposition_of_torus.hasse")

install(TARGETS Hasse_diagram_from_file DESTINATION bin)
