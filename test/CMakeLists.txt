add_executable(
    testunit
    example.cpp
    filamentTest.cpp
)

target_link_libraries(
    testunit
    gtest_main
    TubeSystem_lib
   lodePNG glfw ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})


add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/testunit
)
