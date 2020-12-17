add_executable(
    testunit
    example.cpp
    filamentTest.cpp
)

target_link_libraries(
    testunit
    gtest_main
    gmock_main
    TubeSystem_lib
   lodePNG glfw ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})


add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/testunit
)

# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE "Debug")
# endif()

# add_compile_options(-fno-omit-frame-pointer) 
# add_compile_options(-fsanitize=address) 
# add_compile_options(-lasan)

# # What we want: LDFLAGS="-fsanitize=address"
# link_libraries("-fsanitize=address") 