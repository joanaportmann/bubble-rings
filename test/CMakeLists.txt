add_executable(
    testunit
    example_add.cpp
    )

target_link_libraries(
    testunit
    gtest_main
        )

add_test(
  NAME
    unit
  COMMAND
    ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/testunit
  )
