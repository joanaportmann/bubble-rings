file(GLOB SRCS *.cpp)

ADD_EXECUTABLE(testunit ${SRCS})

TARGET_LINK_LIBRARIES(testunit
    libgtest
    libgmock
)

add_test(NAME unit
         COMMAND ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/testunit)

         