# OpenGL & GLEW library
find_package(OpenGL)
ADD_DEFINITIONS(-DGLEW_STATIC)

# source files
file(GLOB SOURCES ./*.cpp)
file(GLOB HEADERS ./*.h)
file(GLOB SHADERS ./*.vert ./*.frag)

# libraries
add_library(TubeSystem_lib "")


# Make sure the textures and shaders are available
set(TEXTURE_PATH ${CMAKE_SOURCE_DIR}/textures CACHE PATH "location of texture images")
set(SHADER_PATH  ${CMAKE_SOURCE_DIR}/src      CACHE PATH "location of shader programs")
add_definitions("-DTEXTURE_PATH=\"${TEXTURE_PATH}\"")
add_definitions("-DSHADER_PATH=\"${SHADER_PATH}\"")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# executable
add_executable(TubeSystem ${HEADERS} ${SOURCES})
target_link_libraries (TubeSystem Eigen3::Eigen)
target_include_directories(TubeSystem SYSTEM PUBLIC ${GLEW_INCLUDE_DIRS})
# Note: target_link_libraries(glfw) should actually bring in the necessary header files.
target_link_libraries(TubeSystem lodePNG glfw ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})



target_sources(
    TubeSystem_lib
    PRIVATE
        filament.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/filament.h
)

target_include_directories(TubeSystem_lib SYSTEM PUBLIC ${GLEW_INCLUDE_DIRS})
# Note: target_link_libraries(glfw) should actually bring in the necessary header files.
target_link_libraries(TubeSystem_lib lodePNG glfw ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})
