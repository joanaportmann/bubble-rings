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
set(VL_BASE "../lib/VisualizationLibrary")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${VL_BASE}/_INSTALL/cmake")
set(VL_ROOT "${VL_BASE}/_INSTALL")
add_definitions("-DTEXTURE_PATH=\"${TEXTURE_PATH}\"")
add_definitions("-DSHADER_PATH=\"${SHADER_PATH}\"")

# XXX: Questionable...
# include_directories(${VL_ROOT})
# add_subdirectory(${VL_ROOT})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VL COMPONENTS VLCore VLGraphics VLVG VLVolume VLMolecule REQUIRED)

# executable
add_executable(TubeSystem ${HEADERS} ${SOURCES})
include_directories(${VL_INCLUDE_DIRS})
target_include_directories(TubeSystem SYSTEM PUBLIC ${GLEW_INCLUDE_DIRS})
# Note: target_link_libraries(glfw) should actually bring in the necessary header files.
target_link_libraries(TubeSystem lodePNG glfw ${GLEW_LIBRARIES} ${OPENGL_LIBRARIES})
target_link_libraries (TubeSystem Eigen3::Eigen)
target_link_libraries(TubeSystem ${VL_LIBRARIES})



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
