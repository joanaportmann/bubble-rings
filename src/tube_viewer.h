
//=============================================================================
#ifndef TUBE_VIEWER_H
#define TUBE_VIEWER_H
//=============================================================================

#include "sphere.h"
#include "gl.h"
#include "glfw_window.h"
#include "shader.h"
#include "path.h"
#include "frame.h"
#include "tube.h"
#include "filament.h"
#include "planet.h"

/// OpenGL viewer that handles all the rendering for us
class Tube_viewer : public GLFW_window
{
public:
  /// default constructor
  /// \_title the window's title
  /// \_width the window's width
  /// \_height the window's height
  Tube_viewer(const char *_title, int _width, int _height);

protected:
  /// function that is called on the, creation of the widget for the initialisation of OpenGL
  virtual void initialize();

  /// resize function - called when the window size is changed
  virtual void resize(int width, int height);

  /// paint function - called when the window should be refreshed
  virtual void paint();

  virtual void timer();

  /// keyboard interaction
  virtual void keyboard(int key, int scancode, int action, int mods);

  /// function that draws the bubble ring
  /// \param _projection the projection matrix for the scene
  /// \param _view the view matrix for the scene
  void draw_scene(mat4 &_projection, mat4 &_view, vec3 &eye);

  void drawCircle(const std::vector<vec3> pts, float radius);

private:

  // Draw sphere and map underwater image
  Planet background;
  Sphere unit_sphere_;
  vec3 originalFirstVertex;

  // filament object
  Filament filament;

  bool timer_active_ = false;

  // tube object
  Tube tube;


  float thickness = 0.12;
  float circulation = 4;
  float tension = 0;
  float alpha = 0.5;
  float length = 0.1;
  bool renderOnlyPolygon = false;
  bool modifyThickness = true;
  bool backgroundOn = true;
  bool showCoordinateAxis = false;
  bool rungeKutta4 = true;
  bool euler = false;

  /// default color shader (renders only texture)
  Shader color_shader_;

  /// simple shader for visualizing curves (just using solid color).
  Shader background_shader_;

  /// Shader
  Shader test_tube_shader_;

  Shader solid_color_shader_;

  Shader reflection_shader_;

  /// state whether the rendering should be in color or not
  bool greyscale_;

  /// rotation in x direction
  float x_angle_;
  /// rotation in y direction
  float y_angle_;

  /// eye's distance in radii from the observed point
  float dist_factor_;

  std::vector<vec3> verticesOfTube;

  Frame coordinateAxis;
  vec3 center_of_coordinatesystem;


  /// the field of view for the camera
  float fovy_;
  /// the near plane for the virtual camera
  float near_;
  /// the far plane for the virtual camera
  float far_;

  /// current viewport dimension
  int width_, height_;
};

//=============================================================================
#endif
//=============================================================================
