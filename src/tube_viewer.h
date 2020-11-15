//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================
#ifndef TUBE_VIEWER_H
#define TUBE_VIEWER_H
//=============================================================================

#include "gl.h"
#include "glfw_window.h"
#include "shader.h"
#include "path.h"
#include "frame.h"
#include "tube.h"
#include "filament.h"

/// OpenGL viewer that handles all the rendering for us
class Tube_viewer : public GLFW_window
{
public:
  /// default constructor
  /// \_title the window's title
  /// \_width the window's width
  /// \_height the window's height
  Tube_viewer(const char *_title, int _width, int _height);

  // Control polygon for circle to be contoured
  std::vector<FilamentPoint> control_polygon_ = {
      {{-5.2, 3.2, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{-4.4, 3.94, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{-3.7, 4.4, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{-2.85, 4.68, 0.0}, 0.4, 4, vec3(0, 0, 0)},
      {{-1.88, 4.72, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{-0.43, 4.62, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{0.2, 4.18, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{0.87, 3.68, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{1.09, 3.24, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{1.2, 2.7, 0.0}, 0.3, 4, vec3(0, 0, 0)},
      {{1.45, 2.08, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{1.5, 1.32, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{1.3, 0.2, 0.0}, 0.1, 4, vec3(0, 0, 0)},
      {{-0.6, -1.45, 0.0}, 0.5, 4, vec3(0, 0, 0)},
      {{-2.8, -1.5, 0.0}, 0.3, 4, vec3(0, 0, 0)},
      {{-4.95, -0.7, 0.0}, 0.3, 4, vec3(0, 0, 0)},
      {{-5.5, 1.6, 0.5}, 0.1, 3, vec3(0, 0, 0)}};

protected:
  /// function that is called on the creation of the widget for the initialisation of OpenGL
  virtual void initialize();

  /// resize function - called when the window size is changed
  virtual void resize(int width, int height);

  /// paint function - called when the window should be refreshed
  virtual void paint();

  virtual void timer();

  /// keyboard interaction
  virtual void keyboard(int key, int scancode, int action, int mods);

  /// function that draws the planet system
  /// \param _projection the projection matrix for the scene
  /// \param _view the view matrix for the scene
  void draw_scene(mat4 &_projection, mat4 &_view);

  void drawCircle(const std::vector<vec3> pts, float radius);

private:
  // filament object
  Filament filament;

  // tube object
  Tube tube;

  /// default color shader (renders only texture)
  Shader color_shader_;
  /// phong shader (renders texture and basic illumination)
  Shader phong_shader_;

  // sun shader (renders the sun: texture plus an optional shimmer effect)
  Shader tube_shader_;

  /// simple shader for visualizing curves (just using solid color).
  Shader solid_color_shader_;

  /// Shader
  Shader test_tube_shader_;

  /// interval for the animation timer
  bool timer_active_;

  /// state whether the rendering should be in color or not
  bool greyscale_;

  /// rotation in x direction
  float x_angle_;
  /// rotation in y direction
  float y_angle_;

  /// eye's distance in radii from the observed point
  float dist_factor_;

  std::vector<vec3> verticesOfTube;

  Path ship_path_renderer_, ship_path_cp_renderer_, circle2, circle1;
  Frame ship_path_frame_;

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
