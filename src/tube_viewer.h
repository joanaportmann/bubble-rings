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



/// OpenGL viewer that handles all the rendering for us
class Tube_viewer : public GLFW_window
{
public:

    /// default constructor
    /// \_title the window's title
    /// \_width the window's width
    /// \_height the window's height
    Tube_viewer(const char* _title, int _width, int _height);


protected:

    /// function that is called on the creation of the widget for the initialisation of OpenGL
    virtual void initialize();

    /// resize function - called when the window size is changed
    virtual void resize(int width, int height);

    /// paint function - called when the window should be refreshed
    virtual void paint();

    /// keyboard interaction
    virtual void keyboard(int key, int scancode, int action, int mods);

    /// function that draws the planet system
    /// \param _projection the projection matrix for the scene
    /// \param _view the view matrix for the scene
    void draw_scene(mat4& _projection, mat4& _view);


private:

    /// default color shader (renders only texture)
    Shader   color_shader_;
    /// phong shader (renders texture and basic illumination)
    Shader   phong_shader_;

    /// simple shader for visualizing curves (just using solid color).
    Shader   solid_color_shader_;

    /// interval for the animation timer
    bool  timer_active_;
    /// update factor for the animation
    float time_step_;

    /// state whether the rendering should be in color or not
    bool greyscale_;

     /// rotation in x direction 
    float x_angle_;
    /// rotation in y direction 
    float y_angle_;

      /// eye's distance in radii from the observed point
    float dist_factor_;

    Path ship_path_renderer_, ship_path_cp_renderer_, circle1, circle2;
    Frame ship_path_frame_;

    // Control polygon for cubic spline
    std::vector<vec3> control_polygon_ = {
        {2.0, 2.0, 0.0},
        {3.0, 0.0, 0.0},
        {3.0, 0.0, -2.0},
        {0.0, 2.0, -3.0},
        {-2.0, 2.0, -1.0},
        {-3.0, 3.0, 0.0},
        {-1.0, 1.0, 1.0}, 
        {2.0, 2.0, 0.0}
    };

    /// the field of view for the camera
    float fovy_;
    /// the near plane for the virtual camera
    float near_;
    /// the far plane for the virtual camera
    float far_;


    /// current viewport dimension
    int  width_, height_;
};


//=============================================================================
#endif
//=============================================================================
