//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================

#include "tube_viewer.h"
#include "glmath.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <array>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <math.h>

//=============================================================================



Tube_viewer::Tube_viewer(const char* _title, int _width, int _height)
	: GLFW_window(_title, _width, _height)
{
	// start animation
	timer_active_ = true;
	time_step_ = 1.0f / 24.0f; // one hour

	// rendering parameters
	greyscale_ = false;
	fovy_ = 45;
	near_ = 0.01f;
	far_ = 20;

	x_angle_ = -90.0f;
	y_angle_ = 0.0f;
	dist_factor_ = 9.0f;

	srand((unsigned int)time(NULL));
}

//-----------------------------------------------------------------------------

void
Tube_viewer::
keyboard(int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS || action == GLFW_REPEAT)
	{
		switch (key)
		{
		case GLFW_KEY_G:
		{
			greyscale_ = !greyscale_;
			break;
		}

		case GLFW_KEY_C:
			curve_display_mode_ = CurveDisplayMode((int(curve_display_mode_) + 1) % int(CURVE_SHOW_NUM_MODES));
			break;
		
		case GLFW_KEY_ESCAPE:
		{
			glfwSetWindowShouldClose(window_, GL_TRUE);
			break;
		}

		case GLFW_KEY_LEFT:
		{
			y_angle_ -= 10.0;
			break;
		}

		case GLFW_KEY_RIGHT:
		{
			y_angle_ += 10.0;
			break;
		}

		case GLFW_KEY_DOWN:
		{
			x_angle_ += 10.0;
			break;
		}

		case GLFW_KEY_UP:
		{
			x_angle_ -= 10.0;
			break;
		}

			// Key 9 increases and key 8 decreases the `dist_factor_` within the range - 2.5 < `dist_factor_` < 20.0.
		case GLFW_KEY_8:
		{
			if (dist_factor_ >= 3.0) dist_factor_ -= 0.5;
			break;
		}

		case GLFW_KEY_9:
		{
			if (dist_factor_ <= 19.5) dist_factor_ += 0.5;
			break;
		}
		}
	}
}

//-----------------------------------------------------------------------------


void Tube_viewer::resize(int _width, int _height)
{
	width_ = _width;
	height_ = _height;
	glViewport(0, 0, _width, _height);
}

//-----------------------------------------------------------------------------

std::vector<vec3> unitCircleVertices(int n) {
	std::vector<vec3> vertices;

	for (int i=0; i < n; i++) {
	
		float x = cos(2 * M_PI / n * i);
		float y = sin(2 * M_PI / n * i);
		vec4 vertex =  mat4::translate(vec3(2, 2, 0)) * vec4(x, y, 0, 1);
		vec3 vertex_3 = vec3(vertex.x, vertex.y, vertex.z);
		vertices.push_back(vertex_3);
	}

	return vertices;
}

//--------------------------------------------------------------------------------


void Tube_viewer::initialize()
{
	// set initial state
	glClearColor(0, 0, 0, 0);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// setup shaders
	color_shader_.load(SHADER_PATH "/color.vert", SHADER_PATH "/color.frag");
	phong_shader_.load(SHADER_PATH "/phong.vert", SHADER_PATH "/phong.frag");

	solid_color_shader_.load(SHADER_PATH "/solid_color.vert", SHADER_PATH "/solid_color.frag");

	ship_path_renderer_.initialize();
	ship_path_cp_renderer_.initialize();
	ship_path_frame_.initialize();
	unitCircle.initialize();

	ship_path_.set_control_polygon(control_polygon_, true);
	ship_path_renderer_.sample(ship_path_);
	ship_path_cp_renderer_.setPoints(ship_path_.bezier_control_points());

	unitCircle.setPoints(unitCircleVertices(8));
}
//-----------------------------------------------------------------------------


void Tube_viewer::paint()
{
	// clear framebuffer and depth buffer first
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	vec4 eye, center, up;
	float radius, x_rotation, y_rotation;
	mat4 rotation;

	center = vec4(0, 0, 0, 0);

	x_rotation = x_angle_;
	y_rotation = y_angle_;
	rotation = mat4::rotate_y(y_rotation) * mat4::rotate_x(x_rotation);
	eye = center + rotation * vec4(0, 0, -(dist_factor_ ), 0);
	up = rotation * vec4(0, 1, 0, 0);

	
	mat4    view = mat4::look_at(vec3(eye), (vec3)center, (vec3)up);

	mat4 projection = mat4::perspective(fovy_, (float)width_ / (float)height_, near_, far_);
	draw_scene(projection, view);
	}


//-----------------------------------------------------------------------------

		
void Tube_viewer::draw_scene(mat4& _projection, mat4& _view)
{
		ship_path_frame_.draw(solid_color_shader_, _projection * _view, ship_path_(ship_path_param_));
		
		solid_color_shader_.use();
		solid_color_shader_.set_uniform("modelview_projection_matrix", _projection * _view);
		solid_color_shader_.set_uniform("color", vec4(0, 0.8, 0.8, 1.0));
		ship_path_cp_renderer_.draw();
		unitCircle.draw();

		// Bezier curve
		solid_color_shader_.use();
		solid_color_shader_.set_uniform("modelview_projection_matrix", _projection * _view);
		solid_color_shader_.set_uniform("color", vec4(1.0, 0.0, 0.0, 1.0));
		ship_path_renderer_.draw();


	// the sun is centered at the origin and -- for lighting -- considered to be a point, so that is the light position in world coordinates
	vec4 light = vec4(0.0, 0.0, 0.0, 1.0); //in world coordinates
	// convert light into camera coordinates
	light = _view * light;

	// check for OpenGL errors
	glCheckError();
}

//=============================================================================
