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
#include <stdlib.h> /* srand, rand */
#include <time.h>	/* time */
#include <array>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <math.h>

#include <Eigen/Dense>
//#include <iostream>

//#include <Eigen/src/Geometry/Quaternion.h>
//using Eigen::Vector3d;
using namespace std;

//=============================================================================

Tube_viewer::Tube_viewer(const char* _title, int _width, int _height)
	: GLFW_window(_title, _width, _height),
	tube(control_polygon_)
{
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

void Tube_viewer::
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

		case GLFW_KEY_ESCAPE:
		{
			glfwSetWindowShouldClose(window_, GL_TRUE);
			break;
		}

		case GLFW_KEY_LEFT:
		{
			y_angle_ -= 0.2 * M_PI;
			break;
		}

		case GLFW_KEY_RIGHT:
		{
			y_angle_ += 0.1 * M_PI;
			break;
		}

		case GLFW_KEY_DOWN:
		{
			x_angle_ += 0.1 * M_PI;
			break;
		}

		case GLFW_KEY_UP:
		{
			x_angle_ -= 0.1 * M_PI;
			break;
		}

			// Key 9 increases and key 8 decreases the `dist_factor_` within the range - 2.5 < `dist_factor_` < 20.0.
		case GLFW_KEY_8:
		{
			if (dist_factor_ >= 3.0)
				dist_factor_ -= 0.5;
			break;
		}

		case GLFW_KEY_9:
		{
			if (dist_factor_ <= 19.5)
				dist_factor_ += 0.5;
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

std::vector<vec3> circleVertices(int n, vec3 center, vec3 normal, float radius)
{
	std::vector<vec3> vertices;

	for (int i = 0; i < n; i++)
	{
		float x = cos(2 * M_PI / n * i) * radius;
		float y = sin(2 * M_PI / n * i) * radius;

		auto rotation = Eigen::Quaterniond::FromTwoVectors(vec3(0, 0, 1), normal);

		vec3 vertex = (rotation.matrix() * vec3(x, y, 0)) + center;
		vertices.push_back(vertex);
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
	//tube_shader_.load(SHADER_PATH   "/tube.vert", SHADER_PATH   "/tube.frag");

	solid_color_shader_.load(SHADER_PATH "/solid_color.vert", SHADER_PATH "/solid_color.frag");

	ship_path_renderer_.initialize();
	ship_path_cp_renderer_.initialize();
	ship_path_frame_.initialize();

	//ship_path_.set_control_polygon(control_polygon_, true);
	//ship_path_renderer_.sample(ship_path_);
	ship_path_cp_renderer_.setPoints(control_polygon_);
}
//-----------------------------------------------------------------------------

void Tube_viewer::drawCircle(std::vector<vec3> control_polygon_, float radius)
{
	//std::vector<Path> circles;
	for (int i = 0; i < control_polygon_.size(); i++)
	{
		Path circle;
		circle.initialize();
		std::vector<vec3> verticesOfOneCircle = circleVertices(
			7,
			control_polygon_[i],
			control_polygon_[(i + 1) % control_polygon_.size()] - control_polygon_[i],
			radius
		);
		circle.setPoints(verticesOfOneCircle);
		circle.draw();
	};
};

void Tube_viewer::paint()
{
	// clear framebuffer and depth buffer first
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	vec4 eye, center, up;
	float x_rotation, y_rotation;
	mat4 rotation;

	center = vec4(0, 0, 0, 0);

	x_rotation = x_angle_;
	y_rotation = y_angle_;

	Eigen::Affine3f rotation_x;
	rotation_x = Eigen::AngleAxisf(x_rotation, vec3::UnitX().cast<float>());
	mat4 rotation_x_matrix = rotation_x.matrix().cast<double>();

	Eigen::Affine3f rotation_y;
	rotation_y = Eigen::AngleAxisf(y_rotation, vec3::UnitY().cast<float>());
	mat4 rotation_y_matrix = rotation_y.matrix().cast<double>();

	rotation = rotation_y_matrix * rotation_x_matrix;
	eye = center + rotation * vec4(0, 0, -(dist_factor_), 0);
	up = rotation * vec4(0, 1, 0, 0);

	mat4 view;
	view = MatUtils::look_at(
		eye,
		center,
		up);
	mat4 projection;
	projection = MatUtils::perspective(fovy_, (float)width_ / (float)height_, near_, far_);
	draw_scene(projection, view);
}

//-----------------------------------------------------------------------------

void Tube_viewer::draw_scene(mat4 &_projection, mat4 &_view)
{
	// the matrices we need: model, modelview, modelview-projection, normal
	mat4 m_matrix;
	mat4 mv_matrix;
	mat4 mvp_matrix;
	mat3 n_matrix;

	// the sun is centered at the origin and -- for lighting -- considered to be a point, so that is the light position in world coordinates
	vec4 light = vec4(0.0, 0.0, 0.0, 1.0); //in world coordinates
	// convert light into camera coordinates
	light = _view * light;

	// render polygonpath
	mat4 matrix;
	matrix = _projection * _view;
	mv_matrix = matrix;
	mvp_matrix = matrix;

	solid_color_shader_.use();
	solid_color_shader_.set_uniform("modelview_projection_matrix", matrix);
	solid_color_shader_.set_uniform("color", vec4(0.8, 0.8, 0.2, 0.6));
	ship_path_cp_renderer_.draw();
	//tube_shader_.use();
	//tube_shader_.set_uniform("modelview_projection_matrix", mvp_matrix);
	//tube_shader_.set_uniform("tex", 0);
	//tube_shader_.set_uniform("greyscale", (int)greyscale_);
	//sun_.tex_.bind();






//   Eigen::Affine3d r = create_rotation_matrix(1.0, 1.0, 1.0);
//   Eigen::Affine3d t(Eigen::Translation3d(Eigen::Vector3d(1,1,2)));

//   Eigen::Matrix4d m = (t * r).matrix(); // Option 1

//   //Eigen::Matrix4d m = t.matrix(); // Option 2
//   m *= r.matrix();

	// n_matrix << 1, 2, 3,
    //  4, 5, 6,
    //  7, 8, 9;


	// phong_shader_.use();
	// phong_shader_.set_uniform("modelview_projection_matrix", mvp_matrix);
	// phong_shader_.set_uniform("modelview_matrix", mv_matrix);
	// phong_shader_.set_uniform("normal_matrix", n_matrix);
	// //phong_shader_.set_uniform("tex", 0);
	// phong_shader_.set_uniform("greyscale", (int)greyscale_);
	tube.draw();
	
	// render circles around polygonpath
	drawCircle(control_polygon_, 0.3);
	// check for OpenGL errors
	glCheckError();
}

// 	Eigen::Affine3d create_rotation_matrix(double ax, double ay, double az) {
//   Eigen::Affine3d rx =
//       Eigen::Affine3d(Eigen::AngleAxisd(ax, Eigen::Vector3d(1, 0, 0)));
//   Eigen::Affine3d ry =
//       Eigen::Affine3d(Eigen::AngleAxisd(ay, Eigen::Vector3d(0, 1, 0)));
//   Eigen::Affine3d rz =
//       Eigen::Affine3d(Eigen::AngleAxisd(az, Eigen::Vector3d(0, 0, 1)));
//   return rz * ry * rx;
// }
//=============================================================================
