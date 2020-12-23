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
#include "glfw_window.h"
#include <iostream>
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"

#include <Eigen/Dense>

#define numberOfVerticesPerTubeCircle 30

//using Eigen::Vector3d;
using namespace std;

//=============================================================================

Tube_viewer::Tube_viewer(const char *_title, int _width, int _height)
	: GLFW_window(_title, _width, _height),
	  filament(0.12, 4),
	  tube(filament)
{
	// rendering parameters
	greyscale_ = false;
	fovy_ = 70;
	near_ = 0.01f;
	far_ = 20;

	x_angle_ = 0.0f;
	y_angle_ = 0.0f;
	dist_factor_ = 9.0f;

	srand((unsigned int)time(NULL));
}

//----------------------------------------------------------------------------

void Tube_viewer::
	keyboard(int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS || action == GLFW_REPEAT)
	{
		switch (key)
		{

		case GLFW_KEY_ESCAPE:
		{
			glfwSetWindowShouldClose(window_, GL_TRUE);
			break;
		}

		case GLFW_KEY_LEFT:
		{
			y_angle_ -= 0.05 * M_PI;
			break;
		}

		case GLFW_KEY_R:
		{
			filament = Filament(thickness, circulation);
			Tube tube_(filament);
			Tube *tube;
			tube = &tube_;

			break;
		}

		case GLFW_KEY_RIGHT:
		{
			y_angle_ += 0.05 * M_PI;
			break;
		}

		case GLFW_KEY_DOWN:
		{
			x_angle_ += 0.05 * M_PI;
			break;
		}

		case GLFW_KEY_UP:
		{
			x_angle_ -= 0.05 * M_PI;
			break;
		}

		case GLFW_KEY_SPACE:
		{
			timer_active_ = !timer_active_;
			break;
		}

		case GLFW_KEY_S:
		{

			filament.updateSkeleton();
			break;
		}

			// Key 9 increases and key 8 decreases the `dist_factor_` within the range - 2.5 < `dist_factor_` < 20.0.
		case GLFW_KEY_8:
		{
			if (dist_factor_ >= 1.0)
				dist_factor_ -= 0.2;
			break;
		}

		case GLFW_KEY_9:
		{
			if (dist_factor_ <= 19.5)
				dist_factor_ += 0.2;
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
	tube.initialize();
	// set initial state
	glClearColor(0, 0, 0, 0);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// setup shader
	color_shader_.load(SHADER_PATH "/color.vert", SHADER_PATH "/color.frag");
	phong_shader_.load(SHADER_PATH "/phong.vert", SHADER_PATH "/phong.frag");

	solid_color_shader_.load(SHADER_PATH "/solid_color.vert", SHADER_PATH "/solid_color.frag");
	test_tube_shader_.load(SHADER_PATH "/test_tube.vert", SHADER_PATH "/test_tube.frag");

	ship_path_renderer_.initialize();
	ship_path_cp_renderer_.initialize();
	ship_path_frame_.initialize();

	//ship_path_cp_renderer_.setPoints(control_polygon_);
}
//-----------------------------------------------------------------------------

// void Tube_viewer::drawCircle(std::vector<vec3> control_polygon_, float radius)
// {
// 	//std::vector<Path> circles;
// 	for (int i = 0; i < control_polygon_.size(); i++)
// 	{
// 		Path circle;
// 		circle.initialize();
// 		vec3 edgeAfter = control_polygon_[(i + 1) % control_polygon_.size()] - control_polygon_[i];
// 		vec3 edgeBefore = control_polygon_[i] - control_polygon_[(i - 1 + control_polygon_.size()) % control_polygon_.size()];
// 		std::vector<vec3> verticesOfOneCircle = circleVertices(
// 			numberOfVerticesPerTubeCircle,
// 			control_polygon_[i],
// 			(edgeAfter + edgeBefore).normalized(),
// 			radius
// 		);
// 		circle.setPoints(verticesOfOneCircle);
// 		//circle.draw();
// 	};
// };

void Tube_viewer::paint()
{
	// clear framebuffer and depth buffer first
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	vec4 eye, center, up;
	float x_rotation, y_rotation;
	mat4 rotation;

	center = vec4(0, 0, 0, 0);

	x_rotation = -x_angle_;
	y_rotation = -y_angle_;

	Eigen::Affine3f rotation_x;
	rotation_x = Eigen::AngleAxisf(x_rotation, vec3::UnitX().cast<float>());
	mat4 rotation_x_matrix = rotation_x.matrix().cast<double>();

	Eigen::Affine3f rotation_y;
	rotation_y = Eigen::AngleAxisf(y_rotation, vec3::UnitY().cast<float>());
	mat4 rotation_y_matrix = rotation_y.matrix().cast<double>();

	rotation = rotation_y_matrix * rotation_x_matrix;
	eye = center + rotation * vec4((dist_factor_), (dist_factor_), (dist_factor_), 0);
	up = rotation * vec4(0, 1, 0, 0);

	mat4 view;
	view = MatUtils::look_at(
		eye,
		center,
		up);
	mat4 projection;
	projection = MatUtils::perspective(fovy_, (float)width_ / (float)height_, near_, far_);
	draw_scene(projection, view);

	// Our state
	bool show_demo_window = false;
	glClearColor(0.08f, 0.12f, 0.38f, 1.00f);

	// Start the Dear ImGui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	// render GUI
	ImGui::Begin("Settings");
	ImGui::Text("Set start configuration of bubble ring.");
	ImGui::SliderFloat("Thickness", &thickness, 0.0f, 0.5f);
	ImGui::SliderFloat("Circulation", &circulation, 0.0f, 10.0f);
	if (ImGui::Button("Reset"))
	{
		
		filament = Filament(thickness, circulation);
		Tube tube_(filament);
		Tube *tube;
		tube = &tube_;
		timer_active_ = false;
	} 
	if (ImGui::Button("Start"))
	{
		timer_active_ = true;
	} ImGui::SameLine();
	if (ImGui::Button("Stop"))
	{
		timer_active_ = false;
	} ImGui::SameLine();
		if (ImGui::Button("One step"))
	{
			filament.updateSkeleton();
	}
	ImGui::End();

	if (show_demo_window)
		ImGui::ShowDemoWindow(&show_demo_window);

	// Render dear imgui into screen
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	int display_w, display_h;
	glfwGetFramebufferSize(window_, &display_w, &display_h);
	glViewport(0, 0, display_w, display_h);
	glfwSwapBuffers(window_);
}
//-----------------------------------------------------------------------------

void Tube_viewer::timer()
{
	if (timer_active_)
		filament.updateSkeleton();
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

	Eigen::Affine3f rotation_y;
	rotation_y = Eigen::AngleAxisf(0.0f, vec3::UnitY().cast<float>());
	mat4 rotation_y_mat = rotation_y.matrix().cast<double>();

	// render tube
	m_matrix = mat4::Identity();
	mv_matrix = _view * m_matrix;
	mvp_matrix = _projection * mv_matrix;
	mat3 normal_matrix;
	normal_matrix = mat3::Identity();

	test_tube_shader_.use();
	test_tube_shader_.set_uniform("modelview_projection_matrix", mvp_matrix);
	test_tube_shader_.set_uniform("normal_matrix", normal_matrix);
	//test_tube_shader_.set_uniform("light_position", _view * vec4(0, 0, 0, 1));
	// test_tube_shader_.set_uniform("color", vec4(0.8, 0.8, 0.2, 0.6));
	tube.draw();
	center_of_coordinatesystem = vec3(0, 0, 0);

	ship_path_frame_.draw(solid_color_shader_, _projection * _view, center_of_coordinatesystem);

	// render circles around polygonpath
	//drawCircle(control_polygon_, 0.3);
	// check for OpenGL errors
	glCheckError();
}

//=============================================================================
