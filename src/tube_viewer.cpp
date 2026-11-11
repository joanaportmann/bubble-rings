//=============================================================================

#include "tube_viewer.h"
#include "glmath.h"
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <array>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <math.h>
#include <vector>
#include "glfw_window.h"
#include <iostream>
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"
#include <chrono>
#include "filament.h"

#include <Eigen/Dense>

#define numberOfVerticesPerTubeCircle 30
#define numEdges 20

using namespace std;

//=============================================================================

Tube_viewer::Tube_viewer(const char *_title, int _width, int _height)
    : GLFW_window(_title, _width, _height),
      unit_sphere_(50), // level of tessellation
      filament(0.12f, 4.0f, numEdges),
      tube(filament),
      background(0.0f, 0.0f, 50.0f, 0.0f)
{
    // rendering parameters
    fovy_ = 70;
    near_ = 0.01f;
    far_ = 20;

    x_angle_ = 0.0f;
    y_angle_ = -M_PI / 4.0f;
    dist_factor_ = 8.0f;

    originalFirstVertex = filament.getFilamentPoints()[0].position;

    srand((unsigned int)time(NULL));
}

//----------------------------------------------------------------------------

void Tube_viewer::keyboard(int key, int scancode, int action, int mods)
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
            y_angle_ += 0.05f * M_PI;
            break;
        }
        case GLFW_KEY_RIGHT:
        {
            y_angle_ -= 0.05f * M_PI;
            break;
        }

        case GLFW_KEY_R:
        {
            filament = Filament(thickness, circulation, numEdges);
            filament.setAlpha(alpha);
            filament.setResampleLength(length);
            filament.setModifyThickness(modifyThickness);
            filament.setRungeKutta4(rungeKutta4);
            filament.setEuler(euler);

            timer_active_ = false;
            break;
        }

        case GLFW_KEY_DOWN:
        {
            if (x_angle_ > -M_PI / 2 + 0.05 * M_PI)
                x_angle_ -= 0.05f * M_PI;
            break;
        }

        case GLFW_KEY_UP:
        {
            if (x_angle_ < M_PI / 2 - 0.05 * M_PI)
                x_angle_ += 0.05f * M_PI;
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

        case GLFW_KEY_8:
        {
            if (dist_factor_ >= 1.0f)
                dist_factor_ -= 0.2f;
            break;
        }

        case GLFW_KEY_9:
        {
            if (dist_factor_ <= 19.5f)
                dist_factor_ += 0.2f;
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

//--------------------------------------------------------------------------------

void Tube_viewer::initialize()
{
    tube.initialize();
    centerline_path_.initialize();

    // set initial state
    glClearColor(0.08f, 0.12f, 0.38f, 1.00f);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Allocate textures
    background.tex_.init(GL_TEXTURE0, GL_TEXTURE_2D, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR, GL_REPEAT);
    tube.tex_.init(GL_TEXTURE0, GL_TEXTURE_2D, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR, GL_REPEAT);

    // Load/generate textures
    background.tex_.loadPNG(TEXTURE_PATH "/underwaterSphere.png");
    tube.tex_.loadPNG(TEXTURE_PATH "/underwaterSphere.png");

    // setup shader
    background_shader_.load(SHADER_PATH "/background.vert", SHADER_PATH "/background.frag");   // Shader for background img
    bubble_ring_shader_.load(SHADER_PATH "/bubble_ring.vert", SHADER_PATH "/bubble_ring.frag"); // Shader for bubble ring (reflection, refraction, phong)
    solid_color_shader_.load(SHADER_PATH "/solid_color.vert", SHADER_PATH "/solid_color.frag"); // Shader for coordinate axis and centerline

    coordinateAxis.initialize();
}

//-----------------------------------------------------------------------------

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
    eye = center + rotation * vec4(0, 0, dist_factor_, 0);
    up = vec4(0, 1, 0, 0);

    mat4 view = MatUtils::look_at(eye, center, up);
    mat4 projection = MatUtils::perspective(fovy_, (float)width_ / (float)height_, near_, far_);
    vec3 eye_3d = vec3(eye(0), eye(1), eye(2));
    draw_scene(projection, view, eye_3d);

    bool show_demo_window = false;
    glClearColor(0.08f, 0.12f, 0.38f, 1.00f);

    // Start the Dear ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::StyleColorsClassic();

    ImGuiWindowFlags window_flags = ImGuiWindowFlags_AlwaysAutoResize;
    bool open_ptr = true;
    ImGui::Begin("Settings", &open_ptr, window_flags);

    if (ImGui::CollapsingHeader("Bubble ring parameters"))
    {
        ImGui::Text("Set initial configuration and click Reset to apply:");
        ImGui::SliderFloat("Thickness", &thickness, 0.01f, 0.8f);
        ImGui::SliderFloat("Circulation", &circulation, -4.0f, 50.0f);
    }
    if (ImGui::CollapsingHeader("Operations"))
    {
        if (ImGui::Checkbox("Modify thickness", &modifyThickness))
        {
            filament.setModifyThickness(modifyThickness);
        }
        ImGui::Text("Resampling and Catmull-Rom knot parameterization:");
        if (ImGui::SliderFloat("Resample length", &length, 0.02f, 1.0f))
        {
            filament.setResampleLength(length);
        }
        if (ImGui::SliderFloat("Alpha", &alpha, 0.0f, 1.0f))
        {
            filament.setAlpha(alpha);
        }
        if (ImGui::Checkbox("Runge Kutta 4", &rungeKutta4))
        {
            euler = !rungeKutta4;
            filament.setRungeKutta4(rungeKutta4);
            filament.setEuler(euler);
        }
        if (ImGui::Checkbox("Euler", &euler))
        {
            rungeKutta4 = !euler;
            filament.setEuler(euler);
            filament.setRungeKutta4(rungeKutta4);
        }
    }
    if (ImGui::CollapsingHeader("Visualization"))
    {
        ImGui::Checkbox("Recenter", &filament.recenter);
        ImGui::Checkbox("Render only Polygon", &renderOnlyPolygon);
        ImGui::Checkbox("Underwater background", &backgroundOn);
        ImGui::Checkbox("Show coordinate axis", &showCoordinateAxis);
    }

    if (ImGui::Button("Reset"))
    {
        filament = Filament(thickness, circulation, numEdges);
        filament.setAlpha(alpha);
        filament.setResampleLength(length);
        filament.setModifyThickness(modifyThickness);
        filament.setRungeKutta4(rungeKutta4);
        filament.setEuler(euler);

        timer_active_ = false;
    }
    if (ImGui::Button("Start"))
    {
        timer_active_ = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Stop"))
    {
        timer_active_ = false;
    }
    ImGui::SameLine();
    if (ImGui::Button("One step"))
    {
        filament.updateSkeleton();
    }

    std::string text = "Frame: %d";
    text += std::to_string(filament.framecouter);
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + ImGui::GetColumnWidth() - ImGui::CalcTextSize(text.c_str()).x - ImGui::GetScrollX() - 2 * ImGui::GetStyle().ItemSpacing.x);
    ImGui::Text("Frame: %d", filament.framecouter);
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

void Tube_viewer::draw_scene(mat4 &_projection, mat4 &_view, vec3 &eye)
{
    if (backgroundOn)
    {
        // render background
        Eigen::Affine3f scaling;
        scaling = Eigen::Scaling(background.radius_);
        mat4 m_matrix_bg = scaling.matrix().cast<double>();
        mat4 mv_matrix_bg = _view * m_matrix_bg;
        mat4 mvp_matrix_background = _projection * mv_matrix_bg;
        background_shader_.use();
        background_shader_.set_uniform("modelview_projection_matrix", mvp_matrix_background);
        background_shader_.set_uniform("tex", 0);
        background.tex_.bind();
        glCullFace(GL_FRONT);
        unit_sphere_.draw();
        glCullFace(GL_BACK);
    }

    mat4 m_matrix = mat4::Identity();
    mat4 mv_matrix = _view * m_matrix;
    mat4 mvp_matrix = _projection * mv_matrix;

    mat3 mv_matrix_3d;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mv_matrix_3d(i, j) = mv_matrix(i, j);
        }
    }
    mat3 normal_matrix = mv_matrix_3d.inverse().transpose();

    if (!renderOnlyPolygon)
    {
        Eigen::Matrix3d rotationSun1;
        rotationSun1 = Eigen::AngleAxisd(0.9425, vec3(0, 0, 1));

        Eigen::Matrix3d rotationSun2;
        rotationSun2 = Eigen::AngleAxisd(M_PI + 0.8378, vec3(0, 1, 0));

        vec3 light = normal_matrix * rotationSun2 * (rotationSun1 * vec3(0, 1, 0));

        bubble_ring_shader_.use();
        bubble_ring_shader_.set_uniform("modelview_projection_matrix", mvp_matrix);
        bubble_ring_shader_.set_uniform("normal_matrix", normal_matrix);
        bubble_ring_shader_.set_uniform("modelview_matrix", mv_matrix);
        bubble_ring_shader_.set_uniform("tex", 0);
        bubble_ring_shader_.set_uniform("eyePositionW", eye);
        bubble_ring_shader_.set_uniform("light_position", light);

        tube.draw();
    }
    else
    {
        std::vector<FilamentPoint> FilamentPoints = filament.getFilamentPoints();
        vector<vec3> controlPolygonForDebugging;
        controlPolygonForDebugging.reserve(FilamentPoints.size());
        for (size_t i = 0; i < FilamentPoints.size(); i++)
        {
            controlPolygonForDebugging.push_back(FilamentPoints[i].position);
        }

        vec3 offset_ = vec3(0, 0, 0);
        if (filament.recenter && !controlPolygonForDebugging.empty())
            offset_ = originalFirstVertex - controlPolygonForDebugging[0];
        centerline_path_.setPoints(controlPolygonForDebugging, offset_);

        solid_color_shader_.use();
        solid_color_shader_.set_uniform("modelview_projection_matrix", mvp_matrix);
        solid_color_shader_.set_uniform("color", vec4(1.0, 1.0, 1.0, 1.0));
        centerline_path_.draw();
    }

    if (showCoordinateAxis)
    {
        center_of_coordinatesystem = vec3(0, 0, 0);
        coordinateAxis.draw(solid_color_shader_, _projection * _view, center_of_coordinatesystem);
    }

    // check for OpenGL errors
    glCheckError();
}

//=============================================================================
