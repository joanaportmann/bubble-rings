
//=============================================================================
#include "glfw_window.h"
#include <iostream>
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"

// About Desktop OpenGL function loaders:
//  Modern desktop OpenGL doesn't have a standard portable header file to load OpenGL function pointers.
//  Helper libraries are often used for this purpose! Here we are supporting a few common ones (gl3w, glew, glad).
//  You may use another loader/header of your choice (glext, glLoadGen, etc.), or chose to manually implement your own.
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
#include <GL/gl3w.h> // Initialize with gl3wInit()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
#include <GL/glew.h> // Initialize with glewInit()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
#include <glad/glad.h> // Initialize with gladLoadGL()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD2)
#include <glad/gl.h> // Initialize with gladLoadGL(...) or gladLoaderLoadGL()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING2)
#define GLFW_INCLUDE_NONE      // GLFW including OpenGL headers causes ambiguity or multiple definition errors.
#include <glbinding/Binding.h> // Initialize with glbinding::Binding::initialize()
#include <glbinding/gl/gl.h>
using namespace gl;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING3)
#define GLFW_INCLUDE_NONE        // GLFW including OpenGL headers causes ambiguity or multiple definition errors.
#include <glbinding/glbinding.h> // Initialize with glbinding::initialize()
#include <glbinding/gl/gl.h>
using namespace gl;
#else
#include IMGUI_IMPL_OPENGL_LOADER_CUSTOM
#endif

// Include glfw3.h after our OpenGL definitions
#include <GLFW/glfw3.h>

// [Win32] Our example includes a copy of glfw3.lib pre-compiled with VS2010 to maximize ease of testing and compatibility with old VS compilers.
// To link with VS2010-era libraries, VS2015+ requires linking with legacy_stdio_definitions.lib, which we do using this pragma.
// Your own project should not be affected, as you are likely to link with a newer binary of GLFW that is adequate for your version of Visual Studio.
#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

//=============================================================================

GLFW_window *GLFW_window::instance__ = NULL;

//-----------------------------------------------------------------------------

static void glfw_error_callback(int error, const char *description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

//-----------------------------------------------------------------------------

GLFW_window::GLFW_window(const char *_title, int _width, int _height)
{
    // initialize glfw window
    if (!glfwInit())
    {
        std::cerr << "Cannot initialize GLFW!\n";
        exit(EXIT_FAILURE);
    }

    // request core profile and OpenGL version 3.2
    const char *glsl_version = "#version 460";
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

    // try to create window
    window_ = glfwCreateWindow(_width, _height, _title, NULL, NULL);
    if (!window_)
    {
        std::cerr << "Window creation failed!\n";
        std::cerr << "Attempting fall-back (for INF03 machines)" << std::endl;

        // Request OpenGL version 3.1
        // Note: below version 3.2, we must request GLFW_OPENGL_ANY_PROFILE
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_ANY_PROFILE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
        window_ = glfwCreateWindow(_width, _height, _title, NULL, NULL);

        if (!window_)
        {
            glfwTerminate();
            std::cerr << "Window creation failed!\n";
            exit(EXIT_FAILURE);
        }
    }

    // make this window the current one
    glfwMakeContextCurrent(window_);

    // enable vsync
    glfwSwapInterval(1);

    // register glfw callbacks
    glfwSetKeyCallback(window_, keyboard__);
    glfwSetFramebufferSizeCallback(window_, resize__);
    glfwSetErrorCallback(glfw_error_callback);

    // Initialize OpenGL loader
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
    bool err = gl3wInit() != 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
    bool err = glewInit() != GLEW_OK;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
    bool err = gladLoadGL() == 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD2)
    bool err = gladLoadGL(glfwGetProcAddress) == 0; // glad2 recommend using the windowing library loader instead of the (optionally) bundled one.
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING2)
    bool err = false;
    glbinding::Binding::initialize();
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING3)
    bool err = false;
    glbinding::initialize([](const char *name) { return (glbinding::ProcAddress)glfwGetProcAddress(name); });
#else
    bool err = false; // If you use IMGUI_IMPL_OPENGL_LOADER_CUSTOM, your loader is likely to requires some form of initialization.
#endif
    if (err)
    {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        exit(EXIT_FAILURE);
    }

    // Setup Dear ImGui context

    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    // Setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(window_, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // now that we have a GL context, initialize GLEW
    glewExperimental = GL_TRUE;
    GLenum err1 = glewInit();
    if (err1 != GLEW_OK)
    {
        std::cerr << "Error initializing GLEW: " << glewGetErrorString(err1) << std::endl;
        exit(1);
    }

    // debug: print GL and GLSL version
    std::cout << "GLEW   " << glewGetString(GLEW_VERSION) << std::endl;
    std::cout << "GL     " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GLSL   " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;

    // call glGetError once to clear error queue
    GLenum error = glGetError();

    instance__ = this;
}

//-----------------------------------------------------------------------------

GLFW_window::~GLFW_window()
{
    glfwTerminate();
}

//-----------------------------------------------------------------------------

int GLFW_window::run()
{

    // initialize OpenGL
    initialize();

    // query framebuffer width and height
    // call resize to initialize viewport
    int width, height;
    glfwGetFramebufferSize(window_, &width, &height);
    resize(width, height);

    // now run the event loop
    while (!glfwWindowShouldClose(window_))
    {
        // call timer function
        timer();

        // draw scene
        paint();

        // swap buffers
        // glfwSwapBuffers(window_);

        // handle events
        // Poll and handle events (inputs, window resize, etc.)
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        glfwPollEvents();

        // Our state
        bool show_demo_window = true;
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // render your GUI
        ImGui::Begin("Demo window");
        ImGui::Button("Hello!");
        ImGui::End();

        // Render dear imgui into screen
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        int display_w, display_h;
        glfwGetFramebufferSize(window_, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glfwSwapBuffers(window_);
    }

    glfwDestroyWindow(window_);
    glClear(GL_COLOR_BUFFER_BIT);

    return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------

void GLFW_window::error__(int error, const char *description)
{
    fputs(description, stderr);
}

void GLFW_window::keyboard__(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    instance__->keyboard(key, scancode, action, mods);
}

void GLFW_window::resize__(GLFWwindow *window, int width, int height)
{
    instance__->resize(width, height);
    instance__->paint();
    glfwSwapBuffers(instance__->window_);
}

//=============================================================================
