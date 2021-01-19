#include "tube_viewer.h"

#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"


#ifdef _WIN32
#  include <windows.h>
#  include <stdlib.h>
#  include <errhandlingapi.h>
#endif

//=============================================================================



int main(int argc, char *argv[])
{
#ifdef _WIN32
    // This make crashes very visible - without them, starting the
    // application from cmd.exe or powershell can surprisingly hide
    // any signs of a an application crash!
    SetErrorMode(0);
#endif
    IMGUI_CHECKVERSION();
    // TODO: Maximize window
    Tube_viewer window("Tube System", 1200, 640);
    return window.run();
}


