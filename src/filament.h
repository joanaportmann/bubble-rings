
#ifndef FILAMENT_H
#define FILAMENT_H
//=============================================================================

#include "gl.h"
#include <Eigen/Dense>
#include "glmath.h"
#include <vector>
//=============================================================================


struct FilamentPoint
{
    vec3 position;
    float a;
};


class Filament
{

public:

    // Constructor
    Filament(std::vector<FilamentPoint> tubeVertices_);

    std::vector<FilamentPoint> getFilamentPoints();

    std::vector<vec3> getBubbleRingSkeleton();

private:

    std::vector<FilamentPoint> controlPolygon_;

    // Todo
    void updateFilament();

    std::vector<vec3> circleVertices_t(int n, vec3 normal);

};


//=============================================================================
#endif
//=============================================================================
