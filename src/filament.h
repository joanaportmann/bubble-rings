//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================
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

private:

    std::vector<FilamentPoint> controlPolygon_;

    // Todo
    void updateFilament();

};


//=============================================================================
#endif
//=============================================================================
