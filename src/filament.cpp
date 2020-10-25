//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================

#include "filament.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdlib.h> 

using namespace std;

//=============================================================================


Filament::Filament(std::vector<FilamentPoint> CP)
    : controlPolygon_(CP)
{
}

std::vector<FilamentPoint> Filament::getFilamentPoints() {
    return controlPolygon_;
}

//-----------------------------------------------------------------------------