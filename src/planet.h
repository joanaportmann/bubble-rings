//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================
#ifndef PLANET_H
#define PLANET_H
//=============================================================================

#include "texture.h"
#include "glmath.h"

//=============================================================================


/// default class for the planet
class Planet
{
public:

    /// default constructor
    Planet(const float _angle_step_orbit,
           const float _angle_step_self,
           const float _radius,
           const float _distance) :
        pos_(_distance,0,0,1),
        distance_(_distance),
        radius_(_radius),
        angle_orbit_(0),
        angle_self_(0),
        angle_step_orbit_(_angle_step_orbit),
        angle_step_self_(_angle_step_self)
    {}

    
    /// the current position of the planet
    vec4 pos_;
    /// distance of the planet from the orbit center (=origin)
    float distance_;
    /// the radius of the planet
    float radius_;

    /// current rotation around the orbit center
    float angle_orbit_;
    /// current rotation around planet center
    float angle_self_;

    /// how much to sweep around orbit in one day (in radians)
    float angle_step_orbit_;
    /// how much to rotation around planet center in one day (in radians)
    float angle_step_self_;


    /// main diffuse texture for the planet
    Texture tex_;
};

//=============================================================================
#endif // PLANET_H
//=============================================================================
