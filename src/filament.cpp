#include "filament.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdlib.h>

using namespace std;

#define numberOfVerticesPerTubeCircle 30
#define _USE_MATH_DEFINES
#define RM_mu 0.4723665527
#define delta 0.6420127083
#define circulation 4
#define g -9.8
#define kinematic_viscosity 1e-06
#define At -1

//=============================================================================

Filament::Filament(std::vector<FilamentPoint> CP)
    : controlPolygon_(CP)
{
}

std::vector<FilamentPoint> Filament::getFilamentPoints()
{
    return controlPolygon_;
}

//-----------------------------------------------------------------------------

std::vector<vec3> verticesofOneCircle_(int n, vec3 center, vec3 normal, float radius)
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

//----------------------------------------------------------------------------------

std::vector<vec3> Filament::getBubbleRingSkeleton()
{
    std::vector<vec3> verticesOfTube;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 edgeAfter = controlPolygon_[(i + 1) % controlPolygon_.size()].position - controlPolygon_[i].position;
        vec3 edgeBefore = controlPolygon_[i].position - controlPolygon_[(i - 1 + controlPolygon_.size()) % controlPolygon_.size()].position;
        std::vector<vec3> verticesOfOneCircle = verticesofOneCircle_(
            numberOfVerticesPerTubeCircle,
            controlPolygon_[i].position,
            (edgeBefore + edgeAfter).normalized(),
            controlPolygon_[i].a);
        for (int j = 0; j < verticesOfOneCircle.size(); j++)
        {
            verticesOfTube.push_back(verticesOfOneCircle[j]);
        }
    };

    return verticesOfTube;
};

//------------------------------------------------------------------------------------

//Biotsavart velocity

vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a)
{
    float aSqr = a * a * RM_mu * RM_mu;
    vec3 R0_ = R0 - p;
    vec3 R1_ = R1 - p;
    vec3 RPrime = R1 - R0;
    vec3 cross01 = R0_.cross(R1_);
    float r1 = R1_.dot(RPrime) / (sqrt(aSqr + R1_.norm()) * (RPrime.norm() * aSqr + cross01.norm()));
    float r0 = R0_.dot(RPrime) / (sqrt(aSqr + R0_.norm()) * (RPrime.norm() * aSqr + cross01.norm()));

    return Gamma * (r1 - r0) * cross01 / (4 * M_PI);
}

//-------------------------------------------------------------------------------------

vec3 Filament::biotSavart()
{
    vec3 temp_vel = vec3(0, 0, 0);

    for (int j = 1; j < controlPolygon_.size(); j++)
    {
        vec3 R0 = controlPolygon_[j - 1].position - controlPolygon_[j].position;
        vec3 R1 = controlPolygon_[j + 1].position - controlPolygon_[j].position;
        float a = controlPolygon_[j].a;
        float Gamma = circulation;
        vec3 position = controlPolygon_[j].position;

        temp_vel += vec3(1, 1, 1)
        //biotsavartedge(position, R0, R1, Gamma, a);
    }

    return temp_vel;
}

//-------------------------------------------------------------------------------------

// void Filament::updateSkeleton()
// {
//     vec3 velocity;
//     velocity = biotSavart();
// };

//-----------------------------------------------------------------------------------