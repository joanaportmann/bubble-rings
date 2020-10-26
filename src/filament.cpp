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

//=============================================================================


Filament::Filament(std::vector<FilamentPoint> CP)
    : controlPolygon_(CP)
{
}

std::vector<FilamentPoint> Filament::getFilamentPoints() {
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
            controlPolygon_[i].a
        );
        for (int j = 0; j < verticesOfOneCircle.size(); j++)
        {
            verticesOfTube.push_back(verticesOfOneCircle[j]);
        }
    };

    return verticesOfTube;
};

//------------------------------------------------------------------------------------