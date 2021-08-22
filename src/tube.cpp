
//=============================================================================

#include "tube.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "filament.h"
#include <chrono>

// ATTENTION: Keep in sync with the one in filament.cpp
#define numberOfVerticesPerTubeCircle 30

using namespace std;

//=============================================================================

Tube::Tube(Filament &filament)
    : filament_(filament)
{
}

//-----------------------------------------------------------------------------

Tube::~Tube()
{
    if (vbo_)
        glDeleteBuffers(1, &vbo_);
    if (nbo_)
        glDeleteBuffers(1, &nbo_);
    if (tan_bo_)
        glDeleteBuffers(1, &tan_bo_);
    if (bitan_bo_)
        glDeleteBuffers(1, &bitan_bo_);
    if (tbo_)
        glDeleteBuffers(1, &tbo_);
    if (ibo_)
        glDeleteBuffers(1, &ibo_);
    if (vao_)
        glDeleteVertexArrays(1, &vao_);
}

//----------------------------------------------------------------------------

void Tube::createTriangleStruct()
{
    std::vector<FilamentPoint> filamentPoints = filament_.getFilamentPoints();

    tubeVertices.clear();
    tubeVertices = filament_.getBubbleRingSkeleton();

    for (unsigned int v = 0; v < tubeVertices.size(); ++v)
    {
        Triangle triangle1, triangle2;

        //Calculate closest vertices of next circle to get better visualization (without twist)
        // if (v % numberOfVerticesPerTubeCircle == 1)
        // {
        //     shift = 0;
        //     double dist = (tubeVertices[v] - tubeVertices[v + numberOfVerticesPerTubeCircle]).norm();
        //     for (int u = v + numberOfVerticesPerTubeCircle + 1; u < v + (2 * numberOfVerticesPerTubeCircle); u++)
        //     {
        //         double current_dist = (tubeVertices[v] - tubeVertices[u]).norm();
        //         if (current_dist < dist)
        //         {
        //             dist = current_dist;
        //             shift = u % numberOfVerticesPerTubeCircle;
        //         }
        //     }
        // }

        bool lastVertexInCircle = (v % numberOfVerticesPerTubeCircle == numberOfVerticesPerTubeCircle - 1);

        unsigned int i0 = v;
        unsigned int i1 = (lastVertexInCircle ? (v - (numberOfVerticesPerTubeCircle - 1)) : (v + 1)) % tubeVertices.size();
        unsigned int i2 = (v + numberOfVerticesPerTubeCircle) % tubeVertices.size();
        unsigned int i3 = (lastVertexInCircle ? (v + 1) : v + (numberOfVerticesPerTubeCircle + 1)) % tubeVertices.size();

        triangle1.ind0 = i0;
        triangle1.ind1 = i2;
        triangle1.ind2 = i3;
        Tube::triangles_.push_back(triangle1);

        triangle2.ind0 = i0;
        triangle2.ind1 = i3;
        triangle2.ind2 = i1;
        Tube::triangles_.push_back(triangle2);
    }
}

//------------------------------------------------------------------------------------

/** * 
 * Computes triangle's normals 
 */

void Tube::compute_normals()
{
    int c = 0;
    // compute triangle normals
    for (Triangle &t : triangles_)
    {
        const vec3 &p0 = tubeVertices[t.ind0];
        const vec3 &p1 = tubeVertices[t.ind1];
        const vec3 &p2 = tubeVertices[t.ind2];

        t.normal = ((p1 - p0).cross(p2 - p0)).normalized();
    }
}

//--------------------------------------------------------------------------------------

void Tube::initialize()
{
    // generate vertex array object
    glGenVertexArrays(1, &vao_);

    // generate buffers
    glGenBuffers(1, &vbo_);
    glGenBuffers(1, &nbo_);
}

//------------------------------------------------------------------------------------

int Tube::wrap2(int i)
{
    return (i + 2 * tubeVertices.size()) % (2 * tubeVertices.size());
}

//---------------------------------------------------------------------------------

void Tube::updateBuffers()
{
    triangles_.clear();
    Tube::createTriangleStruct();

    // To measure time (debugging)
    // auto start = std::chrono::steady_clock::now();
    // ....
    // auto end = std::chrono::steady_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end - start;
    // std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    Tube::compute_normals();
    std::vector<GLfloat> positions(3 * 3 * triangles_.size());
    std::vector<GLfloat> faceNormals(3 * triangles_.size());
    std::vector<GLfloat> normals(3 * 3 * triangles_.size());

    unsigned int p(0), i(0), t(0);

    // generate triangles
    int count = 0;
    for (Triangle &t : triangles_)
    {
        positions[p++] = tubeVertices[t.ind0](0);
        positions[p++] = tubeVertices[t.ind0](1);
        positions[p++] = tubeVertices[t.ind0](2);

        positions[p++] = tubeVertices[t.ind1](0);
        positions[p++] = tubeVertices[t.ind1](1);
        positions[p++] = tubeVertices[t.ind1](2);

        positions[p++] = tubeVertices[t.ind2](0);
        positions[p++] = tubeVertices[t.ind2](1);
        positions[p++] = tubeVertices[t.ind2](2);
    }

    for (int i = 0; i < tubeVertices.size(); i++)
    {
        vec3 normal0 = triangles_[wrap2(2 * i)].normal;
        vec3 normal1 = triangles_[wrap2(2 * i + 1)].normal;
        vec3 normal2 = triangles_[wrap2(2 * i - numberOfVerticesPerTubeCircle)].normal;
        vec3 normal3 = triangles_[wrap2(2 * i - numberOfVerticesPerTubeCircle - 1)].normal;
        vec3 normal4 = triangles_[wrap2(2 * i - numberOfVerticesPerTubeCircle - 2)].normal;
        vec3 normal5 = triangles_[wrap2(2 * i - 1)].normal;

        vec3 normal = 0.25 * (normal2 - normal5) + 0.125 * (normal0 - normal1 - normal3 + normal4);
        normal.normalize();

        normals[3 * 3 * wrap2(2 * i)] = normal(0);
        normals[3 * 3 * wrap2(2 * i) + 1] = normal(1);
        normals[3 * 3 * wrap2(2 * i) + 2] = normal(2);

        normals[3 * 3 * wrap2(2 * i + 1)] = normal(0);
        normals[3 * 3 * wrap2(2 * i + 1) + 1] = normal(1);
        normals[3 * 3 * wrap2(2 * i + 1) + 2] = normal(2);

        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle) + 3] = normal(0);
        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle) + 3 + 1] = normal(1);
        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle) + 3 + 2] = normal(2);

        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle - 1) + 3] = normal(0);
        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle - 1) + 3 + 1] = normal(1);
        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle - 1) + 3 + 2] = normal(2);

        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle - 2) + 6] = normal(0);
        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle - 2) + 6 + 1] = normal(1);
        normals[3 * 3 * wrap2(2 * i - 2 * numberOfVerticesPerTubeCircle - 2) + 6 + 2] = normal(2);

        normals[3 * 3 * wrap2(2 * i - 1) + 6] = normal(0);
        normals[3 * 3 * wrap2(2 * i - 1) + 6 + 1] = normal(1);
        normals[3 * 3 * wrap2(2 * i - 1) + 6 + 2] = normal(2);
    }

    n_positions = positions.size();

    glBindVertexArray(vao_);

    // vertex positions -> attribute 0

    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(float), &positions[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);
    

    // normal vectors -> attribute 1

    glBindBuffer(GL_ARRAY_BUFFER, nbo_);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), &normals[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);
}

//-----------------------------------------------------------------------------

void Tube::draw(GLenum mode)
{
    //if (n_indices_ == 0)
    if (filament_.updatedFilament)
    {
        updateBuffers();
        filament_.updatedFilament = false;
    }

    glBindVertexArray(vao_);
    glDrawArrays(mode, 0, n_positions);
    glBindVertexArray(0);
}

//=============================================================================
