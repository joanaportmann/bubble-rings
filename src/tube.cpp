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
    tubeVertices.clear();
    tubeVertices = filament_.getBubbleRingSkeleton();

    for (unsigned int v = 0; v < tubeVertices.size(); ++v)
    {
        Triangle triangle1, triangle2;

        bool lastVertexInCircle = (v % numberOfVerticesPerTubeCircle == numberOfVerticesPerTubeCircle - 1);

        unsigned int i0 = v;
        unsigned int i1 = (lastVertexInCircle ? (v - (numberOfVerticesPerTubeCircle - 1)) : (v + 1)) % tubeVertices.size();
        unsigned int i2 = (v + numberOfVerticesPerTubeCircle) % tubeVertices.size();
        unsigned int i3 = (lastVertexInCircle ? (v + 1) : v + (numberOfVerticesPerTubeCircle + 1)) % tubeVertices.size();

        // Counter-clockwise winding when viewed from outside the tube (outward normals)
        triangle1.ind0 = i0;
        triangle1.ind1 = i1;
        triangle1.ind2 = i3;
        Tube::triangles_.push_back(triangle1);

        triangle2.ind0 = i0;
        triangle2.ind1 = i3;
        triangle2.ind2 = i2;
        Tube::triangles_.push_back(triangle2);
    }
}

//------------------------------------------------------------------------------------

/**
 * Computes outward face normals for all triangles
 */
void Tube::compute_normals()
{
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
    Tube::compute_normals();

    // Accumulate and normalize the six incident triangle face normals per tube vertex
    std::vector<vec3> vertexNormals(tubeVertices.size(), vec3::Zero());
    for (const Triangle &tri : triangles_)
    {
        vertexNormals[tri.ind0] += tri.normal;
        vertexNormals[tri.ind1] += tri.normal;
        vertexNormals[tri.ind2] += tri.normal;
    }
    for (vec3 &vn : vertexNormals)
    {
        if (vn.norm() > 1e-12)
        {
            vn.normalize();
        }
    }

    std::vector<GLfloat> positions(3 * 3 * triangles_.size());
    std::vector<GLfloat> normals(3 * 3 * triangles_.size());

    unsigned int p = 0;
    unsigned int n_idx = 0;

    for (const Triangle &tri : triangles_)
    {
        const int indices[3] = {tri.ind0, tri.ind1, tri.ind2};
        for (int k = 0; k < 3; ++k)
        {
            const vec3 &pos = tubeVertices[indices[k]];
            const vec3 &nrm = vertexNormals[indices[k]];
            positions[p++] = static_cast<GLfloat>(pos(0));
            positions[p++] = static_cast<GLfloat>(pos(1));
            positions[p++] = static_cast<GLfloat>(pos(2));

            normals[n_idx++] = static_cast<GLfloat>(nrm(0));
            normals[n_idx++] = static_cast<GLfloat>(nrm(1));
            normals[n_idx++] = static_cast<GLfloat>(nrm(2));
        }
    }

    // Number of vertices (3 floats per vertex) passed to glDrawArrays
    n_positions = static_cast<unsigned int>(positions.size() / 3);

    glBindVertexArray(vao_);

    // vertex positions -> attribute 0
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(GLfloat), positions.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    // normal vectors -> attribute 1
    glBindBuffer(GL_ARRAY_BUFFER, nbo_);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(GLfloat), normals.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(1);
}

//-----------------------------------------------------------------------------

void Tube::draw(GLenum mode)
{
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
