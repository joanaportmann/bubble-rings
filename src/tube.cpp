//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================

#include "tube.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdlib.h> 

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
    tubeVertices = filament_.getBubbleRingSkeleton();

    for (unsigned int v = 0; v < tubeVertices.size(); ++v)
    {
        Triangle triangle1, triangle2;
       
        bool lastVertexInCircle = v % numberOfVerticesPerTubeCircle == numberOfVerticesPerTubeCircle - 1;

        unsigned int i0 = v;
        unsigned int i1 = (lastVertexInCircle ? v - (numberOfVerticesPerTubeCircle - 1) : v + 1) % tubeVertices.size();
        unsigned int i2 = (v + numberOfVerticesPerTubeCircle) % tubeVertices.size();
        unsigned int i3 = (lastVertexInCircle ? v + 1 : v + (numberOfVerticesPerTubeCircle + 1)) % tubeVertices.size();


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

void Tube::initialize(){
     // generate vertex array object
    glGenVertexArrays(1, &vao_);
  
    // generate buffers
    glGenBuffers(1, &vbo_);
    glGenBuffers(1, &nbo_);
    glGenBuffers(1, &ibo_);
}

//------------------------------------------------------------------------------------

void Tube::updateBuffers()
{

    Tube::createTriangleStruct();
    Tube::compute_normals();

    std::vector<GLfloat> positions(3 * tubeVertices.size());
    std::vector<GLuint> indices(3 * triangles_.size());
    std::vector<GLfloat> normals(3 * triangles_.size());
    std::vector<GLfloat> texcoords(2 * tubeVertices.size());

    unsigned int p(0), i(0), n(0), t(0);

    // generate vertices
    for (int k = 0; k < tubeVertices.size(); k++)
    {
        positions[p++] = tubeVertices[k](0);
        positions[p++] = tubeVertices[k](1);
        positions[p++] = tubeVertices[k](2);

        texcoords[t++] = 0.8;
        texcoords[t++] = 0.5;
    }


    // generate triangles
    int count = 0;
    for (Triangle &t : triangles_)
    {
        indices[i++] = t.ind0;
        indices[i++] = t.ind1;
        indices[i++] = t.ind2;

        if (count++ % 2 == 1) continue;

        normals[n++] = t.normal(0);
        normals[n++] = t.normal(1);
        normals[n++] = t.normal(2);
    }

    n_indices_ = 3 * triangles_.size();

    
 

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


    // triangle indices

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices[0], GL_STATIC_DRAW);
}

//-----------------------------------------------------------------------------

void Tube::draw(GLenum mode)
{
   //if (n_indices_ == 0)
    updateBuffers();

    glBindVertexArray(vao_);
    glDrawElements(mode, n_indices_, GL_UNSIGNED_INT, NULL);
    glBindVertexArray(0);
}

//=============================================================================
