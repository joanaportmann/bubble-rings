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

//=============================================================================


Tube::Tube(std::vector<vec3> tubeVertices){}


//-----------------------------------------------------------------------------


Tube::~Tube()
{
    if (vbo_)  glDeleteBuffers(1, &vbo_);
    if (nbo_)  glDeleteBuffers(1, &nbo_);
    if (tan_bo_)  glDeleteBuffers(1, &tan_bo_);
    if (bitan_bo_)  glDeleteBuffers(1, &bitan_bo_);
    if (tbo_)  glDeleteBuffers(1, &tbo_);
    if (ibo_)  glDeleteBuffers(1, &ibo_);
    if (vao_)  glDeleteVertexArrays(1, &vao_);
}


//-----------------------------------------------------------------------------


void Tube::initialize()
{
    const unsigned int n_vertices   = 3 * tubeVertices.size();
    // Hardcoded: Momentan hat Kreis 7 vertces
    const unsigned int n_triangles  = 2 * tubeVertices.size();

    std::vector<GLfloat> positions(3*n_vertices);
    std::vector<GLuint >   indices(3*n_triangles);

    unsigned int p(0), i(0);
    //unsigned int t(0), n(0), tan(0), bitan(0);

    // generate vertices
    for (unsigned int i=0; i<tubeVertices.size(); ++i)
    {
            vec3 currentVector = tubeVertices[i];
            positions[p++] = currentVector(0);
            positions[p++] = currentVector(1);
            positions[p++] = currentVector(2);

            //normals[n++] = x;
            //normals[n++] = y;
            //normals[n++] = z;

            //texcoords[t++] = 1.0-u;
            //texcoords[t++] = 1.0-v;
    }


    // generate triangles
    for (unsigned int v=0; v<tubeVertices.size()-1; ++v)
    {       
        bool lastVertexInCircle = v % 7 == 6;
            unsigned int i0 = v;
            unsigned int i1 = (lastVertexInCircle ? v - 6 : v + 1) % tubeVertices.size();
            unsigned int i2 = (v+7) % tubeVertices.size();
            unsigned int i3 = (lastVertexInCircle ? v + 1 : v + 8) % tubeVertices.size();

            indices[i++] = i0;
            indices[i++] = i1;
            indices[i++] = i2;

            indices[i++] = i1;
            indices[i++] = i2;
            indices[i++] = i3;
        
    }
    n_indices_ = 3*n_triangles;


    // generate vertex array object
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);


    // vertex positions -> attribute 0
    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, 3*n_vertices*sizeof(float), &positions[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    // triangle indices
    glGenBuffers(1, &ibo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*n_triangles*sizeof(GLuint), &indices[0], GL_STATIC_DRAW);
}


//-----------------------------------------------------------------------------


void Tube::draw(GLenum mode)
{
    if (n_indices_ == 0) initialize();

    glBindVertexArray(vao_);
    glDrawElements(mode, n_indices_, GL_UNSIGNED_INT, NULL);
    glBindVertexArray(0);
}


//=============================================================================
