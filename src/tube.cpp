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
//#include "tube_viewer.h"

using namespace std;

//=============================================================================

Tube::Tube(std::vector<vec3> CP)
    : controlPolygon_(CP)
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
std::string debugString0("000");
std::string debugString1("111");
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

std::vector<vec3> verticesOfAllCircles(const std::vector<vec3> &controlPolygon, float radius)
{
    std::vector<vec3> verticesOfTube;

    for (int i = 0; i < controlPolygon.size(); i++)
    {
        std::vector<vec3> verticesOfOneCircle = verticesofOneCircle_(
            7,
            controlPolygon[i],
            controlPolygon[(i + 1) % controlPolygon.size()] - controlPolygon[i],
            radius
        );
        for (int j = 0; j < verticesOfOneCircle.size(); j++)
        {
            verticesOfTube.push_back(verticesOfOneCircle[j]);
        }
    };

    return verticesOfTube;
};

//------------------------------------------------------------------------------------

// void angleWeights(const vec3 &p0, const vec3 &p1, const vec3 &p2,
//                   double &w0, double &w1, double &w2)
// {
//     // compute angle weights
//     const vec3 e01 = (p1 - p0).normalized();
//     const vec3 e12 = (p2 - p1).normalized();
//     const vec3 e20 = (p0 - p2).normalized();
//     w0 = acos(std::max(-1.0, std::min(1.0, e01.dot(-e20))));
//     w1 = acos(std::max(-1.0, std::min(1.0, e12.dot(-e01))));
//     w2 = acos(std::max(-1.0, std::min(1.0, e20.dot(-e12))));
// }

//------------------------------------------------------------------------------------

void Tube::createTriangleAndVertexStructs()
{

    tubeVertices = verticesOfAllCircles(controlPolygon_, 0.3);

    for (unsigned int v = 0; v < tubeVertices.size(); ++v)
    {
        Triangle triangle1, triangle2;
        Vertex vertex0, vertex1, vertex2, vertex3;
        bool lastVertexInCircle = v % 7 == 6;

        unsigned int i0 = v;
        unsigned int i1 = (lastVertexInCircle ? v - 6 : v + 1) % tubeVertices.size();
        unsigned int i2 = (v + 7) % tubeVertices.size();
        unsigned int i3 = (lastVertexInCircle ? v + 1 : v + 8) % tubeVertices.size();

        vertex0.position = tubeVertices[i0];
        Tube::vertices_.push_back(vertex0);
        vertex1.position = tubeVertices[i1];
        Tube::vertices_.push_back(vertex1);
        vertex2.position = tubeVertices[i2];
        Tube::vertices_.push_back(vertex2);
        vertex3.position = tubeVertices[i3];
        Tube::vertices_.push_back(vertex3);

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
 * Computes vertex normals by averaging the normals of their incident triangles.
 * Stores the vertex normals in the Vertex::normal member variable.
 * Weights the normals by their triangles' angles.
 */

void Tube::compute_normals()
{

    // compute triangle normals
    for (Triangle &t : triangles_)
    {
        const vec3 &p0 = vertices_.at(t.ind0).position;
        const vec3 &p1 = vertices_.at(t.ind1).position;
        const vec3 &p2 = vertices_.at(t.ind2).position;

        t.normal = ((p1 - p0).cross(p2 - p0)).normalized();
    }
    // initialize vertex normals to zero
    // for (Vertex &v : vertices_)
    // {
    //     v.normal = vec3(0, 0, 0);
    // }

    // for (Triangle &t : triangles_)
    // {
    //     const vec3 &p0 = vertices_.at(t.ind0).position;
    //     const vec3 &p1 = vertices_.at(t.ind1).position;
    //     const vec3 &p2 = vertices_.at(t.ind2).position;

    //     double w0, w1, w2;

    //     // Weigh the normals by their triangles' angles.
    //     angleWeights(p0, p1, p2, w0, w1, w2);
   
    //     //adding the normals all together

    //     vertices_.at(t.ind0).normal += t.normal * w0;
    //     vertices_.at(t.ind1).normal += t.normal * w1;
    //     vertices_.at(t.ind2).normal += t.normal * w2;
    // }


    // for (Vertex &v : vertices_)
    // {
    //     v.normal = (v.normal).normalized();
    //     //cout << v.position << "\n"; 
    // } 
}

//--------------------------------------------------------------------------------------

void Tube::initialize()
{

    Tube::createTriangleAndVertexStructs();
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

 
    int counter = 0;
    // generate triangles
    for (Triangle &t : triangles_)
    {
        //if (counter++ > 2) break;

        indices[i++] = t.ind0;
        indices[i++] = t.ind1;
        indices[i++] = t.ind2;

        normals[n++] = t.normal(0);
        normals[n++] = t.normal(1);
        normals[n++] = t.normal(2);
    }

    n_indices_ = 3 * triangles_.size();

    
    // generate vertex array object
    glGenVertexArrays(1, &vao_);
  
    // generate buffers
    glGenBuffers(1, &vbo_);
    glGenBuffers(1, &ibo_);
    glGenBuffers(1, &nbo_);

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

    // texture coordinates -> attribute 2
    glGenBuffers(1, &tbo_);
    glBindBuffer(GL_ARRAY_BUFFER, tbo_);
    glBufferData(GL_ARRAY_BUFFER, texcoords.size() * sizeof(float), &texcoords[0], GL_STATIC_DRAW);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(2);

    // triangle indices

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices[0], GL_STATIC_DRAW);
}

//-----------------------------------------------------------------------------

void Tube::draw(GLenum mode)
{
    if (n_indices_ == 0)
        initialize();

    glBindVertexArray(vao_);
    glDrawElements(mode, n_indices_, GL_UNSIGNED_INT, NULL);
    glBindVertexArray(0);
}

//=============================================================================
