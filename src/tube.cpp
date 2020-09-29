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
//#include "tube_viewer.h"

//=============================================================================


Tube::Tube(std::vector<vec3> CP)
    : controlPolygon_(CP)
{
}


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

std::vector<vec3> circleVertices_t(int n, vec3 center, vec3 normal, float radius)
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

std::vector<vec3> verticesOfAllCircles( std::vector<vec3> controlPolygon, float radius, std::vector<vec3> verticesOfTube)
{
	for (int i = 1; i < controlPolygon.size() - 1; i++)
	{
		std::vector<vec3> verticesOfOneCircle = circleVertices_t(12, controlPolygon[i], controlPolygon[i + 1] - controlPolygon[i], radius);
		for (int i = 0; i < verticesOfOneCircle.size() - 1; i++)
		{
			verticesOfTube.push_back(verticesOfOneCircle[i]);
		}
		
	};
	return verticesOfTube;
};

//------------------------------------------------------------------------------------


void Tube::initialize()
{
    Tube::tubeVertices = verticesOfAllCircles(controlPolygon_ , 0.3, Tube::tubeVertices);
    const unsigned int n_vertices   = 3 * Tube::tubeVertices.size();
    // Hardcoded: Momentan hat Kreis 7 vertces
    const unsigned int n_triangles  = 2 * Tube::tubeVertices.size();

    std::vector<GLfloat> positions(3*n_vertices);
    std::vector<GLuint >   indices(3*n_triangles);

    unsigned int p(0), i(0);
    //unsigned int t(0), n(0), tan(0), bitan(0);

    // generate vertices
    for (unsigned int i=0; i<tubeVertices.size() - 1; ++i)
    {
            vec3 currentVector = Tube::tubeVertices[i];
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
    for (unsigned int v=0; v<Tube::tubeVertices.size()-1; ++v)
    {       
        bool lastVertexInCircle = v % 7 == 6;
            unsigned int i0 = v;
            unsigned int i1 = (lastVertexInCircle ? v - 6 : v + 1) % Tube::tubeVertices.size();
            unsigned int i2 = (v+7) % Tube::tubeVertices.size();
            unsigned int i3 = (lastVertexInCircle ? v + 1 : v + 8) % Tube::tubeVertices.size();

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
