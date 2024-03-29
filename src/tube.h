
//=============================================================================
#ifndef TUBE_H
#define TUBE_H
//=============================================================================

#include "gl.h"
#include <Eigen/Dense>
#include "glmath.h"
#include <vector>
#include "filament.h"
#include "texture.h"
//=============================================================================


class Tube
{

public:

    Tube(Filament &filament);

    /// destructor
    ~Tube();

    /// render mesh of the sphere
    void draw(GLenum mode=GL_TRIANGLES);

    void initialize();
    Texture tex_;
    




private:

    void updateBuffers();

    /// Compute normal vectors for triangles and vertices
    void compute_normals();

    std::vector<vec3> circleVertices_t(int n, vec3 normal);

    // generate triangle structs
    void createTriangleStruct();
    int wrap2(int i);

    /// indices of the triangle vertices
    unsigned int n_positions = 0;

    

    // vertex array object
    GLuint vao_ = 0;
    /// vertex buffer object
    GLuint vbo_ = 0;
    /// normals buffer object
    GLuint nbo_ = 0;
    /// tangents buffer object
    GLuint tan_bo_ = 0;
    /// bitangents buffer object
    GLuint bitan_bo_ = 0;
    /// texture coordinates buffer object
    GLuint tbo_ = 0;
    /// index buffer object
    GLuint ibo_ = 0;

    GLuint lineVAO, lineVBO;

    /// a triangle is specified by three indices and a normal
    struct Triangle
    {
        /// index of first vertex 
        int ind0;
        /// index of second vertex 
        int ind1;
        /// index of third vertex 
        int ind2;
        /// triangle normal
        vec3 normal;
    };

    /// Array of triangles
    std::vector<Triangle> triangles_;

    std::vector<vec3> tubeVertices;
    Filament &filament_;
    
    

};


//=============================================================================
#endif
//=============================================================================
