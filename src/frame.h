#ifndef FRAME_H
#define FRAME_H

#include "shader.h"
#include <Eigen/Geometry>

struct Frame {
    Frame() {}

   
    // Initialize the OpenGL arrays/buffers with the arrow geometry for frame
    // visualization
    void initialize() {
        // generate vertex array object
        glGenVertexArrays(1, &m_arrowVAO);
        glBindVertexArray(m_arrowVAO);

        std::vector<GLfloat> coords;
        std::vector<GLint> idx;

        // Arrow geometry along the x axis.
        // Shaft cylinder
        const size_t nsubdiv = 20;
        const float radius = 0.02;
        for (size_t i = 0; i < nsubdiv; ++i) {
            float theta = (2 * M_PI * i) / nsubdiv;
            coords.push_back(0);
            coords.push_back( radius * cos(theta));
            coords.push_back(-radius * sin(theta));
            coords.push_back(0.75);
            coords.push_back( radius * cos(theta));
            coords.push_back(-radius * sin(theta));
            idx.push_back(2 * i                   + 0);
            idx.push_back(2 * i                   + 1);
            idx.push_back(2 * ((i + 1) % nsubdiv) + 1);
            idx.push_back(2 * i                   + 0);
            idx.push_back(2 * ((i + 1) % nsubdiv) + 1);
            idx.push_back(2 * ((i + 1) % nsubdiv) + 0);
        }
        // Arrow head cone
        size_t tipIdx = coords.size()  / 3;
        coords.push_back(1.0); coords.push_back(0); coords.push_back(0);

        size_t offset = tipIdx + 1;
        for (size_t i = 0; i < nsubdiv; ++i) {
            float theta = (2 * M_PI * i) / nsubdiv;
            coords.push_back(0.75);
            coords.push_back( radius * 1.5 * cos(theta));
            coords.push_back(-radius * 1.5 * sin(theta));

            idx.push_back(offset + i);
            idx.push_back(tipIdx);
            idx.push_back(offset + (i + 1) % nsubdiv);
        }

        // vertex positions -> attribute 0
        glGenBuffers(1, &m_arrowVBO);
        glBindBuffer(GL_ARRAY_BUFFER, m_arrowVBO);
        glBufferData(GL_ARRAY_BUFFER, coords.size() * sizeof(GLfloat), coords.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(0);

        // triangle indices
        glGenBuffers(1, &m_arrowIBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_arrowIBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(GLuint), idx.data(), GL_STATIC_DRAW);

        m_nidx = idx.size();

        glBindVertexArray(0);
    }

    void draw(Shader &s, const mat4 &vp, const vec3 &pos) {
        glBindVertexArray(m_arrowVAO);

        Eigen::Affine3f scaling;
        scaling = Eigen::Scaling(1.0f);
        mat4 scaling_matrix = scaling.matrix().cast<double>();

        Eigen::Affine3f rotation_z;
        rotation_z = Eigen::AngleAxisf(-0.5 * M_PI, vec3::UnitZ().cast<float>());
        mat4 rotation_z_matrix = rotation_z.matrix().cast<double>();

        Eigen::Affine3f rotation_y;
        rotation_y = Eigen::AngleAxisf(0.5 * M_PI, vec3::UnitY().cast<float>());
        mat4 rotation_y_matrix = rotation_y.matrix().cast<double>();

        mat4 M_1, M_2;
        mat4 M_0 = vp * scaling_matrix;
        M_1 = vp * rotation_z_matrix * scaling_matrix;
        M_2 = vp * rotation_y_matrix * scaling_matrix;
        s.use();
        s.set_uniform<mat4>("modelview_projection_matrix", M_0); 
        s.set_uniform("color", vec4(1.0, 0.0, 0.0, 1.0)); 
        glDrawElements(GL_TRIANGLES, m_nidx, GL_UNSIGNED_INT, NULL); // x-Achse

        s.set_uniform("color", vec4(0.0, 1.0, 0.0, 1.0));
        s.set_uniform("modelview_projection_matrix", M_1);
        glDrawElements(GL_TRIANGLES, m_nidx, GL_UNSIGNED_INT, NULL); // y-Achse

        s.set_uniform("color", vec4(0.0, 0.0, 1.0, 1.0));
        s.set_uniform("modelview_projection_matrix", M_2);
        glDrawElements(GL_TRIANGLES, m_nidx, GL_UNSIGNED_INT, NULL); // z -Achse

        glBindVertexArray(0);
    }

    ~Frame() {
        if (m_arrowVAO != 0) glDeleteVertexArrays(1, &m_arrowVAO);
        if (m_arrowVBO != 0) glDeleteBuffers(1, &m_arrowVBO);
        if (m_arrowIBO != 0) glDeleteBuffers(1, &m_arrowIBO);
    }
private:
    GLuint m_arrowVAO = 0,
           m_arrowVBO = 0,
           m_arrowIBO = 0;
    size_t m_nidx = 0;
    bool m_useParallelTransport = false;
};

#endif /* end of include guard: FRAME_H */
