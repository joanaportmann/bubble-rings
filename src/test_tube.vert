
#version 140
#extension GL_ARB_explicit_attrib_location : enable

layout (location = 0) in vec4 v_position;
layout (location = 1) in vec3 v_normal;
layout (location = 2) in vec2 v_texcoord;

//out vec2 v2f_texcoord;
out vec3 v2f_normal;
//out vec3 v2f_light;
//out vec3 v2f_view;

uniform mat4 modelview_projection_matrix;
uniform mat3 normal_matrix;
uniform vec4 light_position; //in eye space coordinates already


void main() {
    // Compute vertices' normalized device coordinates
    gl_Position = modelview_projection_matrix * v_position;
    	
    v2f_normal = normalize(normal_matrix * v_normal);
}