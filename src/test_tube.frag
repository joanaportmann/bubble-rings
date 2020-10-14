//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================

#version 140

flat in vec3 v2f_normal;

out vec4 f_color;

void main() {
    // f_color = color;
    float normal_to_light;
    normal_to_light = 0.5 * dot(vec3(0.2, 0.3, 0.1), v2f_normal); 
    vec3 color;
    color.x = 0.5 + normal_to_light;
    color.y = 0.5 + normal_to_light;
    color.z = 0.5 + normal_to_light;
    f_color = vec4(color, 1);
}
