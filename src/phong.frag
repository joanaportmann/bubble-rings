//=============================================================================
//
//   Exercise code for the lecture "Introduction to Computer Graphics"
//     by Prof. Mario Botsch, Bielefeld University
//
//   Copyright (C) by Computer Graphics Group, Bielefeld University
//
//=============================================================================

#version 140

in vec3 v2f_normal;

out vec4 f_color;

void main() {
    vec3 normal = normalize(v2f_normal);
    f_color.x = normal.x;
    f_color.y = normal.y;
    f_color.z = normal.z;
    f_color.w = 0.3;
}
