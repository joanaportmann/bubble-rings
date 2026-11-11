#version 140
#extension GL_ARB_explicit_attrib_location : enable

#define M_PI 3.1415926535897932384626433832795

layout (location = 0) in vec4 v_position;
layout (location = 1) in vec3 v_normal;

out vec2 v2f_texcoord_reflect;
out vec2 v2f_texcoord_refract;
out vec3 v2f_normal;
out vec3 v2f_normal_viewspace;
out vec3 v2f_light;
out vec3 v2f_view;

uniform vec3 eyePositionW;
uniform mat4 modelview_projection_matrix;
uniform mat3 normal_matrix;
uniform vec3 light_position;
uniform mat4 modelview_matrix;

void main() {
    // Compute vertices' normalized device coordinates
    gl_Position = modelview_projection_matrix * v_position;
    v2f_normal = normalize(v_normal);

    vec3 incident = normalize(vec3(v_position) - eyePositionW);

    // reflection
    vec3 reflected_normal = normalize(reflect(incident, v2f_normal));
    vec2 reflected_normal_xz = normalize(vec2(reflected_normal.x, reflected_normal.z));
    float angle = reflected_normal.z <= 0.0 ? acos(-reflected_normal_xz.x) : 2.0 * M_PI - acos(-reflected_normal_xz.x);
    v2f_texcoord_reflect = vec2(angle / (2.0 * M_PI), -asin(reflected_normal.y) / M_PI + 0.5);

    // refraction (with total internal reflection fallback)
    vec3 refracted_raw = refract(incident, v2f_normal, 1.0 / 1.33);
    vec3 refracted_normal = length(refracted_raw) < 1e-5 ? reflected_normal : normalize(refracted_raw);
    vec2 refracted_normal_xz = normalize(vec2(refracted_normal.x, refracted_normal.z));
    float angleRefract = refracted_normal.z <= 0.0 ? acos(-refracted_normal_xz.x) : 2.0 * M_PI - acos(-refracted_normal_xz.x);
    v2f_texcoord_refract = vec2(angleRefract / (2.0 * M_PI), -asin(refracted_normal.y) / M_PI + 0.5);

    // Phong lighting vectors in view space
    v2f_light = normalize(light_position);
    v2f_view = normalize(-vec3(modelview_matrix * v_position));
    v2f_normal_viewspace = normalize(normal_matrix * v_normal);
}
