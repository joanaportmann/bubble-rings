
#version 140
#extension GL_ARB_explicit_attrib_location : enable

#define M_PI 3.1415926535897932384626433832795

layout (location = 0) in vec4 v_position;
layout (location = 1) in vec3 v_normal;
layout (location = 2) in vec2 v_texcoord;

out vec2 v2f_texcoord_reflect;
out vec2 v2f_texcoord_refract;
out vec3 v2f_normal;
out vec3 v2f_normal_viewspace;
out vec2 v2f_texcoord;
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

    // reflection
    vec3 reflected_normal = normalize(reflect(eyePositionW - vec3(v_position), v2f_normal));
    vec2 reflected_normal_xz = normalize(vec2(reflected_normal.x, reflected_normal.z));
    float angle = reflected_normal.z <= 0.0 ? acos(-reflected_normal_xz.x) : 2 * M_PI - acos(-reflected_normal_xz.x);
    v2f_texcoord_reflect = vec2(angle / M_PI / 2, -asin(reflected_normal.y) / M_PI + 0.5);

    // refraction
    vec3 refracted_normal = normalize(refract(normalize(eyePositionW - vec3(v_position)), -v2f_normal, 1.0 / 1.3));
    vec2 refracted_normal_xz = normalize(vec2(refracted_normal.x, refracted_normal.z));
    float angleRefract = refracted_normal.z <= 0.0 ? acos(-refracted_normal_xz.x) : 2 * M_PI - acos(-refracted_normal_xz.x);
    v2f_texcoord_refract = vec2(angleRefract / M_PI / 2, -asin(refracted_normal.y) / M_PI + 0.5);

    // Phong 
    v2f_texcoord = v_texcoord;
	v2f_light = normalize(vec3(light_position ));
    v2f_view = normalize(vec3(modelview_matrix * v_position));
	v2f_normal_viewspace = normalize(normal_matrix * v_normal);
	
}
