#version 460

in vec3 v2f_normal;
in vec2 v2f_texcoord_reflect;
in vec2 v2f_texcoord_refract;


out vec4 f_color;

uniform sampler2D tex;

void main() {

    vec3 normal_normalized = normalize(v2f_normal);
    vec3 color = 0.01 * vec3(normal_normalized.x, 0.0, 0.0);
    color = 0.6 * texture(tex, v2f_texcoord_reflect.st).rgb;
    color += 0.4 * texture(tex, v2f_texcoord_refract.st).rgb;
//     color = 0.01*texture(tex, v2f_texcoord.st).rgb;
//    color += normal_normalized;
    f_color = vec4(color, 0.7);
}
