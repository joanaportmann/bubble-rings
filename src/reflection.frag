#version 460

flat in vec3 v2f_normal;
in vec2 v2f_texcoord;
in vec3 eye;


out vec4 f_color;

uniform sampler2D tex;

void main() {

    vec3 color = vec3(0.0,0.0,0.0);

	color += texture(tex, v2f_texcoord.st).rgb;
    f_color = vec4(color, 1);
}
