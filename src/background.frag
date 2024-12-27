#version 140

in  vec2 v2f_texcoord;
out vec4 f_color;

uniform sampler2D tex;

void main()
{
    // fetch color (rgb and alpha) from texture
    vec4 color = texture(tex, v2f_texcoord.st);
    f_color = color;
}
