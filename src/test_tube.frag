#version 460

flat in vec3 v2f_normal;

out vec4 f_color;

void main() {
    float normal_to_light;
    vec3 normal_normalized = normalize(v2f_normal);
    normal_to_light = dot(normalize(vec3(0, 1, 0)), normal_normalized);
    vec3 color;
    color.x = 0.1 + 0.45 * (normal_normalized.x + 1) * 0.45 * (normal_to_light + 1);
    color.y = 0.1 + 0.45 * (normal_normalized.y + 1) * 0.45 * (normal_to_light + 1);
    color.z = 0.1 + 0.45 * (normal_normalized.z + 1) * 0.45 * (normal_to_light + 1);
    f_color = vec4(color, 1);
}
