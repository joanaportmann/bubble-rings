#version 140

in vec3 v2f_normal;
in vec3 v2f_normal_viewspace;
in vec2 v2f_texcoord_reflect;
in vec2 v2f_texcoord_refract;
in vec3 v2f_light;
in vec3 v2f_view;

out vec4 f_color;

const float shininess = 8.0;
uniform sampler2D tex;
const vec3 sunlight = vec3(1.0, 0.941, 0.898);

void main() {
    vec3 I_in = sunlight;

    // Equal blend of reflected and refracted environment lookups
    vec3 c_reflect = texture(tex, v2f_texcoord_reflect.st).rgb;
    vec3 c_refract = texture(tex, v2f_texcoord_refract.st).rgb;
    vec3 color = 0.5 * c_reflect + 0.5 * c_refract;

    // Add diffuse light
    vec3 n_view = normalize(v2f_normal_viewspace);
    vec3 l_view = normalize(v2f_light);
    vec3 v_view = normalize(v2f_view);
    float nl = dot(n_view, l_view);
    if (nl > 0.0) {
        color += 0.3 * I_in * color * nl;
    }

    // Add specular highlight
    vec3 r = reflect(-l_view, n_view);
    float rv = dot(v_view, r);
    if (nl > 0.0 && rv > 0.0) {
        color += 0.4 * sunlight * pow(rv, shininess);
    }

    f_color = vec4(color, 0.8);
}
