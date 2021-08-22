#version 460

in vec3 v2f_normal;
in vec3 v2f_normal_viewspace;
in vec2 v2f_texcoord_reflect;
in vec2 v2f_texcoord_refract;
in vec2 v2f_texcoord;
in vec3 v2f_light;
in vec3 v2f_view;

out vec4 f_color;

const float shininess = 8.0;
// const vec3  sunlight = vec3(1.0, 0.941, 0.898);
uniform sampler2D tex;

const vec3  sunlight = vec3(1.0, 0.941, 0.898);

void main() {


	vec3 c_material = texture(tex, v2f_texcoord.st).rgb;
	vec3 I_a = 0.2 * sunlight;
	vec3 I_in = sunlight;

	//TODO weight refraction and reflection based on angle

    vec3 normal_normalized = normalize(v2f_normal);
    vec3 color = 0.01 * vec3(normal_normalized.x, 0.0, 0.0);   
    color =  0.5 * texture(tex, v2f_texcoord_reflect.st).rgb;
    color += 0.5 * texture(tex, v2f_texcoord_refract.st).rgb;

	// Add ambient light
	color += 0 * I_a * c_material;

	// Add diffuse light
    float nl = dot(v2f_normal_viewspace, v2f_light);
	if (nl > 0) {
		color += 0.3 * I_in * c_material * nl;
	}

	// Add specular light
	vec3 r = reflect(v2f_light, v2f_normal_viewspace);
	float rv = dot(v2f_view, r);
	if (nl > 0 && rv > 0) {
		color +=  0.4 * sunlight * c_material * pow(dot(r, v2f_view), shininess);
	}

 
    // add required alpha value
    f_color = vec4(color, 0.8);

}
