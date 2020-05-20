#version 120

attribute vec4 vPosition;
attribute vec4 vColor;
attribute vec4 vNormal;
varying vec4 color;

uniform mat4 cubie_ctm; // Uniform current transformation matrix
uniform mat4 global_ctm;
uniform mat4 model_view; // Uniform model view matrix
uniform mat4 projection; // Uniform projection matrix

uniform vec4 light_position = vec4(0, 0, 1.5, 1);
vec4 ambient, diffuse, specular, reflection;

void main()
{
	ambient = vColor * 0.3;
	vec4 N = normalize(projection * model_view * global_ctm * cubie_ctm * vNormal);
	vec4 L_temp = projection * model_view * (light_position - (global_ctm * cubie_ctm * vPosition));
	vec4 L = normalize(L_temp);
	diffuse = max(dot(L,N), 0.0) * vColor;

	vec4 eye_point = vec4(0.0, 0.0, 1.5, 1.0);
	vec4 V = normalize(projection * model_view * (eye_point - (global_ctm * cubie_ctm * vPosition)));
	vec4 H = normalize(L + V);
	specular = pow(max(dot(N, H), 0.0), 80) * vec4(1.0, 1.0, 1.0, 1.0);

	float dist = length(L_temp);
	float attenuation = 1 / (0.2 + (0.2 * dist) + (0.25 * dist * dist));

	gl_Position = projection * model_view * global_ctm * cubie_ctm * vPosition / vPosition.w;
	color = ambient + (attenuation * (diffuse + specular));
}
