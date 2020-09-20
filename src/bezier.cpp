#include "bezier.h"
#include <utility>
#include <cmath>
#include <array>

// Calculate a point one of the Bezier curve segments
// @param bp_offset  index into bezier_control_points_ of the first of four
//                   control points defining the Bezier segment we want to evaluate.
// @param t          parametric distance along the curve at which to evaluate
vec3 PiecewiseBezier::eval_bezier(int bp_offset, float t) const {
	//Casteljau-Algorithmus
	vec3 c00 = mix(bezier_control_points_[bp_offset], bezier_control_points_[bp_offset + 1], t);
	vec3 c01 = mix(bezier_control_points_[bp_offset + 1], bezier_control_points_[bp_offset + 2], t);
	vec3 c02 = mix(bezier_control_points_[bp_offset + 2], bezier_control_points_[bp_offset + 3], t);

	vec3 c10 = mix(c00, c01, t);
	vec3 c11 = mix(c01, c02, t);

	return mix(c10, c11, t);
}

// Calculate a tangent at point at one of the Bezier curve segments
// @param bp_offset  index into bezier_control_points_ of the first of four
//                   control points defining the Bezier segment we want to compute
//                   the tangent at
// @param t          parametric distance along the curve at which to evaluate
vec3 PiecewiseBezier::eval_bezier_tangent(int bp_offset, float t) const {

	// The vector between the two last points of the Casteljau Alg is the tangent divided by n
	vec3 c00 = mix(bezier_control_points_[bp_offset], bezier_control_points_[bp_offset + 1], t);
	vec3 c01 = mix(bezier_control_points_[bp_offset + 1], bezier_control_points_[bp_offset + 2], t);
	vec3 c02 = mix(bezier_control_points_[bp_offset + 2], bezier_control_points_[bp_offset + 3], t);

	vec3 c10 = mix(c00, c01, t);
	vec3 c11 = mix(c01, c02, t);

	return 3 * (c11 - c10);
}

std::vector<vec3> PiecewiseBezier::control_polygon_to_bezier_points(std::vector<vec3> const& cp) {
	std::vector<vec3> bezier_pts;
	size_t numSegments = cp.size() - 3;

	//Debugging
	std::cout << numSegments << "segments\n";

	for (auto&& p : bezier_pts) {
		std::cout << p.x << ", " << p.y << ", " << p.z << "\n";
	}
	std::cout << "----------------\n";

	for (int i = 1; i <= cp.size() - 2; i++) {
		vec3 prev_2_thirds = mix(cp[i - 1], cp[i], 2.0f / 3.0f);
		vec3 next_1_third = mix(cp[i], cp[i + 1], 1.0f / 3.0f);

		if (i != 1) bezier_pts.emplace_back(prev_2_thirds);
		bezier_pts.emplace_back(mix(prev_2_thirds, next_1_third, 0.5));
		if (i != cp.size() - 2) bezier_pts.emplace_back(next_1_third);
	}

	//Debugging
	for (auto&& p : bezier_pts) {
		std::cout << p.x << ", " << p.y << ", " << p.z << "\n";
	}

	return bezier_pts;
}

vec3 PiecewiseBezier::eval_piecewise_bezier_curve(float t) const {
	if (t == 1) return eval_bezier(3 * (num_segments() - 1), 1.f);

	float t_ = t * num_segments();
	int segment_index = floor(t_);
	return eval_bezier(3 * segment_index, t_ - float(segment_index));
}

vec3 PiecewiseBezier::operator()(float t) const {
	return eval_piecewise_bezier_curve(t);
}

vec3 PiecewiseBezier::tangent(float t) const {
	if (t == 1) return eval_bezier(3 * (num_segments() - 1), 1.f);

	float t_ = t * num_segments();
	int segment_index = floor(t_);
	return eval_bezier_tangent(3 * segment_index, t_ - float(segment_index));
}

void PiecewiseBezier::set_control_polygon(const std::vector<vec3>& control_polygon, bool loop) {
	control_polygon_ = control_polygon;
	if (loop) {
		control_polygon_.push_back(control_polygon[0]);
		control_polygon_.push_back(control_polygon[1]);
		control_polygon_.push_back(control_polygon[2]);
	}

	std::cout << control_polygon.size() << " (original)\n";
	std::cout << control_polygon_.size() << " (extended)\n";
	bezier_control_points_ = control_polygon_to_bezier_points(control_polygon_);
}
