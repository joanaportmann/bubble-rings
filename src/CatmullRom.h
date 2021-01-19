#ifndef CATMULL_ROM_H
#define CATMULL_ROM_H

#include <Eigen/Dense>
#include "glmath.h"
#include <vector>

class CatmullRom
{
public:
	CatmullRom();
	~CatmullRom();

	void add_way_point(vec3 point);
	vec3 node(int i) const { return _nodes[i]; }
	double length_from_starting_point(int i) const { return _distances[i]; }
	bool has_next_node(int i) const { return static_cast<int>(_nodes.size()) > i; }
	int node_count() const {  return static_cast<int>(_nodes.size()); }
	bool is_empty() { return _nodes.empty(); }
	double total_length() const
	{
		assert(!_distances.empty());
		return _distances[_distances.size() - 1];
	}

	void increment_steps(int steps) { _steps+=steps; }
	void set_steps(int steps) { _steps = steps; }

private:
	void _on_way_point_added();
	void add_node( vec3 node);
	void clear();

	vec3 interpolate(double u, const vec3& P0, const vec3& P1, const vec3& P2, const vec3& P3);
	std::vector<vec3> _way_points;
	std::vector<vec3> _nodes;
	std::vector<double> _distances;
	int _steps;


};

#endif