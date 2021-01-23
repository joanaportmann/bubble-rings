#include "CatmullRom.h"
#include <iostream>
#include <Eigen/Dense>
#include "glmath.h"
#include <vector>

using namespace std;

CatmullRom::CatmullRom()
	: _steps(3)
{
}

CatmullRom::~CatmullRom() {}

void CatmullRom::clear()
{
	_nodes.clear();
	_way_points.clear();
	_distances.clear();
}

void CatmullRom::add_way_point(vec3 point)
{
	_way_points.push_back(point);
	_on_way_point_added();
}

void CatmullRom::add_node(vec3 node)
{
	_nodes.push_back(node);

	if (_nodes.size() == 1)
	{
		_distances.push_back(0);
	}
	else
	{
		int new_node_index = _nodes.size() - 1;

		double segment_distance = (_nodes[new_node_index] - _nodes[new_node_index - 1]).norm();
		_distances.push_back(segment_distance + _distances[new_node_index - 1]);
	}
}

void CatmullRom::_on_way_point_added()
{
	if (_way_points.size() < 4)
	{
		return;
	}

	int new_control_point_index = _way_points.size() - 1;
	int pt = new_control_point_index - 2;
	for (int i = 0; i <= _steps; i++)
	{
		double u = (double)i / (double)_steps;

		add_node(interpolate(u, _way_points[pt - 1], _way_points[pt], _way_points[pt + 1], _way_points[pt + 2]));
	}
}

vec3 CatmullRom::interpolate(double u, const vec3 &P0, const vec3 &P1, const vec3 &P2, const vec3 &P3)
{
	vec3 point;
	point = u * u * u * ((-1) * P0 + 3 * P1 - 3 * P2 + P3) / 2;
	point += u * u * (2 * P0 - 5 * P1 + 4 * P2 - P3) / 2;
	point += u * ((-1) * P0 + P2) / 2;
	point += P1;

	return point;
}