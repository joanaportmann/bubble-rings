#pragma once

#include <iostream>
#include <cassert>
#include <math.h>
#include <Eigen/Dense>

using vec4 = Eigen::Vector4d;
using vec3 = Eigen::Vector3d;
using mat3 = Eigen::Matrix3d;
using mat4 = Eigen::Matrix4d;

class MatUtils
{
public:
    /// return frustum matrix
    static mat4 frustum(float left, float right, float bottom, float top, float near, float far);
    /// return matrix for perspective projection (special case of frustum matrix)
    static mat4 perspective(float fovy, float aspect, float near, float far);
    /// return look-at camera matrix
    static mat4 look_at(const vec4& eye, const vec4& center, const vec4& up);
};

#ifdef _WIN32
#define M_PI 3.14159265
#include <algorithm>
#endif


