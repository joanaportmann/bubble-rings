#pragma once

#include <iostream>
#include <cassert>
#include <math.h>
#include <Eigen/Core>

using vec4 = Eigen::Vector4d;
using vec3 = Eigen::Vector3d;
using mat3 = Eigen::Matrix3d;
using mat4 = Eigen::Matrix4d;

#ifdef _WIN32
#define M_PI 3.14159265
#include <algorithm>
#endif


