#include "glmath.h"

//=============================================================================

#include <iostream>
#include <cassert>
#include <math.h>
#include <Eigen/Dense> 
//#include <Eigen/src/Geometry.h>


mat4 MatUtils::look_at(const vec4& eye_4d, const vec4& center_4d, const vec4& up_4d)
{
    vec3 eye, center, up;
    eye = eye_4d.head<3>();
    center = center_4d.head<3>();
    up = up_4d.head<3>();
    vec3 z = (eye-center).normalized();
    vec3 x = (up.cross(z)).normalized();
    vec3 y = (z.cross(x)).normalized();

    mat4 m;
    m(0,0)=x[0]; m(0,1)=x[1]; m(0,2)=x[2]; m(0,3)=-x.dot(eye);
    m(1,0)=y[0]; m(1,1)=y[1]; m(1,2)=y[2]; m(1,3)=-y.dot(eye);
    m(2,0)=z[0]; m(2,1)=z[1]; m(2,2)=z[2]; m(2,3)=-z.dot(eye);
    m(3,0)=0.0;  m(3,1)=0.0;  m(3,2)=0.0;  m(3,3)=1.0;

    return m;
}

//-----------------------------------------------------------------------------

mat4 MatUtils::frustum(float l, float r, float b, float t, float n, float f)
{
    mat4 m = mat4::Zero();
    m(0,0) = 2.0f * n / (r-l);
    m(1,1) = 2.0f * n / (t-b);
    m(0,2) = (r+l) / (r-l);
    m(1,2) = (t+b) / (t-b);
    m(2,2) = (n-f) / (f-n);
    m(3,2) = -1.0f;
    m(2,3) = -2.0f * f * n / (f-n);

    return m;
}

//-----------------------------------------------------------------------------


mat4 MatUtils::perspective(float fovy, float aspect, float near, float far)
{
    float t = near * tan( fovy * (float)M_PI / 360.0f );
    float b = -t;
    float l = b * aspect;
    float r = t * aspect;

    return MatUtils::frustum(l, r, b, t, near, far);
}


