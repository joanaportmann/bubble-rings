#include "filament.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdlib.h>

using namespace std;

#define numberOfVerticesPerTubeCircle 30
#define _USE_MATH_DEFINES
#define RM_mu 0.4723665527
#define delta 0.6420127083
#define circulation 4
#define gravity -9.8
#define kinematic_viscosity 1e-06
#define At -1
#define nu 1e-06

//=============================================================================

Filament::Filament()
{

    // Set filament circle
    for (float i = 0; i <= 2 * M_PI; i += 0.5)
    {
        controlPolygon_.push_back({{cos(i), sin(i), 0},
                                   0.05,
                                   1,
                                   vec3(0, 0, 0)});
        // cout << controlPolygon_[i-1].position << "_____________________________ \n" ;
    }

    // cout << controlPolygon_.size();
}

std::vector<FilamentPoint> Filament::getFilamentPoints()
{
    return controlPolygon_;
}

//-----------------------------------------------------------------------------

std::vector<vec3> verticesofOneCircle_(int n, vec3 center, vec3 normal, float radius)
{
    std::vector<vec3> vertices;

    for (int i = 0; i < n; i++)
    {
        float x = cos(2 * M_PI / n * i) * radius;
        float y = sin(2 * M_PI / n * i) * radius;

        auto rotation = Eigen::Quaterniond::FromTwoVectors(vec3(0, 0, 1), normal);

        vec3 vertex = (rotation.matrix() * vec3(x, y, 0)) + center;

        vertices.push_back(vertex);
    }

    return vertices;
}

//----------------------------------------------------------------------------------

std::vector<vec3> Filament::getBubbleRingSkeleton()
{
    std::vector<vec3> verticesOfTube;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 edgeAfter = controlPolygon_[(i + 1) % controlPolygon_.size()].position - controlPolygon_[i].position;
        vec3 edgeBefore = controlPolygon_[i].position - controlPolygon_[(i - 1 + controlPolygon_.size()) % controlPolygon_.size()].position;
        std::vector<vec3> verticesOfOneCircle = verticesofOneCircle_(
            numberOfVerticesPerTubeCircle,
            controlPolygon_[i].position,
            (edgeBefore + edgeAfter).normalized(),
            controlPolygon_[i].a);
        for (int j = 0; j < verticesOfOneCircle.size(); j++)
        {
            verticesOfTube.push_back(verticesOfOneCircle[j]);
        }
    };

    return verticesOfTube;
};

//------------------------------------------------------------------------------------

//Biotsavart velocity

vec3 Filament::biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a)
{
    float aSqr = a * a * RM_mu * RM_mu;
    vec3 R0_ = R0 - p;
    vec3 R1_ = R1 - p;
    vec3 RPrime = R1 - R0;
    vec3 cross01 = R0_.cross(R1_);
    float r1 = R1_.dot(RPrime) / (sqrt(aSqr + R1_.norm()) * (RPrime.norm() * aSqr + cross01.norm()));
    float r0 = R0_.dot(RPrime) / (sqrt(aSqr + R0_.norm()) * (RPrime.norm() * aSqr + cross01.norm()));

    return Gamma * (r1 - r0) * cross01 / (4 * M_PI);
}

//--------------------------------------------------------------------------------------

int Filament::wrap(int i)
{
    return (i + controlPolygon_.size()) % controlPolygon_.size();
}

// Calculating u_LIA
vec3 Filament::localizedInduction(int j, std::vector<FilamentPoint> temp_controlPolygon_)
{
    // Grab data
    vec3 e_next = temp_controlPolygon_[wrap(j + 1)].position - temp_controlPolygon_[j].position;
    vec3 e_prev = temp_controlPolygon_[j].position - temp_controlPolygon_[wrap(j - 1)].position;
    float l_prev = e_prev.norm();
    float l_next = e_next.norm();
    float a_prev = temp_controlPolygon_[j].a;
    float a_next = temp_controlPolygon_[wrap(j + 1)].a;

    // Curvature
    vec3 kB = 2.0 * e_prev.normalized().cross(e_next.normalized()) / (e_prev + e_next).norm();

    // Circulation
    float C = 0.5 * (temp_controlPolygon_[j].C, temp_controlPolygon_[wrap(j + 1)].C);

    // Log term
    float logTerm = log(l_prev * l_next / (a_prev * a_next * delta * delta));

    // Compute
    return C / (4 * M_PI) * 0.5 * logTerm * kB;
}

//-------------------------------------------------------------------------------------

vec3 Filament::biotSavartAndLocalizedInduction(int j, std::vector<FilamentPoint> temp_controlPolygon_)
{
    vec3 temp_vel = vec3(0, 0, 0);
    vec3 position = temp_controlPolygon_[j].position;

    for (int j = 1; j < temp_controlPolygon_.size(); j++)
    {
        vec3 R0 = temp_controlPolygon_[j - 1].position;
        vec3 R1 = temp_controlPolygon_[wrap(j + 1)].position;
        float a = temp_controlPolygon_[j].a;
        float Gamma = circulation;
        temp_vel += biotsavartedge(position, R0, R1, Gamma, a);
    }

    temp_vel += localizedInduction(j, temp_controlPolygon_);
    return temp_vel;
}

//-------------------------------------------------------------------------------------

/** The last two terms of Eq. (18) are evaluated on each edge, where a_j and T_j are both defined  
 * 
 * "Using the fact that(−16πν+CT×)−1=(256π2ν2+C2)−1(−16πν−CT×) is in the plane orthogonal to T, 
 * we split this equation into normal and tangential differential equations for γ" Eq. (13a) (13b)
 * 
 * The tangential part (Eq. (13b)) equation does not change the shape of the curve and can 
 * be reduced to Burgers’ equation for the cross sectional area A=πa2 on a fixed curve.
 **/
vec3 Filament::boussinesq_on_edge(int i, std::vector<FilamentPoint> temp_controlPolygon_)
{
    //  Boussinesq on edges
    // Read of edge
    float a = temp_controlPolygon_[i].a;
    float C = temp_controlPolygon_[i].C;

    // Coefficients are defined as constants above
    vec3 g = vec3(0, 0, gravity);

    // Get points and tangents
    vec3 srcP = temp_controlPolygon_[i].position;
    vec3 dstP = temp_controlPolygon_[wrap(i + 1)].position;
    vec3 edge = dstP - srcP;
    vec3 T = edge.normalized();

    // Compute parameters
    float drag_t = 8 * M_PI * nu;
    float drag_n = drag_t * 2;
    vec3 Atg = At * g;
    vec3 Atg_t = Atg.dot(T) * T;
    vec3 Atg_n = Atg - Atg_t;
    float denom = drag_n * drag_n + C * C;

    //normal part
    return (
               drag_n * M_PI * a * a * Atg_n + C * M_PI * a * a * T.cross(Atg_n)) /
           denom;
};

vec3 Filament::oneStepOfRungeKutta(int i, std::vector<FilamentPoint> temp_controlPolygon_)
{
    vec3 v_temp;
    float time_step_;
    time_step_ = 1.0f;

    // Calculating u_BS per vertex of filament

    v_temp = biotSavartAndLocalizedInduction(i, temp_controlPolygon_);

    // Calculating and adding normal flow velocity γ_normal and averaging to vertices
    vec3 y_normal = (boussinesq_on_edge(i, temp_controlPolygon_) + boussinesq_on_edge((wrap(i + 1)) % temp_controlPolygon_.size(), temp_controlPolygon_)) / 2;

    v_temp += y_normal;
    v_temp *= time_step_;
    return v_temp;
};

void Filament::updateFilament()
{
    std::vector<FilamentPoint> temp_polygon1, temp_polygon2, temp_polygon3;
    temp_polygon1 = controlPolygon_;
    std::vector<vec3> K1, K2, K3, K4;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, controlPolygon_);
        K1.push_back(temp_K);
        temp_polygon1[i].position += temp_K * 0.5;
    }

    temp_polygon2 = temp_polygon1;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon1);
        K2.push_back(temp_K);
        temp_polygon2[i].position += temp_K * 0.5;
    }

    temp_polygon3 = temp_polygon2;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon2);
        K3.push_back(temp_K);
        temp_polygon3[i].position += temp_K;
    }

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon2);
        K4.push_back(temp_K);
    }

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        controlPolygon_[i].position += (K1[i] + 2 * K2[i] + 2 * K3[i] + K4[i]) / 6 * 0.002;
    }
};

void Filament::updateSkeleton()
{
    updateFilament();
};

//-----------------------------------------------------------------------------------