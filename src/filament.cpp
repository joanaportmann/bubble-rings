#include "filament.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "gl.h"

using namespace std;
using namespace Eigen;

#define numberOfVerticesPerTubeCircle 30
#define _USE_MATH_DEFINES
#define RM_mu 0.4723665527
#define delta 0.6420127083
//#define circulation 4
#define gravity -9.8
#define kinematic_viscosity 1e-06
#define At -1
#define nu 1e-06
#define NUMBER_OF_POINTS 17
//=============================================================================

Filament::Filament()
{

    filamentPositions.resize(NUMBER_OF_POINTS, 3);

    filamentPositions << -5.2, 3.2, 0.0,
        -4.4, 3.94, 0.0,
        -3.7, 4.4, 0.0,
        -2.85, 4.68, 0.0,
        -1.88, 4.72, 0.0,
        -0.43, 4.62, 0.0,
        0.2, 4.18, 0.0,
        0.87, 3.68, 0.0,
        1.09, 3.24, 0.0,
        1.2, 2.7, 0.0,
        1.45, 2.08, 0.0,
        1.5, 1.32, 0.0,
        1.3, 0.2, 0.0,
        -0.6, -1.45, 0.0,
        -2.8, -1.5, 0.0,
        -4.95, -0.7, 0.0,
        -5.5, 1.6, 0.5;

    thickness.resize(NUMBER_OF_POINTS);
    thickness << 0.3, 0.4, 0.3, 0.5, 0.4, 0.4, 0.5, 0.4, 0.3, 0.4, 0.4, 0.5, 0.3, 0.2, 0.2, 0.1, 0.1;

    circulation.resize(NUMBER_OF_POINTS);
    circulation << 0.3, 0.4, 0.3, 0.5, 0.4, 0.4, 0.5, 0.4, 0.3, 0.4, 0.4, 0.5, 0.3, 0.2, 0.2, 0.1, 0.1;
}

MatrixXf Filament::getFilamentPoints()
{
    return filamentPositions;
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

    for (int i = 0; i < filamentPositions.rows(); i++)
    {
        vec3 edgeAfter = filamentPositions.col((i + 1) % filamentPositions.rows()) - filamentPositions.col(i);
        vec3 edgeBefore = filamentPositions.col(i) - filamentPositions.col((i - 1 + filamentPositions.rows()) % filamentPositions.rows());
        std::vector<vec3> verticesOfOneCircle = verticesofOneCircle_(
            numberOfVerticesPerTubeCircle,
            filamentPositions.col(i),
            (edgeBefore + edgeAfter).normalized(),
            thickness(i));
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
    return (i + filamentPositions.rows()) % filamentPositions.rows();
}

// Calculating u_LIA
vec3 Filament::localizedInduction(int j, MatrixXd temp_filamentPoints)
{
    // Grab data
    vec3 e_next = temp_filamentPoints.row(wrap(j + 1)) - temp_filamentPoints.row(j);
    vec3 e_prev = temp_filamentPoints.row(j) - temp_filamentPoints.row(wrap(j - 1));
    float l_prev = e_prev.norm();
    float l_next = e_next.norm();
    float a_prev = thickness(j);
    float a_next = thickness(wrap(j + 1));

    // Curvature
    vec3 kB = 2.0 * e_prev.normalized().cross(e_next.normalized()) / (e_prev + e_next).norm();

    // Circulation
    float C = 0.5 * (circulation(j) + circulation(wrap(j + 1)));

    // Log term
    float logTerm = log(l_prev * l_next / (a_prev * a_next * delta * delta));

    // Compute
    return C / (4 * M_PI) * 0.5 * logTerm * kB;
}

//-------------------------------------------------------------------------------------

vec3 Filament::biotSavartAndLocalizedInduction(int j, std::vector<FilamentPoint> temp_filamentPoints)
{
    vec3 temp_vel = vec3(0, 0, 0);
    vec3 position = temp_filamentPoints[j].position;

    for (int j = 1; j < temp_filamentPositions.rows(); j++)
    {
        vec3 R0 = temp_filamentPoints[j - 1].position;
        vec3 R1 = temp_filamentPoints[wrap(j + 1)].position;
        float a = temp_filamentPoints[j].a;
        float Gamma = circulation(j);
        temp_vel += biotsavartedge(position, R0, R1, Gamma, a);
    }

    temp_vel += localizedInduction(j, temp_filamentPoints);
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
vec3 Filament::boussinesq_on_edge(int i, MatrixXd temp_filamentPoints)
{
    //  Boussinesq on edges
    // Read of edge
    float a = ;
    float C = ;

    // Coefficients are defined as constants above
    vec3 g = vec3(0, 0, gravity);

    // Get points and tangents
    vec3 srcP = temp_filamentPositions.col(i);
    vec3 dstP = temp_filamentPoints[wrap(i + 1)].position;
    vec3 edge = dstP - srcP;
    MatrixXd edges;
    edges.resize(NUMBER_OF_POINTS)

    for( int j = 0; j < filamentPositions.rows(); j++)
    {
       edges.row(i) =  temp_filamentPositions.col(i) - temp_filamentPositions.col(wrap(i+1));
    }
    
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

vec3 Filament::oneStepOfRungeKutta(int i, std::vector<FilamentPoint> temp_filamentPoints)
{
    vec3 v_temp;
    float time_step_;
    time_step_ = 0.001f;

    // Calculating u_BS per vertex of filament

    v_temp = biotSavartAndLocalizedInduction(i, temp_filamentPoints);

    // Calculating and adding normal flow velocity γ_normal and averaging to vertices
    vec3 y_normal = (boussinesq_on_edge(i, temp_filamentPoints) + boussinesq_on_edge((wrap(i + 1)) % temp_filamentPositions.rows()(), temp_filamentPoints)) / 2;

    v_temp += y_normal;
    v_temp *= time_step_;
    return v_temp;
};

void Filament::updateFilament()
{
    std::vector<FilamentPoint> temp_polygon1, temp_polygon2, temp_polygon3;
    temp_polygon1 = controlPolygon_;
    std::vector<vec3> K1, K2, K3, K4;

    for (int i = 0; i < filamentPositions.rows(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, controlPolygon_);
        K1.push_back(temp_K);
        temp_polygon1[i].position += temp_K * 0.5;
    }

    temp_polygon2 = temp_polygon1;

    for (int i = 0; i < filamentPositions.rows(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon1);
        K2.push_back(temp_K);
        temp_polygon2[i].position += temp_K * 0.5;
    }

    temp_polygon3 = temp_polygon2;

    for (int i = 0; i < filamentPositions.rows(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon2);
        K3.push_back(temp_K);
        temp_polygon3[i].position += temp_K;
    }

    for (int i = 0; i < filamentPositions.rows(); i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon2);
        K4.push_back(temp_K);
    }

    for (int i = 0; i < filamentPositions.rows(); i++)
    {
        filamentPositions.col(i) += (K1[i] + 2 * K2[i] + 2 * K3[i] + K4[i]) / 6;
    }
};

void Filament::updateSkeleton()
{
    updateFilament();
};

//-----------------------------------------------------------------------------------