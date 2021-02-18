#include "filament.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include "CatmullRom.h"
#include <assert.h>

using namespace std;

typedef Eigen::Triplet<double> T;

// ATTENTION: Keep in sync with the one in tube.cpp
#define numberOfVerticesPerTubeCircle 30

#define _USE_MATH_DEFINES
#define RM_mu 0.4723665527f
#define delta 0.6420127083
#define circulation 4
#define gravity -9.8
#define kinematic_viscosity 1e-06
#define At -1
#define nu 1e-06


//=============================================================================

// Create Filament of radius 0.6 and add random float [0, 0.02] to x of vertices. Outcommit random adding for testing.
Filament::Filament(float thickness_, float circulation_, int numEdges)
{
    // int unsigned time_seed;
    // time_seed = static_cast<unsigned>(time(0));
    // srand(time_seed);
    // cout << "Time seed " << time_seed << "\n";
    // srand(static_cast<unsigned>(time(0)));
    unsigned int seed = 1612214691;

    // Seed from Houdini
    //unsigned int seed = 1612215069;
    srand(seed);
    // Set filament circle
    float angle = 2 * M_PI / numEdges;
    for (int i = 0; i < numEdges; i ++)
    {
        float r0 = 0;
        float r1 = 0;
        float r2 = 0;
        // r0 = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / 0.02));
        // r1 = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / 0.02));
        // r2 = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / 0.02));

        controlPolygon_.push_back({{0.6 * cos(i * angle) + r0, 0.6 * sin(i * angle) + r1, 0.0 + r2},
                                   thickness_,
                                   circulation_});

        vec3 position = {0.6 * cos(i * angle) + r0, 0.6 * sin(i * angle) + r1, 0.0 + r2};
        cout << "startposition: " << position << "\n" << "\n";
    }
    originalControlPolygon_ = controlPolygon_;
}

//---------------------------------------------------------------------------

Filament::~Filament() {}

//----------------------------------------------------------------------------

std::vector<FilamentPoint> Filament::getFilamentPoints()
{
    return controlPolygon_;
}

//-----------------------------------------------------------------------------

std::vector<vec3> Filament::verticesofOneCircle_(int n, vec3 center, vec3 normal, vec3 up, float radius, bool reCenter)
{
    Eigen::Matrix3d rotation;
    double angle;
    if (normal(0) == 0 && normal(1) == 0)
    {
        vec3 normalized_normal = normal.normalized();
        if (normalized_normal(2) == -1)
        {
            vec3 axis = vec3(0, 0, 1);
            rotation = Eigen::AngleAxisd(M_PI, axis);
        }
        else
        {
            rotation = Eigen::Matrix3d::Identity();
        }
    }
    else
    {
        angle = acos(normal(2) / normal.norm());
        vec3 axis = vec3(-normal(1), normal(0), 0);
        rotation = Eigen::AngleAxisd(angle, axis.normalized());
    }

    // Rotation of first vertex of ring in space
    vec3 vertexNotOriented = rotation * vec3(radius, 0, 0);

    // Determine direction of rotation (since acos is always positive)
    int signUpCrossVertex = up.cross(vertexNotOriented).dot(normal) < 0 ? 1 : -1;

    double angleToUp = acos(vertexNotOriented.dot(up) / (vertexNotOriented.norm() * up.norm()));
    Eigen::Matrix3d rotationToOrientRing;
    rotationToOrientRing = Eigen::AngleAxisd(signUpCrossVertex * angleToUp, normal.normalized());

    std::vector<vec3> vertices;
    for (int i = 0; i < n; i++)
    {
        float x = cos(2 * M_PI / n * i) * radius;
        float y = sin(2 * M_PI / n * i) * radius;
        vec3 vertex = rotationToOrientRing * rotation * vec3(x, y, 0);

        vertex += center;
        if (recenter)
        {

            // TODO!!!!!!!!!!!!!!!
            vec3 shift = originalControlPolygon_[0].position - controlPolygon_[0].position;
            vertex += shift;
        }
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
        vec3 v1 = controlPolygon_[0].position - controlPolygon_[controlPolygon_.size() * 2 / 3].position;
        vec3 v2 = controlPolygon_[0].position - controlPolygon_[controlPolygon_.size() * 1 / 3].position;
        vec3 up = v1.cross(v2);
        std::vector<vec3> verticesOfOneCircle = verticesofOneCircle_(
            numberOfVerticesPerTubeCircle,
            controlPolygon_[i].position,
            (edgeBefore + edgeAfter).normalized(),
            up,
            controlPolygon_[i].a,
            recenter);
        for (int j = 0; j < verticesOfOneCircle.size(); j++)
        {
            verticesOfTube.push_back(verticesOfOneCircle[j]);
        }
    };

    return verticesOfTube;
};

//------------------------------------------------------------------------------------

// Biotsavart velocity
// biotsavartedge
// evaluates the Biot-Savart field at a point p generated by an edge with end points R0 and R1 with circulation Gamma.
// a is the core radius, and RM_mu encodes the cutoff model (Saffman P.213).

vec3 Filament::biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a)
{
    float aSqr = a * a * RM_mu * RM_mu;
    vec3 R0_ = R0 - p;
    vec3 R1_ = R1 - p;
    vec3 RPrime = R1 - R0;
    vec3 cross01 = R0_.cross(R1_);
    float r1 = R1_.dot(RPrime) / (sqrt(aSqr + R1_.norm() * R1_.norm()) * (RPrime.norm() * RPrime.norm() * aSqr + cross01.norm() * cross01.norm()));
    float r0 = R0_.dot(RPrime) / (sqrt(aSqr + R0_.norm() * R0_.norm()) * (RPrime.norm() * RPrime.norm() * aSqr + cross01.norm() * cross01.norm()));

    return Gamma * (r1 - r0) * cross01 / (4 * M_PI);
}

//--------------------------------------------------------------------------------------

int Filament::wrap(int i)
{
    return (i + controlPolygon_.size()) % controlPolygon_.size();
}

// Calculating u_LIA
vec3 Filament::localizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_)
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

    // cout << "kB: " << kB << "\n";
    // Circulation
    float C = 0.5 * (temp_controlPolygon_[j].C + temp_controlPolygon_[wrap(j + 1)].C);

    // Log term
    float logTerm = log(l_prev * l_next / (a_prev * a_next * delta * delta));

    //cout << "logTerm: " << logTerm << "\n";

    // Compute
    //cout << "localizedInduction: " << C / (4 * M_PI) * 0.5 * logTerm * kB << "\n";
    return C / (4 * M_PI) * 0.5 * logTerm * kB;
}

//-------------------------------------------------------------------------------------

vec3 Filament::biotSavartAndLocalizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_)
{
    vec3 temp_vel = vec3(0, 0, 0);
    vec3 position = temp_controlPolygon_[j].position;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 R0 = temp_controlPolygon_[wrap(i)].position;
        vec3 R1 = temp_controlPolygon_[wrap(i + 1)].position;
        float a = temp_controlPolygon_[i].a;
        float Gamma = temp_controlPolygon_[i].C;
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
vec3 Filament::boussinesq_on_edge(int i, const std::vector<FilamentPoint> &temp_controlPolygon_)
{
    //  Boussinesq on edges
    // Read of edge
    float a = temp_controlPolygon_[i].a;
    float C = temp_controlPolygon_[i].C;

    // Coefficients are defined as constants above
    vec3 g = vec3(0, gravity, 0);

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

vec3 Filament::oneStepOfRungeKutta(int i, const std::vector<FilamentPoint> &temp_controlPolygon_)
{
    vec3 v_temp = vec3(0, 0, 0);

    // Calculating u_BS per vertex of filament

    v_temp = biotSavartAndLocalizedInduction(i, temp_controlPolygon_);

    // Calculating and adding normal flow velocity γ_normal and averaging to vertices
    vec3 y_normal = (boussinesq_on_edge(wrap(i - 1), temp_controlPolygon_) + boussinesq_on_edge(i, temp_controlPolygon_)) / 2;

    v_temp += y_normal;
    v_temp *= time_step_;
    return v_temp;
};

// Runge Kutta evaluation
void Filament::BiotSavartAndLocalizedInduction()
{
    std::vector<FilamentPoint> temp_polygon1, temp_polygon2, temp_polygon3;
    temp_polygon1 = controlPolygon_;
    std::vector<vec3> K1, K2, K3, K4;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, controlPolygon_);
        K1.push_back(temp_K);
        temp_polygon1[i].position += temp_K * 0.5;
    }

    temp_polygon2 = controlPolygon_;
    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon1);
        K2.push_back(temp_K);
        temp_polygon2[i].position += temp_K * 0.5;
    }

    temp_polygon3 = controlPolygon_;

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon2);
        K3.push_back(temp_K);
        temp_polygon3[i].position += temp_K;
    }

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon3);
        K4.push_back(temp_K);
    }

    for (int k = 0; k < controlPolygon_.size(); k++)
    {
        vec3 update = (K1[k] + 2 * K2[k] + 2 * K3[k] + K4[k]) / 6;
        // cout << "update BiotSavart " << update << "\n"
        //      << "\n";
        controlPolygon_[k].position += update;
    }
};

/* Precomputations for Burger's equation. 
* _e: Per edge (in Houdini)
* _v: Per vertes (in Houdini)
*/
void Filament::preComputations()
{
    // edge computations
    edges_e.clear();
    tangents_e.clear();
    lengths_e.clear();
    areas_e.clear();
    effectiveGravities_e.clear();
    point_lengths_v.clear();
    flux_v.clear();

    for (int i = 0; i < controlPolygon_.size(); i++)
        edges_e.push_back(controlPolygon_[wrap(i + 1)].position - controlPolygon_[i].position);
    for (int i = 0; i < edges_e.size(); i++)
        tangents_e.push_back(edges_e[i].normalized());
    for (int i = 0; i < edges_e.size(); i++)
        lengths_e.push_back(edges_e[i].norm());
    for (int i = 0; i < controlPolygon_.size(); i++)
        areas_e.push_back(std::pow(controlPolygon_[i].a, 2) * M_PI);
    for (int i = 0; i < controlPolygon_.size(); i++)
        effectiveGravities_e.push_back((vec3(0, gravity, 0) * At).dot(tangents_e[i]));

    //TODO: Average a  and C (just for stats in Houdini?)

    for (int i = 0; i < controlPolygon_.size(); i++)
        point_lengths_v.push_back((lengths_e[i] + lengths_e[wrap(i - 1)]) / 2);

    // compute point flux as in Godunov's method

    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        /* compute (gravity(prevPrim) = gravity(i-1))
        * minus -> area and gravity from previous vertex
        * plus -> area and gravity from next edge
        */

        float minus = effectiveGravities_e[wrap(i - 1)] * areas_e[wrap(i - 1)];
        float plus = effectiveGravities_e[i] * areas_e[i];

        if (minus > std::max(0.0f, -plus))
        {
            // positive case
            flux_v.push_back(1.0 / (8 * M_PI) * minus * areas_e[wrap(i - 1)]);
             AreaUsed_v = areas_e[wrap(i - 1)];
        }
        else if (plus < std::min(0.0f, (-minus)))
        {
            // negative case
            flux_v.push_back(1.0 / (8 * M_PI) * plus * areas_e[i]);
            AreaUsed_v = areas_e[i];
        }
        else
        {
            // neutral
            flux_v.push_back(0.0f);
            AreaUsed_v = 0;
        }
    };
};

Eigen::VectorXd Filament::doBurgerStepOnBubbleRing()
{
    // Create A
    Eigen::VectorXd A(controlPolygon_.size());
    for (int j = 0; j < controlPolygon_.size(); j++)
    {
        A(j) = areas_e[j];
    }

    // Create F (flux * nu)
    Eigen::VectorXd F(controlPolygon_.size());
    for (int j = 0; j < controlPolygon_.size(); j++)
    {
        F(j) = flux_v[j];
    }

    //-----------------------------------------------------------------------

    // Prep for Laplacian (L = -d^T * star1 * d^T)

    // Create d
    std::vector<T> trp_d;
    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        trp_d.push_back(T(i, i, 1));
        trp_d.push_back(T(i, (i - 1 + controlPolygon_.size()) % controlPolygon_.size(), -1));
    }

    Eigen::SparseMatrix<double> d(controlPolygon_.size(), controlPolygon_.size()); // default is column major
    d.setFromTriplets(trp_d.begin(), trp_d.end());

    // d.coeffRef(controlPolygon_.size()-1, 0) = 0;
    // d.coeffRef(0, controlPolygon_.size()-1) = -1;

    //cout << "d: -.-.-.-.-.-." << d << "\n" << endl;

    Eigen::SparseMatrix<double> d_transpose = d.transpose();

    //d.makeCompressed(); // optional

    std::vector<T> trp_C_square_div_pointLength;
    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        double entry = std::pow(controlPolygon_[i].C, 2) / point_lengths_v[i];
        trp_C_square_div_pointLength.push_back(T(i, i, entry));
    }

    Eigen::SparseMatrix<double> star1(controlPolygon_.size(), controlPolygon_.size()); // default is column major
    star1.setFromTriplets(trp_C_square_div_pointLength.begin(), trp_C_square_div_pointLength.end());

    //cout << "star1: -.-.-.-.-.-." << star1 << "\n" << endl;

    // Build Laplacian L
    Eigen::SparseMatrix<double> L = -d.transpose() * star1 * d;

    //cout << "L: -.-.-.-.-.-." << L << "\n" << endl;

    //-----------------------------------------------------------------------

    // Create M (edgeLength diagonal matrix) Mass matrix (star0 in Houdini) multiplied with nu / time_step_
    std::vector<T> trp_lengths;
    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        trp_lengths.push_back(T(i, i, lengths_e[i]));
    }
    Eigen::SparseMatrix<double> M(controlPolygon_.size(), controlPolygon_.size()); // default is column major
    M.setFromTriplets(trp_lengths.begin(), trp_lengths.end());

    //-------------------------------------------------------------------------

    // Constants
    double coef = 1.0 / (64. * M_PI * M_PI);
    double nuIdt = nu / time_step_;

    //------------------------------------------------------------------------

    // Backward Euler  (Ax = b)
    Eigen::SparseMatrix<double> LHS = nuIdt * M - (0.5 * coef * L);
    Eigen::MatrixXd RHS = nuIdt * M * A + d.transpose() * F;

    // cout << "LHS " << LHS << "\n" << "\n";
    // cout << "RHS " << RHS << "\n" << "\n";
    // SCALE DUE TO PRECISION
    double scale = 1.0 / RHS.norm();

    // fill A = LHS and b = RHS
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
    cg.setTolerance(1e-5);
    cg.compute(LHS);
    // cout << "estimated error: " << cg.error() << "\n";
    // cout << "tolerance: " << cg.tolerance() << "\n";
    return (cg.solve(RHS * scale)) / scale;
}

float Filament::totalLengthOfControlpolygon()
{
    float length = 0;
    for (int i = 0; i < controlPolygon_.size(); i++)
    {
        length += (controlPolygon_[i].position - controlPolygon_[wrap(i + 1)].position).norm();
    }

    return length;
}

void Filament::resample(float resampleLength)
{
    float segmentnumber_ = std::round(totalLengthOfControlpolygon() / resampleLength);
    int segmentnumber = segmentnumber_;
    float actualResampleLength = totalLengthOfControlpolygon() / segmentnumber;
    std::vector<FilamentPoint> newPoints;

    newPoints.push_back(controlPolygon_[0]);
    for (int i = 1; i < segmentnumber; i++)
    {
        float diffFromStart = i * actualResampleLength;

        float segmentLength = (controlPolygon_[0].position - controlPolygon_[1].position).norm();
        int baseVertex = 0;

        while (diffFromStart > segmentLength)
        {
            diffFromStart -= segmentLength;
            baseVertex++;
            segmentLength = (controlPolygon_[baseVertex].position - controlPolygon_[wrap(baseVertex + 1)].position).norm();
        }
        vec3 v = controlPolygon_[wrap(baseVertex + 1)].position - controlPolygon_[baseVertex].position;
        vec3 position = controlPolygon_[baseVertex].position + v.normalized() * diffFromStart;
        float weight = 1 - (diffFromStart / v.norm());
        float a = (sqrt(controlPolygon_[baseVertex].a) * weight + sqrt(controlPolygon_[wrap(baseVertex + 1)].a) * (1 - weight));
        a *= a;
        float C = controlPolygon_[baseVertex].C * weight + controlPolygon_[wrap(baseVertex + 1)].C * (1 - weight);

        newPoints.push_back({{position(0), position(1), position(2)}, a, C});
    }

    controlPolygon_.assign(newPoints.begin(), newPoints.end());
}

// uniform Catmull-Rom Splines
vec3 Filament::uniformCatmullRom(float u, vec3 &P0, vec3 &P1, vec3 &P2, vec3 &P3)
{
    vec3 point;
    point = u * u * u * ((-1) * P0 + 3 * P1 - 3 * P2 + P3) / 2;
    point += u * u * (2 * P0 - 5 * P1 + 4 * P2 - P3) / 2;
    point += u * ((-1) * P0 + P2) / 2;
    point += P1;

    return point;
}

// Catmull-Rom Splines (alpha = 0 for uniform; alpha = 0.5 for centripetal; alpha = 1 for chordal)
//
// Centripetal will not form loop or self-intersection within a curve segment,
// and cusp will never occur within a curve segment
// and it follows the control points more tightly than other Catmull-Rom Splines
//
// For tension value 0 is a good choice. These values can range from 0 to 1, where 1 means linear interpolation.
vec3 Filament::generalCatmullRom(float tension_, float alpha__, float t, vec3 &P0, vec3 &P1, vec3 &P2, vec3 &P3)
{
    float t01 = pow(vec3(P0 - P1).norm(), alpha__);
    float t12 = pow(vec3(P1 - P2).norm(), alpha__);
    float t23 = pow(vec3(P2 - P3).norm(), alpha__);

    vec3 m1 = (1.0f - tension_) *
              (P2 - P1 + t12 * ((P1 - P0) / t01 - (P2 - P0) / (t01 + t12)));
    vec3 m2 = (1.0f - tension_) *
              (P2 - P1 + t12 * ((P3 - P2) / t23 - (P3 - P1) / (t12 + t23)));

    vec3 a = 2.0f * (P1 - P2) + m1 + m2;
    vec3 b = -3.0f * (P1 - P2) - m1 - m1 - m2;
    vec3 c = m1;
    vec3 d = P1;

    return a * t * t * t +
           b * t * t +
           c * t +
           d;
};

void Filament::resampleCatMullRomWithWeight(float resampleLength)
{
    float segmentnumber_ = std::round(totalLengthOfControlpolygon() / resampleLength);
    int segmentnumber = segmentnumber_;

    float actualResampleLength = totalLengthOfControlpolygon() / segmentnumber;
    std::vector<FilamentPoint> newPoints;

    newPoints.push_back(controlPolygon_[0]);
    for (int i = 1; i < segmentnumber; i++)
    {
        float diffFromStart = i * actualResampleLength;

        float segmentLength = (controlPolygon_[0].position - controlPolygon_[1].position).norm();
        int baseVertex = 0;

        while (diffFromStart > segmentLength)
        {
            diffFromStart -= segmentLength;
            baseVertex++;
            segmentLength = (controlPolygon_[baseVertex].position - controlPolygon_[wrap(baseVertex + 1)].position).norm();
        }
        vec3 v = controlPolygon_[wrap(baseVertex + 1)].position - controlPolygon_[baseVertex].position;
        float weight = (diffFromStart / v.norm());
        //vec3 position = controlPolygon_[baseVertex].position + v.normalized() * diffFromStart;

        // uniform Catmull-Rom Splines
        // vec3 position = uniformCatmullRom(weight,
        //                                   controlPolygon_[wrap(baseVertex - 1)].position,
        //                                   controlPolygon_[baseVertex].position,
        //                                   controlPolygon_[wrap(baseVertex + 1)].position,
        //                                   controlPolygon_[wrap(baseVertex + 2)].position);

        // Centripedal Catmull-Rom Spline (alpha = 0.5)
        vec3 position = generalCatmullRom(tension,
                                          alpha,
                                          weight,
                                          controlPolygon_[wrap(baseVertex - 1)].position,
                                          controlPolygon_[baseVertex].position,
                                          controlPolygon_[wrap(baseVertex + 1)].position,
                                          controlPolygon_[wrap(baseVertex + 2)].position);

        float a = (sqrt(controlPolygon_[baseVertex].a) * (1 - weight) + sqrt(controlPolygon_[wrap(baseVertex + 1)].a) * weight);
        a *= a;
        float C = controlPolygon_[baseVertex].C * (1 - weight) + controlPolygon_[wrap(baseVertex + 1)].C * weight;

        newPoints.push_back({{position(0), position(1), position(2)}, a, C});
    }

    controlPolygon_.assign(newPoints.begin(), newPoints.end());
}

float Filament::totalVolume()
{
    float volume = 0;
    for(int i = 0; i < controlPolygon_.size(); i++) 
    {
        volume += (controlPolygon_[i].position - controlPolygon_[wrap(i+1)].position).norm() * controlPolygon_[i].a * controlPolygon_[i].a * M_PI;
    }
    return volume; 
}

void Filament::updateSkeleton()
{
    BiotSavartAndLocalizedInduction();
    if (modifyThickness)
    {
        preComputations(); // for Burger's equation
        Eigen::VectorXd x = doBurgerStepOnBubbleRing();
        for (int i = 0; i < controlPolygon_.size(); i++)
        {
            cout << x.minCoeff() << "\n" "\n"; 
            controlPolygon_[i].a = sqrt(sqrt(std::pow(x(i) / (M_PI), 2)));
        }
    }
    //resample(resampleLength_);
    resampleCatMullRomWithWeight(resampleLength_);
    updatedFilament = true;
    framecouter++;

    cout << "Volume: " << totalVolume() << "\n";
};

//-----------------------------------------------------------------------------------