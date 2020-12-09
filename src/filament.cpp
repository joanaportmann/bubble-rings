#include "filament.h"
#include "glmath.h"
#include <vector>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <Eigen/Sparse>

using namespace std;

typedef Eigen::Triplet<double> T;

// ATTENTION: Keep in sync with the one in tube.cpp
#define numberOfVerticesPerTubeCircle 20
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
    for (float i = 0; i <= 2 * M_PI; i += 0.17)
    {
        controlPolygon_.push_back({{cos(i), sin(i), 0},
                                   0.12,
                                   4.2,
                                   vec3(0, 0, 0)});
    }

    size = controlPolygon_.size();
    cout << "size: " << size << endl;
}

//----------------------------------------------------------------------------

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

    for (int i = 0; i < size; i++)
    {
        vec3 edgeAfter = controlPolygon_[(i + 1) % size].position - controlPolygon_[i].position;
        vec3 edgeBefore = controlPolygon_[i].position - controlPolygon_[(i - 1 + size) % size].position;
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
    return (i + size) % size;
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

    // Circulation
    float C = 0.5 * (temp_controlPolygon_[j].C, temp_controlPolygon_[wrap(j + 1)].C);

    // Log term
    float logTerm = log(l_prev * l_next / (a_prev * a_next * delta * delta));

    // Compute
    return C / (4 * M_PI) * 0.5 * logTerm * kB;
}

//-------------------------------------------------------------------------------------

vec3 Filament::biotSavartAndLocalizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_)
{
    vec3 temp_vel = vec3(0, 0, 0);
    vec3 position = temp_controlPolygon_[j].position;

    for (int j = 0; j < size; j++)
    {
        vec3 R0 = temp_controlPolygon_[wrap(j - 1)].position;
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
vec3 Filament::boussinesq_on_edge(int i, const std::vector<FilamentPoint> &temp_controlPolygon_)
{
    //  Boussinesq on edges
    // Read of edge
    float a = temp_controlPolygon_[i].a;
    float C = circulation;

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
    vec3 v_temp;

    // Calculating u_BS per vertex of filament

    v_temp = biotSavartAndLocalizedInduction(i, temp_controlPolygon_);

    // Calculating and adding normal flow velocity γ_normal and averaging to vertices
    vec3 y_normal = (boussinesq_on_edge(i, temp_controlPolygon_) + boussinesq_on_edge((wrap(i + 1)), temp_controlPolygon_)) / 2;

    v_temp += y_normal;
    v_temp *= time_step_;
    return v_temp;
};

void Filament::updateFilament()
{
    std::vector<FilamentPoint> temp_polygon1, temp_polygon2, temp_polygon3;
    temp_polygon1 = controlPolygon_;
    std::vector<vec3> K1, K2, K3, K4;

    for (int i = 0; i < size; i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, controlPolygon_);
        K1.push_back(temp_K);
        temp_polygon1[i].position += temp_K * 0.5;
    }

    temp_polygon2 = controlPolygon_;
    for (int i = 0; i < size; i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon1);
        K2.push_back(temp_K);
        temp_polygon2[i].position += temp_K * 0.5;
    }

    temp_polygon3 = controlPolygon_;

    for (int i = 0; i < size; i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon2);
        K3.push_back(temp_K);
        temp_polygon3[i].position += temp_K;
    }

    for (int i = 0; i < size; i++)
    {
        vec3 velocity;
        vec3 temp_K;
        temp_K = Filament::oneStepOfRungeKutta(i, temp_polygon3);
        K4.push_back(temp_K);
    }

    for (int i = 0; i < size; i++)
    {
        controlPolygon_[i].position += (K1[i] + 2 * K2[i] + 2 * K3[i] + K4[i]) / 6;
    }
};

void Filament::preComputations()
{   
    edges_.clear();
    tangents_.clear();
    lengths_.clear();
    areas_.clear();
    effectiveGravities_.clear();
    point_lengths_.clear();
    flux_.clear();

    for (int i = 0; i < size; i++)
        edges_.push_back(controlPolygon_[wrap(i + 1)].position - controlPolygon_[i].position);
    for (int i = 0; i < edges_.size(); i++)
        tangents_.push_back(edges_[i].normalized());
    for (int i = 0; i < edges_.size(); i++)
        lengths_.push_back(edges_[i].norm());
    for (int i = 0; i < size; i++)
        areas_.push_back(std::pow(controlPolygon_[i].a, 2) * M_PI);
    for (int i = 0; i < size; i++)
        effectiveGravities_.push_back((vec3(0, gravity, 0) * At).dot(tangents_[i]));

    //TODO: Average a and C (edge -> point)

    for (int i = 0; i < size; i++)
    point_lengths_.push_back( (lengths_[i] + lengths_[wrap(i - 1)]) / 2);

    // compute point flux as in Godunov's method

    for (int i = 0; i < size; i++)
    {
        // compute (gravity(prevPrim) = gravity(i-1))
        float minus = effectiveGravities_[wrap(i - 1)] * areas_[wrap(i - 1)];
        float plus = effectiveGravities_[i] * areas_[i];

        if (minus > std::max(0.0f, (-plus)))
        {
            // positive case
            flux_.push_back(1.0 / (8 * M_PI) * minus * areas_[wrap(i - 1)]);
            AreaUsed_ = areas_[wrap(i - 1)];
        }
        else if (plus < std::min(0.0f, (-minus)))
        {
            // negative case
            flux_.push_back(1.0 / (8 * M_PI) * plus * areas_[i]);
            AreaUsed_ = areas_[i];
        }
        else
        {
            // neutral
            flux_.push_back(0.0f);
            AreaUsed_ = 0;
        }
    };
};

void Filament::doBurgerStepOnBubbleRing()
{

    preComputations();

    // Create A
    Eigen::VectorXd A(size);
    for (int j = 0; j < size; j++)
    {
        A(j) = areas_[j];
    }

    // Create F (flux * nu)
    Eigen::VectorXd F(size);
    for (int j = 0; j < size; j++)
    {
        F(j) = flux_[j];
    }

    //-----------------------------------------------------------------------

    // Prep for Laplacian (L = -d^T * star1 * d^T)

    // Create d
    std::vector<T> trp_d;
    for (int i = 0; i < size; i++)
    {
        trp_d.push_back(T(i, i, -1));
        trp_d.push_back(T(i, (i + 1) % size, 1));
    }

    Eigen::SparseMatrix<double> d(size, size); // default is column major
    d.setFromTriplets(trp_d.begin(), trp_d.end());

    Eigen::SparseMatrix<double> d_transpose = d.transpose();

    //cout << "d: " << d << endl;
    //cout << "d_transpose: " << d_transpose << endl;

    //d.makeCompressed(); // optional

    std::vector<T> trp_C_square_div_pointLength;
    for (int i = 0; i < size; i++)
    {
        
        double entry = std::pow(controlPolygon_[i].C, 2) / point_lengths_[i];
        trp_C_square_div_pointLength.push_back(T(i, i, entry));
    }

    Eigen::SparseMatrix<double> star1(size, size); // default is column major
    star1.setFromTriplets(trp_C_square_div_pointLength.begin(), trp_C_square_div_pointLength.end());

    

    // Build Laplacian L
    Eigen::SparseMatrix<double> L = -d.transpose() * star1 * d;

    //cout << "L: " << L << endl;

    // (Check sizes of matrices) cout << d.size() << "------------" << star1.size() << "++++++++++++++++++++++";
    //-----------------------------------------------------------------------

    // Create M (edgeLength diagonal matrix) Mass matrix (star0 in Houdini) multiplied with nu / time_step_
 std::vector<T> trp_lengths;
    for (int i = 0; i < size; i++)
    {
        trp_lengths.push_back(T(i, i, lengths_[i] * nu / time_step_));
    }
    Eigen::SparseMatrix<double> M(size, size); // default is column major
    M.setFromTriplets(trp_lengths.begin(), trp_lengths.end());
    cout << "M: " << M << endl;

    //-------------------------------------------------------------------------

    // Constants
    double coef = 1.0/(64. * M_PI * M_PI);

    //------------------------------------------------------------------------

    // Backward Euler  (Ax = b)
    Eigen::SparseMatrix<double> LHS = M - (0.5 * coef * L);
    Eigen::MatrixXd RHS = M * A + d.transpose() * F ;
   
    // SCALE DUE TO PRECISION
    double scale = 1.0/RHS.norm(); 

    
Eigen::VectorXd x(size);

// fill A = LHS and b = RHS
Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
cg.compute(LHS);
x = cg.solve(RHS * scale);
x = x / scale;
//std::cout << "#iterations:     " << cg.iterations() << std::endl;
//std::cout << "estimated error: " << cg.error()      << std::endl;

for (int i = 0; i < size; i++) 
{
    //cout << "update: --------------- + " << std::pow(x[i] / (2 * M_PI),2) << endl;
   controlPolygon_[i].a += sqrt(sqrt(std::pow(x[i] / (2 * M_PI), 2)));
   //controlPolygon_[i].a = 2;
};


};

void Filament::updateSkeleton()
{
    updateFilament();
    updatedFilament = true;

    doBurgerStepOnBubbleRing();
};

//-----------------------------------------------------------------------------------