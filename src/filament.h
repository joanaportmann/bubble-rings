
#ifndef FILAMENT_H
#define FILAMENT_H
//=============================================================================

#include "gl.h"
#include <Eigen/Dense>
#include "glmath.h"
#include <vector>
//=============================================================================


struct FilamentPoint
{
    vec3 position;
    float a;
    float C;
};

class Filament
{

public:

    // Constructor
    Filament();

    ~Filament();

    std::vector<FilamentPoint> getFilamentPoints();

    std::vector<vec3> getBubbleRingSkeleton();

    float time_step_ = 0.01f;
    bool updatedFilament = true;

    // Start configuration of filament
    float circulation = 4;
    float thickness = 0.12;
    
    // Todo
    void updateSkeleton();
    

    

    friend class FilamentTest;

private:

   Eigen::VectorXd doBurgerStepOnBubbleRing();

   // Thickness flow: Burger's equation
    
    std::vector<vec3> edges_e;
    std::vector<vec3> tangents_e;
    std::vector<float> lengths_e;
    std::vector<float> point_lengths_v;
    std::vector<float> areas_e;
    std::vector<float> effectiveGravities_e;
    std::vector<float> flux_v;
    float AreaUsed_v;
    // Eigen::SparseMatrix<double> L_matrix;
    std::vector<FilamentPoint> controlPolygon_ ;
    int size;

    

    int wrap(int i);


    void BiotSavartAndLocalizedInduction(); 

    // Biotsavart velocity
    vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a);
    vec3 biotSavartAndLocalizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 localizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 boussinesq_on_edge(int i, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 oneStepOfRungeKutta(int i, const std::vector<FilamentPoint> &temp_controlPolygon_);

    std::vector<vec3> circleVertices_t(int n, vec3 normal);



    void preComputations ();
    
   

};


//=============================================================================
#endif
//=============================================================================
