
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
    vec3 K;
};

// struct Triplet
// {
//   int  one_, two_, three_;
// };


class Filament
{

public:

    // Constructor
    Filament();

    std::vector<FilamentPoint> getFilamentPoints();

    std::vector<vec3> getBubbleRingSkeleton();

    float time_step_ = 0.01f;
    bool updatedFilament = true;
    
    // Todo
    void updateSkeleton();
    
private:

    std::vector<FilamentPoint> controlPolygon_ ;
    int size;

    int wrap(int i);


    void updateFilament(); 

    // Biotsavart velocity
    vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a);
    vec3 biotSavartAndLocalizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 localizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 boussinesq_on_edge(int i, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 oneStepOfRungeKutta(int i, const std::vector<FilamentPoint> &temp_controlPolygon_);

    std::vector<vec3> circleVertices_t(int n, vec3 normal);


    // Thickness flow: Burger's equation
    
    std::vector<vec3> edges_e;
    std::vector<vec3> tangents_e;
    std::vector<float> lengths_e;
    std::vector<float> point_lengths_v;
    std::vector<float> areas_e;
    std::vector<float> effectiveGravities_e;
    std::vector<float> flux_v;
    float AreaUsed_v;

    void preComputations ();
    
    void doBurgerStepOnBubbleRing();

};


//=============================================================================
#endif
//=============================================================================
