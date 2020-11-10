
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


class Filament
{

public:

    // Constructor
    Filament(std::vector<FilamentPoint> tubeVertices_);

    std::vector<FilamentPoint> getFilamentPoints();

    std::vector<vec3> getBubbleRingSkeleton();

    float time_step_ = 1.0f;

    
    
    // Todo
    void updateSkeleton();
private:

    std::vector<FilamentPoint> controlPolygon_;


    void updateFilament(); 

    // Biotsavart velocity
    static vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a);
    static vec3 biotSavartAndLocalizedInduction(int j, std::vector<FilamentPoint> temp_controlPolygon_);
    static vec3 localizedInduction(int j, std::vector<FilamentPoint> temp_controlPolygon_);
    static vec3 boussinesq_on_edge(int i, std::vector<FilamentPoint> temp_controlPolygon_);
    static vec3 oneStepOfRungeKutta(int i, std::vector<FilamentPoint> temp_controlPolygon_);

    std::vector<vec3> circleVertices_t(int n, vec3 normal);

};


//=============================================================================
#endif
//=============================================================================
