
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
    Filament();

    std::vector<FilamentPoint> getFilamentPoints();

    std::vector<vec3> getBubbleRingSkeleton();

    float time_step_ = 0.0000000001f;

    
    
    // Todo
    void updateSkeleton();
private:

    std::vector<FilamentPoint> controlPolygon_ = {
      {{-5.2, 3.2, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{-4.4, 3.94, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{-3.7, 4.4, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{-2.85, 4.68, 0.0}, 0.4, 1, vec3(0, 0, 0)},
      {{-1.88, 4.72, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{-0.43, 4.62, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{0.2, 4.18, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{0.87, 3.68, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{1.09, 3.24, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{1.2, 2.7, 0.0}, 0.3, 1.1, vec3(0, 0, 0)},
      {{1.45, 2.08, 0.0}, 0.5, 1, vec3(0, 0, 0)},
      {{1.5, 1.32, 0.0}, 0.5, 1.1, vec3(0, 0, 0)},
      {{1.3, 0.2, 0.0}, 0.1, 0.9, vec3(0, 0, 0)},
      {{-0.6, -1.45, 0.0}, 0.5,1.1, vec3(0, 0, 0)},
      {{-2.8, -1.5, 0.0}, 0.3, 0.8, vec3(0, 0, 0)},
      {{-4.95, -0.7, 0.0}, 0.3, 0.7, vec3(0, 0, 0)},
      {{-5.5, 1.6, 0.5}, 0.1, 1.5, vec3(0, 0, 0)}};;

    int wrap(int i);


    void updateFilament(); 

    // Biotsavart velocity
    vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a);
    vec3 biotSavartAndLocalizedInduction(int j, std::vector<FilamentPoint> temp_controlPolygon_);
    vec3 localizedInduction(int j, std::vector<FilamentPoint> temp_controlPolygon_);
    vec3 boussinesq_on_edge(int i, std::vector<FilamentPoint> temp_controlPolygon_);
    vec3 oneStepOfRungeKutta(int i, std::vector<FilamentPoint> temp_controlPolygon_);

    std::vector<vec3> circleVertices_t(int n, vec3 normal);

};


//=============================================================================
#endif
//=============================================================================
