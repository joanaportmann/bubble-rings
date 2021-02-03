
#ifndef FILAMENT_H
#define FILAMENT_H
//=============================================================================

#include "gl.h"
#include <Eigen/Dense>
#include "glmath.h"
#include <vector>
#include "CatmullRom.h"
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
    Filament(float thickness, float circulation);

    // Deconstructor
    ~Filament();

    // Variables
    std::vector<FilamentPoint> getFilamentPoints();
    std::vector<vec3> getBubbleRingSkeleton();
    float time_step_;
    float resampleLength_;
    bool updatedFilament = true;

    // Methods
    void updateSkeleton();

    // Setters
    void setTension(float t) { tension = t; }
    void setAlpha(float alpha_) { alpha = alpha_; }
    void setResampleLength(float length) { resampleLength_ = length; }

    friend class FilamentTest;

private:
    // Variables
    std::vector<vec3> edges_e;
    std::vector<vec3> tangents_e;
    std::vector<float> lengths_e;
    std::vector<float> point_lengths_v;
    std::vector<float> areas_e;
    std::vector<float> effectiveGravities_e;
    std::vector<float> flux_v;
    float AreaUsed_v;
    std::vector<FilamentPoint> controlPolygon_;
    std::vector<vec3> circleVertices_t(int n, vec3 normal);
    CatmullRom curve;
    float tension;
    float alpha;
   
    // Methods

    int wrap(int i);
    float totalLengthOfControlpolygon();
    void resample(float resampleLength);
    vec3 uniformCatmullRom(float u,  vec3 &P0,  vec3 &P1,  vec3 &P2,  vec3 &P3);
    vec3 generalCatmullRom(float tension, float alpha, float u, vec3 &P0, vec3 &P1, vec3 &P2, vec3 &P3);
    void resampleCatMullRomWithWeight(float resampleLength);
    void resampleCatmullRom(float resampleL);

    // Biotsavart velocity
    vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float Gamma, float a);
    vec3 biotSavartAndLocalizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 localizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 boussinesq_on_edge(int i, const std::vector<FilamentPoint> &temp_controlPolygon_);
    vec3 oneStepOfRungeKutta(int i, const std::vector<FilamentPoint> &temp_controlPolygon_);
    void BiotSavartAndLocalizedInduction();

    // Thickness flow: Burger's equation
    void preComputations();
    Eigen::VectorXd doBurgerStepOnBubbleRing();
};

//=============================================================================
#endif
//=============================================================================
