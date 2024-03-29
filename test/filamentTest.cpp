#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "filament.h"
#include <vector>
#include <iostream>
#include <string>
#include <stdlib.h>

using namespace std;

class FilamentTest : public ::testing::Test
{

protected:
    // Setters
    void setControlPolygon(std::vector<FilamentPoint> filament_points_) { filament.controlPolygon_ = filament_points_; }
    void setCirculations(int C)
    {
        for (int i = 0; i < filament.controlPolygon_.size(); i++)
        {
            filament.controlPolygon_[i].C = C;
        }
    }
    void setLengths(std::vector<float> lengths_) { filament.lengths_e = lengths_; }
    void setPointLengths(std::vector<float> pointLengths_) { filament.point_lengths_v = pointLengths_; }
    void setAreas(std::vector<float> areas_) { filament.areas_e = areas_; }
    void setFlux(std::vector<float> flux_) { filament.flux_v = flux_; }

    // Access private functions of filament.cpp

    void preComputations()
    {
        filament.preComputations();
    }
    vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float gamma, float a)
    {
        return filament.biotsavartedge(p, R0, R1, gamma, a);
    }

    vec3 localizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_)
    {
        return filament.localizedInduction(j, temp_controlPolygon_);
    }

    vec3 biotSavartAndLocalizedInduction(int j, const std::vector<FilamentPoint> &temp_controlPolygon_)
    {
        return filament.biotSavartAndLocalizedInduction(j, temp_controlPolygon_);
    }

    Eigen::VectorXd doBurgerStepOnBubbleRing()
    {
        return filament.doBurgerStepOnBubbleRing();
    }

    vec3 boussinesq_on_edge(int i, const std::vector<FilamentPoint> &temp_controlPolygon_)
    {
        return filament.boussinesq_on_edge(i, temp_controlPolygon_);
    }

    vec3 velUpdateForEdgeI(int i, const std::vector<FilamentPoint> &temp_controlPolygon_)
    {
        return filament.velUpdateForEdgeI(i, temp_controlPolygon_);
    }

    void BiotSavartAndLocalizedInduction()
    {
        filament.BiotSavartAndLocalizedInduction();
    }

    void updateSkeleton()
    {
        filament.updateSkeleton();
    }

    int totalLengthOfControlpolygon()
    {
        return filament.totalLengthOfControlpolygon();
    }

    void resample(float length)
    {
        filament.resample(length);
    }

    void resampleCatmullRom(float length)
    {
        filament.resampleCatmullRom(length);
    }

    vec3 catMullRom(float tension_, float alpha__, float k, vec3 P0, vec3 P1, vec3 P2, vec3 P3)
    {
        return filament.generalCatmullRom(tension_, alpha__, k, P0, P1, P2, P3);
    }

    // SetUp and TearDown
    void SetUp() override {}
    void TearDown() override {}

    Filament filament = Filament(0.12, 4, 6);
    std::vector<FilamentPoint> filamentPoints_;
};

using ::testing::ElementsAre;

// Test all edges
TEST_F(FilamentTest, doBurgerStepOnBubbleRingTest)
{
    // Simulate a control polygon of size 26 (positions don't matter)
    filamentPoints_.resize(26, FilamentPoint{vec3(0, 0, 0), 0.12f, 4.0f});
    setControlPolygon(filamentPoints_);

    float lengths[] = {0.1443928, 0.14378038, 0.14596489, 0.14522029, 0.14391315, 0.1442105,
                       0.14567757, 0.14529918, 0.14378813, 0.14460708, 0.144272, 0.14492798,
                       0.14524907, 0.14507067, 0.1436004, 0.14522156, 0.14491448, 0.14476931,
                       0.14478274, 0.14644849, 0.14147678, 0.14195484, 0.14766575, 0.14749822,
                       0.14232254, 0.14339279};
    std::vector<float> lengths__(std::begin(lengths), std::end(lengths));
    setLengths(lengths__);

    float point_lengths[] = {0.14389279, 0.1440866, 0.14487264, 0.1455926, 0.14456671, 0.14406183, 0.14494404, 0.14548838, 0.14454365, 0.14419761,
                             0.14443955, 0.14459999, 0.14508852, 0.14515987, 0.14433554, 0.14441098, 0.14506802, 0.14484189,
                             0.14477602, 0.14561561, 0.14396264, 0.14171581, 0.14481029, 0.14758199,
                             0.14491038, 0.14285767};
    std::vector<float> point_lengths__(std::begin(point_lengths), std::end(point_lengths));
    setPointLengths(point_lengths__);

    float areas[] = {0.04523184, 0.04522017, 0.04521775, 0.04522768, 0.04523836, 0.04525384,
                     0.04525717, 0.04525525, 0.04524418, 0.04522046, 0.04521942, 0.0452452,
                     0.04524633, 0.0452289, 0.04522326, 0.04522177, 0.04523537, 0.04529226,
                     0.04529799, 0.04525774, 0.04522141, 0.04514195, 0.04515675, 0.04528871,
                     0.04531785, 0.04525284};
    std::vector<float> areas__(std::begin(areas), std::end(areas));
    setAreas(areas__);

    float flux[] = {7.93794286e-04, 7.93024839e-04, 7.51929183e-04, 6.55231357e-04,
                    5.28870209e-04, 3.41711711e-04, 1.78574759e-04, -1.39399795e-04,
                    -3.77572520e-04, -5.44049835e-04, -6.72839989e-04, -7.47612328e-04,
                    -7.85265816e-04, -7.96709326e-04, -7.32710992e-04, -6.32955751e-04,
                    -5.38048625e-04, -4.00285702e-04, -1.78429633e-04, 0.0,
                    9.69532266e-06, 1.95918590e-04, 3.69355985e-04, 5.13897103e-04,
                    6.35017001e-04, 7.51230458e-04};
    std::vector<float> flux__(std::begin(flux), std::end(flux));
    setFlux(flux__);

    float nu = 1e-06;
    float time_step_ = 0.01;

    Eigen::VectorXd result_burger;
    result_burger = doBurgerStepOnBubbleRing();

    double expected_result[] = {0.04800835, 0.05702037, 0.06561341, 0.07314099, 0.0791778, 0.08307021, 0.08512572, 0.08354294, 0.0792586, 0.07309467,
                                0.06545621, 0.05695891, 0.0480035, 0.03891299, 0.0306024, 0.0234225, 0.01729339, 0.01274469, 0.01072873, 0.01074678, 0.01086909, 0.01306779,
                                0.01729234, 0.02327722, 0.03053587, 0.03900016};

    std::vector<float> expected_results(std::begin(expected_result), std::end(expected_result));

    for (int i = 0; i < 26; i++)
        EXPECT_NEAR(result_burger(i), expected_results[i], 0.00000001);
}

// Test (for filament with 6 vertices) fourth edge (i = 3)
TEST_F(FilamentTest, biotSavartAndLocalizedInduction)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    filament.setModifyThicknessRungeKutta(false);

    vec3 temp_vel_0_filamentPoint;
    temp_vel_0_filamentPoint = biotSavartAndLocalizedInduction(3, filamentPoints_);

    EXPECT_NEAR(temp_vel_0_filamentPoint(0), 0.0, 0.00001);
    EXPECT_NEAR(temp_vel_0_filamentPoint(1), 0.0, 0.00001);
    EXPECT_NEAR(temp_vel_0_filamentPoint(2), 1.83756, 0.00001);
}

// Test (for filament with 6 vertices) fourth edge (i = 3)
TEST_F(FilamentTest, localizedInduction)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});

    vec3 temp_vel_0_filamentPoint;
    temp_vel_0_filamentPoint = localizedInduction(3, filamentPoints_);

    EXPECT_NEAR(temp_vel_0_filamentPoint(0), 0.0, 0.00001);
    EXPECT_NEAR(temp_vel_0_filamentPoint(1), 0.0, 0.00001);
    EXPECT_NEAR(temp_vel_0_filamentPoint(2), 1.08717, 0.00001);
}

TEST_F(FilamentTest, biotSavartEdge)
{
    vec3 result = biotsavartedge(
        vec3(1.0f, 1.0f, 1.0f),
        vec3(2.0f, 0.0f, 1.0f),
        vec3(0.9f, 0.0f, 0.8f),
        4.0f,
        0.12f);

    EXPECT_NEAR(result(0), 0.0454265, 0.000001);
    EXPECT_NEAR(result(1), 0.0454265, 0.000001);
    EXPECT_NEAR(result(2), -0.249846, 0.000001);
}

// Test (for filament with 6 vertices) third edge (i = 2)
TEST_F(FilamentTest, boussinesqOnEdge)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});

    vec3 result = boussinesq_on_edge(2, filamentPoints_);

    EXPECT_NEAR(result(0), -6.08488e-07, 0.000001);
    EXPECT_NEAR(result(1), 3.57705e-07, 0.000001);
    EXPECT_NEAR(result(2), -0.056169, 0.000001);
}

// Test (for filament with 6 vertices) third edge (i = 2)
TEST_F(FilamentTest, velUpdateForEdgeI)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    vec3 result = velUpdateForEdgeI(2, filamentPoints_);

    BiotSavartAndLocalizedInduction();

    EXPECT_NEAR(result(0), -3.04244e-09, 0.000001);
    EXPECT_NEAR(result(1), 8.75252e-09, 0.000001);
    EXPECT_NEAR(result(2), 0.0174899, 0.000001);
}

// Test Runge Kutta for first vertex
TEST_F(FilamentTest, BiotSavartAndLocalizedInduction)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);
    BiotSavartAndLocalizedInduction();

    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();

    vec3 result = controlPolygon__[3].position;

    EXPECT_NEAR(result(0), -0.586357, 0.00001);
    EXPECT_NEAR(result(1), 1.53235e-05, 0.00001);
    EXPECT_NEAR(result(2), 0.018367, 0.00001);

    // all values:
    // 0.611781	1.37976e-05	0.018357
    // 0.31931	0.519628	0.0176254
    // -0.280898	0.519631	0.0174899
    // -0.586357	1.53235e-05	0.018367
    // -0.293264	-0.5196	0.0192374
    // 0.313136	-0.519601	0.0191793
}

TEST_F(FilamentTest, a_after_RungeKuttaAndBurger)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    BiotSavartAndLocalizedInduction();
    preComputations(); // for Burger's equation
    Eigen::VectorXd x = doBurgerStepOnBubbleRing();
    Eigen::VectorXd a;
    a.resize(filamentPoints_.size());
    for (int i = 0; i < filamentPoints_.size(); i++)
        a(i) = sqrt(sqrt(std::pow(x(i) / (M_PI), 2)));

    // Checking thickness
    double expected_result[] = {0.133562, 0.168164, 0.133802, 0.0863786, 0.0864597, 0.0864255};
    std::vector<float> expected_results(std::begin(expected_result), std::end(expected_result));

    for (int i = 0; i < filamentPoints_.size(); i++)
    {
        EXPECT_NEAR(a(i), expected_results[i], 0.000001);
        cout << "thickness: " << a(i) << "\n";
    }

    // Equivalent to updateSkeleton:
    //
    //     for (int i = 0; i < 1; i++)
    //         filament.updateSkeleton();
    //     std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();

    //     double expected_result[] =
    //         {0.133562, 0.168164, 0.133802, 0.0863786, 0.0864597, 0.0864255};
    //     std::vector<float> expected_results(std::begin(expected_result), std::end(expected_result));
    //     for (int i = 0; i < 6; i++)
    //     {
    //         EXPECT_NEAR(controlPolygon__[i].a, expected_result[i], 0.01);
    //         cout << "thickness: update Filament" << controlPolygon__[i].a << "\n";
    //     }
}

TEST_F(FilamentTest, position_afterUpdateFilament)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    BiotSavartAndLocalizedInduction();
    preComputations(); // for Burger's equation
    Eigen::VectorXd x = doBurgerStepOnBubbleRing();
    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();
    vec3 result = vec3(-0.293264, -0.5196, 0.0192374);

    // All values:
    // 0.611781	1.37976e-05	0.018357
    // 0.31931	0.519628	0.0176254
    // -0.280898	0.519631	0.0174899
    // -0.586357	1.53235e-05	0.018367
    // -0.293264	-0.5196	0.0192374
    // 0.313136	-0.519601	0.0191793

    for (int i = 0; i < 3; i++)
        EXPECT_NEAR(controlPolygon__[4].position(i), result(i), 0.000001);
}

// Setup (Houdini and C++)
// Startpositions:
// 0.611779	0.0	0.0
// 0.319311	0.519615 0.0
// -0.280896 0.519615 0.0
// -0.586356 -5.24537e-08 0.0
// -0.293264 -0.519615 0.0
// 0.313135	-0.519615 0.0

// frame 156 in Houdini
// a = 0.12 (constant)
// C = 4 (constant)
// Biot-Savart field and Biot-Savart loc. ind. and normal flow
TEST_F(FilamentTest, 155_Iterations_OfRungeKutta)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);
        filament.setModifyThicknessRungeKutta(false);

    for (int i = 0; i < 155; i++)
        BiotSavartAndLocalizedInduction();
    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();
    vec3 result = vec3(0.322385, -0.178072, 2.93021);

    // All Positions in 156th frame:
    // 0.628901	0.341356	2.80164
    // 0.325663	0.859853	2.67128
    // -0.292737	0.863521	2.67299
    // -0.606598	0.343823	2.80012
    // -0.302304	-0.175722	2.9266
    // 0.322385	-0.178072	2.93021

    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[5].position(i), result(i), 0.000001);
}

TEST_F(FilamentTest, 100_Iterations_OfRungeKutta_bigControlpolygon)
{
    filamentPoints_.push_back({{0.6, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.559483, 0.216745, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.443405, 0.404217, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.267443, 0.537098, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.055361, 0.597441, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.164198, 0.577095, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.361581, 0.47881, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.51013, 0.315859, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.589784, 0.11025, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.589784, -0.11025, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.51013, -0.315859, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.361581, -0.47881, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.164198, -0.577095, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.0553609, -0.597441, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.267443, -0.537098, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.443405, -0.404218, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.559483, -0.216745, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);
        filament.setModifyThicknessRungeKutta(false);

    for (int i = 0; i < 100; i++)
        BiotSavartAndLocalizedInduction();
    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();
    vec3 result = vec3(-0.166883, 0.74118, 1.71075);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[5].position(i), result(i), 0.000001);
}

TEST_F(FilamentTest, 792_Iterations_OfRungeKutta)
{
    filamentPoints_.push_back({{0.611779, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);
        filament.setModifyThicknessRungeKutta(false);

    for (int i = 0; i < 792; i++)
        BiotSavartAndLocalizedInduction();
    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();
    vec3 result = vec3(0.476367, 5.3476, 11.4972);

    // All Positions in 156th frame:
    // 0.960341	5.86888	10.8438
    // 0.485476	6.38932	10.1874
    // -0.481976	6.39004	10.1812
    // -0.973293	5.86918	10.8351
    // -0.498115	5.34976	11.4941
    // 0.476367	5.3476	11.4972
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[5].position(i), result(i), 0.0001);
}

// Setup (Houdini and C++)
// Startpositions:
// 0.611779	0.0	0.0
// 0.319311	0.519615	0.0
// -0.280896	0.519615	0.0
// -0.586356	-5.24537e-08	0.0
// -0.293264	-0.519615	0.0
// 0.313135	-0.519615	0.0
//
// frame 241 in Houdini
// a = 0.12 (constant)
// C = 4 (constant)
// Biot-Savart field and Biot-Savart loc. ind. and normal flow and thickness flow and modify thickness
// TEST_F(FilamentTest, updateSkeleton)
// {
//     filamentPoints_.push_back({{0.611779, 0.0, 0.0},
//                                0.12,
//                                4});
//     filamentPoints_.push_back({{0.319311, 0.519615, 0.0},
//                                0.12,
//                                4});
//     filamentPoints_.push_back({{-0.280896, 0.519615, 0.0},
//                                0.12,
//                                4});
//     filamentPoints_.push_back({{-0.586356, -5.24537e-08, 0.0},
//                                0.12,
//                                4});
//     filamentPoints_.push_back({{-0.293264, -0.519615, 0.0},
//                                0.12,
//                                4});
//     filamentPoints_.push_back({{0.313135, -0.519615, 0.0},
//                                0.12,
//                                4});
//     setControlPolygon(filamentPoints_);

//     //for (int i = 0; i < 86; i++)
//     //for (int i = 0; i < 1; i++)
//     filament.updateSkeleton();
//     std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();
//     int size = controlPolygon__.size();
//     EXPECT_EQ(size, 48);
// }

TEST_F(FilamentTest, resampleSquare)
{
    filamentPoints_.push_back({{0, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{1, 0, 0},
                               0.12,
                               4});
    filamentPoints_.push_back({{1, 1, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0, 1, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    resample(0.85);

    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();

    vec3 expected_result_0 = vec3(0, 0, 0);
    vec3 expected_result_1 = vec3(0.8, 0, 0);
    vec3 expected_result_2 = vec3(1, 0.6, 0);
    vec3 expected_result_3 = vec3(0.6, 1, 0);
    vec3 expected_result_4 = vec3(0, 0.8, 0);

    EXPECT_EQ(controlPolygon__.size(), 5);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[0].position(i), expected_result_0[i], 0.0000001);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[1].position(i), expected_result_1[i], 0.0000001);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[2].position(i), expected_result_2[i], 0.0000001);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[3].position(i), expected_result_3[i], 0.0000001);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[4].position(i), expected_result_4[i], 0.0000001);
}

TEST_F(FilamentTest, resampleSquareFine)
{
    filamentPoints_.push_back({{0, 0.0, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{1, 0, 0},
                               0.12,
                               4});
    filamentPoints_.push_back({{1, 1, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0, 1, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    resample(0.128);

    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();

    vec3 expected_result_0 = vec3(0, 0, 0);
    vec3 expected_result_1 = vec3(0.129032258, 0, 0);
    vec3 expected_result_8 = vec3(1, 0.032, 0);

    EXPECT_EQ(controlPolygon__.size(), 31);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[0].position(i), expected_result_0[i], 0.0000001);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[1].position(i), expected_result_1[i], 0.000001);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[8].position(i), expected_result_8[i], 0.01);
}

TEST_F(FilamentTest, resampleRandomPolygon)
{
    filamentPoints_.push_back({{1.24156, -4.66619, 0},
                               0.12,
                               4});
    filamentPoints_.push_back({{-1.92, -2.636, 0.5},
                               0.12,
                               4});
    filamentPoints_.push_back({{-10.23, 2.11, 2},
                               0.12,
                               4});
    filamentPoints_.push_back({{1.11, 8.44, 0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    resample(0.5);

    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();

    vec3 expected_result_0 = vec3(1.24156, -4.66619, 0);
    vec3 expected_result_1 = vec3(-1.92 - 0.2222785 * 0.858145, -2.63 + 0.2222785 * 0.489483, 0.5 + 0.2 * 0.1549);

    // Polygonlength: 3.8 + 9.68 + 13.14 + 13.1 = 39.72
    EXPECT_EQ(controlPolygon__.size(), 79);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[0].position(i), expected_result_0[i], 0.0000001);
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(controlPolygon__[8].position(i), expected_result_1[i], 0.01);
}

TEST_F(FilamentTest, resampleSquareThicknessConstant)
{
    filamentPoints_.push_back({{0, 0.0, 0.0},
                               0.12,
                               1});
    filamentPoints_.push_back({{1, 0, 0},
                               0.12,
                               4});
    filamentPoints_.push_back({{1, 1, 0.0},
                               0.12,
                               4});
    filamentPoints_.push_back({{0, 1, 0.0},
                               0.12,
                               4});
    setControlPolygon(filamentPoints_);

    resample(0.85);

    std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();

    EXPECT_EQ(controlPolygon__.size(), 5);
    EXPECT_NEAR(controlPolygon__[1].C, 0.8 * 4 + 0.2 * 1, 0.0001);
    for (int i = 0; i < 5; i++)
        EXPECT_NEAR(controlPolygon__[i].a, 0.12, 0.0000001);
    for (int i = 2; i < 4; i++)
        EXPECT_NEAR(controlPolygon__[i].C, 4, 0.00001);
}

TEST_F(FilamentTest, catMullRom)
{
    vec3 result1 = catMullRom(0.0f, 0.0f, 0.0f, {0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0});
 cout << result1 << "\n";
    vec3 result2 = catMullRom(0.0f, 0.0f, 1.0f, {0, 10, 0}, {1, 20, 0}, {2, 2, 0}, {3, 0, 0});

    EXPECT_EQ(result1[0], 1);
}

// TEST_F(FilamentTest, CatmullRomResampleSquare)
// {
    // filamentPoints_.push_back({{0, 0.0, 0.0},
    //                            0.12,
    //                            4});
    // filamentPoints_.push_back({{1, 0, 0},
    //                            0.12,
    //                            4});
    // filamentPoints_.push_back({{1, 1, 0.0},
    //                            0.12,
    //                            4});
    // filamentPoints_.push_back({{0, 1, 0.0},
    //                            0.12,
    //                            4});
    // setControlPolygon(filamentPoints_);

    // resampleCatmullRom(0.85);

    // std::vector<FilamentPoint> controlPolygon__ = filament.getFilamentPoints();

    // EXPECT_EQ(controlPolygon__.size(), 8);

    // EXPECT_NEAR(controlPolygon__[3].position(0), 1.09375, 0.0001);
    // EXPECT_NEAR(controlPolygon__[3].position(1), 0.507812, 0.0001);
    // EXPECT_NEAR(controlPolygon__[3].position(2), 0, 0.0001);
// }