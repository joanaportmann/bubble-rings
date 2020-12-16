#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "filament.h"
#include <vector>

class FilamentTest : public ::testing::Test
{

protected:

    // Setters
    void setSize(int size) { filament.size = size; }
    void setCirculations(int C)
    {
        for (int i = 0; i < filament.size; i++)
        {
            filament.controlPolygon_[i].C = C;
        }
    }
    void setLengths(std::vector<float> lengths_) {filament.lengths_e = lengths_;}
    void setPointLengths(std::vector<float> pointLengths_) {filament.point_lengths_v = pointLengths_;}
    void setAreas(std::vector<float> areas_) {filament.areas_e = areas_;}
    void setFlux(std::vector<float> flux_) {filament.flux_v = flux_;}

    // Access private functions of filament.cpp
    vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float gamma, float a)
    {
        filament.biotsavartedge(p, R0, R1, gamma, a);
    }

    Eigen::VectorXd doBurgerStepOnBubbleRing() 
    {
        filament.doBurgerStepOnBubbleRing();
    }

    // SetUp and TearDown
    void SetUp() override {}
    void TearDown() override {}

    Filament filament;
};

using ::testing::ElementsAre;

TEST_F(FilamentTest, doBurgerStepOnBubbleRingTest)
{
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

    // ASSERT_EQ(result, vec3(0, 1, 0));
    // ASSERT_EQ(result, ElementsAre(0, 1, 0));

    // ASSERT_TRUE(result.isApprox(vec3(0, 1, 0)));

    ASSERT_EQ(result_burger(0), 0.04800835);
}
