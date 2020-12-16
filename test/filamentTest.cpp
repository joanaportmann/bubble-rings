#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "filament.h"


class FilamentTest : public ::testing::Test
{

protected:

    vec3 biotsavartedge(vec3 p, vec3 R0, vec3 R1, float gamma, float a)
    {
        filament.biotsavartedge(p, R0, R1, gamma, a);
    }

    Eigen::VectorXd doBurgerStepOnBubbleRing()
    {
        filament.doBurgerStepOnBubbleRing();
    }


    void SetUp() override {}

    void TearDown() override {}

private:

    Filament filament;
};


using ::testing::ElementsAre;

TEST_F(FilamentTest, doBurgerStepOnBubbleRingTest)
{
    
    Eigen::VectorXd result_burger = doBurgerStepOnBubbleRing();
    vec3 result = biotsavartedge(vec3(0, 0, 1), vec3(1, 1, 1), vec3(2, 1, 0), 2.2, 0.12);
    // ASSERT_EQ(result, vec3(0, 1, 0));
    // ASSERT_EQ(result, ElementsAre(0, 1, 0));

    // ASSERT_TRUE(result.isApprox(vec3(0, 1, 0)));

    ASSERT_EQ(result(0), 0);
    ASSERT_EQ(result(1), 1);
    ASSERT_EQ(result(2), 0);
}
