#include "testunit.h"

FooTest::FooTest() {
}

FooTest::~FooTest() {};

void FooTest::SetUp() {};

void FooTest::TearDown() {};

TEST_F(FooTest, SevenIsSeven) {
    EXPECT_EQ(7, 7);
}

TEST_F(FooTest, SixIsEight) {
    EXPECT_EQ(6, 8);
}
