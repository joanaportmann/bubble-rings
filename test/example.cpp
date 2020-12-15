#include "gtest/gtest.h"

TEST(example, SevenIsSeven)
{
    
    ASSERT_EQ(7, 7);
}

TEST(example, SevenIsNotSix)
{
    
    ASSERT_EQ(7, 6);
}
