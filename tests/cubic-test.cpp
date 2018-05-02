#include "motion/profile.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace {

typedef motion::profile::cubic<float> MP;
class MPTest : public ::testing::Test {
 protected:
  MP mp;
};

TEST_F(MPTest, InitBoundaryTest) {
  EXPECT_EQ(0.0f, mp.q());
  EXPECT_EQ(0.0f, mp.v());

  EXPECT_EQ(0.0f, mp.q_at(-1));
  EXPECT_EQ(0.0f, mp.v_at(-1));

  EXPECT_EQ(0.0f, mp.q_at(1));
  EXPECT_EQ(0.0f, mp.v_at(1));
};

}  // namespace
