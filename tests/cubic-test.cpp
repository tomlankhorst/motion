#include "motion/profile.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace {

using MPCf = motion::profile::cubic<float>;
class MPTest : public ::testing::Test {
 protected:
    const float tol = 1e-4f;
    MPCf mp;
};

TEST_F(MPTest, InitBoundaryTest) {
  EXPECT_NEAR(0.0f, mp.q(), tol);
  EXPECT_NEAR(0.0f, mp.v(), tol);

  EXPECT_NEAR(0.0f, mp.q_at(-1), tol);
  EXPECT_NEAR(0.0f, mp.v_at(-1), tol);

  EXPECT_NEAR(0.0f, mp.q_at(1), tol);
  EXPECT_NEAR(0.0f, mp.v_at(1), tol);
};

TEST_F(MPTest, UnitProfileTest) {
  mp.set( 1, 1, 0 );
  EXPECT_NEAR(0.0f, mp.q_at(0), tol);
  EXPECT_NEAR(0.0f, mp.v_at(0), tol);

  // Due to symmetry, expect half-way at half-time
  EXPECT_NEAR(0.5f, mp.q_at(0.5), tol);

  EXPECT_NEAR(1.0f, mp.q_at(1), tol);
  EXPECT_NEAR(0.0f, mp.v_at(1), tol);

  // Expect that out of range evaluates to exactly the last in range value
  EXPECT_EQ( mp.q_at(1), mp.q_at(2) );
  EXPECT_EQ( mp.q_at(-1), mp.q_at(0) );
}

TEST_F(MPTest, UnitProfileOffsetTest) {
  mp.set( 1,2,3,4,5,6 );
  EXPECT_NEAR(3, mp.q_at(1), tol);
  EXPECT_NEAR(5, mp.v_at(1), tol);
  
  EXPECT_NEAR(3.3813, mp.q_at(1.5), 1e-1);

  // Some random value (fair diceroll)
  EXPECT_NEAR(3.3817, mp.q_at(1.1010), 1e-1);

  EXPECT_NEAR(4, mp.q_at(2), tol);
  EXPECT_NEAR(6, mp.v_at(2), tol);

  // Expect that out of range evaluates to exactly the last in range value
  EXPECT_EQ( mp.q_at(0), mp.q_at(1) );
  EXPECT_EQ( mp.q_at(2), mp.q_at(3) );
}


using MPCd = motion::profile::cubic<double>;
class MPDTest : public ::testing::Test {
protected:
    const double tol = 1e-9;
    MPCd mp;
};

TEST_F(MPDTest, UnitProfileTest) {
  mp.set( 1, 1, 0 );
  EXPECT_NEAR(0.0f, mp.q_at(0), tol);
  EXPECT_NEAR(0.0f, mp.v_at(0), tol);

  // Due to symmetry, expect half-way at half-time
  EXPECT_NEAR(0.5f, mp.q_at(0.5), tol);

  EXPECT_NEAR(1.0f, mp.q_at(1), tol);
  EXPECT_NEAR(0.0f, mp.v_at(1), tol);
}

}  // namespace
