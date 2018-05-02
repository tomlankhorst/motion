#include "motion/profile.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

namespace {

    using MPLf = motion::profile::linear<float>;

    class MPLTest : public ::testing::Test {
    protected:
        const float tol = 1e-5f;
        MPLf mp;
    };

    TEST_F(MPLTest, InitBoundaryTest) {
      EXPECT_NEAR(0.0f, mp.q(), tol);

      EXPECT_NEAR(0.0f, mp.q_at(-1), tol);

      EXPECT_NEAR(0.0f, mp.q_at(1), tol);
    };

    TEST_F(MPLTest, LinearProfileTest) {
      mp.set(2,1);
      EXPECT_NEAR(0.0f, mp.q_at(0), tol);

      // Due to symmetry, expect half-way at half-time
      EXPECT_NEAR(0.25f, mp.q_at(0.5), tol);

      EXPECT_NEAR(0.5f, mp.q_at(1), tol);

      // Expect that out of range evaluates to exactly the last in range value
      EXPECT_EQ(mp.q_at(2), mp.q_at(3));
      EXPECT_EQ(mp.q_at(-1), mp.q_at(0));
    }

} // namespace
