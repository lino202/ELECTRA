#include "IMP/vectors/vec.hpp"
#include <gtest/gtest.h>

namespace IMP_test {

class VecTest : public ::testing::Test {

protected:

  void SetUp() {}

  IMP::Vec<1, unsigned> vec_u1_;
  IMP::Vec<2, int> vec_i2_;
  IMP::Vec<3, float> vec_f3_;
  IMP::Vec<4, double> vec_d4_;

};

} // End of namespace IMP_test.
