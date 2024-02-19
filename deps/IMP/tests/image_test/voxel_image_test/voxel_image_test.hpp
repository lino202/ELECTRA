#include "IMP/image/voxel_image.hpp"
#include "IMP/image/voxel_image_properties.hpp"
#include <gtest/gtest.h>

#include <vector>

namespace IMP_test {


struct triplet {
    int x_, y_, z_;

    void Set(int x, int y, int z) { x_ = x; y_ = y; z_ = z; }

};


class VoxelImageTest : public ::testing::Test {

protected:

  void SetUp() {
      // Number of voxels with value one.
      this->val_one_.resize(28);

      // Voxel coordinates with value one.
      this->val_one_.at(0).Set(2, 3, 6); this->val_one_.at(14).Set(2, 3, 5);
      this->val_one_.at(1).Set(3, 3, 6); this->val_one_.at(15).Set(3, 3, 5);
      this->val_one_.at(2).Set(4, 3, 6); this->val_one_.at(16).Set(4, 3, 5);
      this->val_one_.at(3).Set(5, 3, 6); this->val_one_.at(17).Set(5, 3, 5);
      this->val_one_.at(4).Set(6, 3, 6); this->val_one_.at(18).Set(6, 3, 5);

      this->val_one_.at(5).Set(3, 4, 6); this->val_one_.at(19).Set(3, 4, 5);
      this->val_one_.at(6).Set(4, 4, 6); this->val_one_.at(20).Set(4, 4, 5);
      this->val_one_.at(7).Set(5, 4, 6); this->val_one_.at(21).Set(5, 4, 5);
      this->val_one_.at(8).Set(6, 4, 6); this->val_one_.at(22).Set(6, 4, 5);

      this->val_one_.at(9).Set(4, 5, 6); this->val_one_.at(23).Set(4, 5, 5);
      this->val_one_.at(10).Set(5, 5, 6); this->val_one_.at(24).Set(5, 5, 5);
      this->val_one_.at(11).Set(6, 5, 6); this->val_one_.at(25).Set(6, 5, 5);

      this->val_one_.at(12).Set(5, 6 ,6); this->val_one_.at(26).Set(5, 6, 5);
      this->val_one_.at(13).Set(6, 6 ,6); this->val_one_.at(27).Set(6, 6, 5);

      // Number of voxels with value two.
      this->val_two_.resize(16);

      // Voxel coordinates with value two.
      this->val_two_.at(0).Set(6, 3, 4); this->val_two_.at(8).Set(6, 3, 3);
      this->val_two_.at(1).Set(7, 3, 4); this->val_two_.at(9).Set(7, 3, 3);
      this->val_two_.at(2).Set(8, 3, 4); this->val_two_.at(10).Set(8, 3, 3);

      this->val_two_.at(3).Set(6, 4, 4); this->val_two_.at(11).Set(6, 4, 3);
      this->val_two_.at(4).Set(7, 4, 4); this->val_two_.at(12).Set(7, 4, 3);

      this->val_two_.at(5).Set(6, 5, 4); this->val_two_.at(13).Set(6, 5, 3);
      this->val_two_.at(6).Set(7, 5, 4); this->val_two_.at(14).Set(7, 5, 3);

      this->val_two_.at(7).Set(6, 6, 4); this->val_two_.at(15).Set(6, 6, 3);
  }

  IMP::VoxelImage<unsigned> uint_image_;
  IMP::VoxelImage<int> int_image_;
  IMP::VoxelImage<float> float_image_;
  IMP::VoxelImage<double> double_image_;

  std::vector<IMP_test::triplet> val_one_;
  std::vector<IMP_test::triplet> val_two_;

};

} // End of namespace IMP_test.
