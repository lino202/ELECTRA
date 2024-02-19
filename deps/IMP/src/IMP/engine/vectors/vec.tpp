/*
 * IMP. Image and Mesh Processing library.
 * Copyright (C) 2016  <Konstantinos A. Mountris> <konstantinos.mountris@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef IMP_ENGINE_VECTORS_VEC_TPP_
#define IMP_ENGINE_VECTORS_VEC_TPP_


#include "IMP/engine/vectors/vec.hpp"


namespace IMP {

template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::Vec()
    : data_(DIM)
{}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::Vec(const Vec<DIM, DATATYPE> &vec)
    : data_(vec.data_)
{}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::Vec(Vec<DIM, DATATYPE> &&vec) noexcept
   : data_(std::move(vec.data_))
{}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::Vec(std::initializer_list<DATATYPE> val_list)
   : data_(DIM)
{
    this->Set(val_list);
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::Vec(const std::vector<DATATYPE> &values)
   : data_(DIM)
{
    this->Set(values);
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::Vec(std::vector<DATATYPE> &&values) noexcept
   : data_(DIM)
{
    this->Set(values);
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::Vec(DATATYPE *data_ptr, std::size_t ptr_len)
   : data_(DIM)
{
    this->Set(data_ptr, ptr_len);
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE>::~Vec()
{}


template <int DIM, typename DATATYPE>
void Vec<DIM, DATATYPE>::Set(std::initializer_list<DATATYPE> val_list)
{
    // Check for size consistency.
    if (val_list.size() != this->dim_) {
        throw std::invalid_argument(Logger::Error("Could not set vector from std::initializer_list."
                                                  " Array size not compatible with vector dimensions."));
    }

    // Assign values to vector data.
    int id = 0;
    for(const auto &val : val_list)  this->data_[id++] = val;

}


template <int DIM, typename DATATYPE>
void Vec<DIM, DATATYPE>::Set(const std::vector<DATATYPE> &values)
{
    // Check for size consistency.
    if (values.size() != this->dim_) {
        throw std::invalid_argument(Logger::Error("Could not set vector from std::vector container."
                                                  " Container size not compatible with Vec vector dimensions.").c_str());
    }

    // Assign values to vector data.
    this->data_ = values;
}


template <int DIM, typename DATATYPE>
void Vec<DIM, DATATYPE>::Set(std::vector<DATATYPE> &&values)
{
    // Check for size consistency.
    if (values.size() != this->dim_) {
        throw std::invalid_argument(Logger::Error("Could not set vector from std::vector container."
                                                  " Container size not compatible with Vec vector dimensions.").c_str());
    }

    // Move values to vector data.
    this->data_ = std::move(values);

}


template <int DIM, typename DATATYPE>
void Vec<DIM, DATATYPE>::Set(DATATYPE *data_ptr, std::size_t ptr_len)
{
    // Check for size consistency.
    if (ptr_len != this->dim_) {
        throw std::invalid_argument(Logger::Error("Could not set vector from dynamic memory array."
                                                  " Array size not compatible with vector dimensions.").c_str());
    }

    // Set the new values to the data.
    this->data_.assign(data_ptr, data_ptr+ptr_len);
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::CwiseMul(const Vec<DIM, DATATYPE> &vec) const
{
    // Component-wise multiplication result.
    Vec<DIM, DATATYPE> res;
    std::transform(this->data_.begin(), this->data_.end(), vec.data_.begin(),
                   res.data_.begin(), std::multiplies<DATATYPE>());
    return res;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::CwiseDiv(const Vec<DIM, DATATYPE> &vec) const
{
    // Component-wise division result.
    Vec<DIM, DATATYPE> res;
    std::transform(this->data_.begin(), this->data_.end(), vec.data_.begin(),
                   res.data_.begin(), std::divides<DATATYPE>());
    return res;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::CwiseAbs() const
{
    // Component-wise absolute value result.
    Vec<DIM, DATATYPE> res;
    std::transform(this->data_.begin(), this->data_.end(), res.data_.begin(), static_cast<DATATYPE (*)(DATATYPE)>(&std::abs));

    return res;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::Sign() const
{
    // Compute the sign vector.
    Vec<DIM, DATATYPE> sign;
    for (int i = 0; i != DIM; ++i) {
        sign[i] = static_cast<DATATYPE>(ALGORITHMS::Sign(this->data_[i]));
    }

    return sign;
}


template <int DIM, typename DATATYPE>
DATATYPE Vec<DIM, DATATYPE>::Length2() const
{
    // Compute squared length of the vector.
    DATATYPE sq_length = static_cast<DATATYPE>(0);
    for (const auto &val : this->data_)  sq_length += val*val;
    return sq_length;
}


template <int DIM, typename DATATYPE>
DATATYPE Vec<DIM, DATATYPE>::Norm() const
{
    // Compute the Euclidean norm of the vector.
    return std::sqrt(this->Length2());
}


template <int DIM, typename DATATYPE>
DATATYPE Vec<DIM, DATATYPE>::Dot(const Vec<DIM, DATATYPE> &vec) const
{
    auto val = DATATYPE{0};
    for (int d = 0; d != DIM; ++d) {
        val += this->data_[d]*vec[d];
    }
    return val;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::Cross(const Vec<DIM, DATATYPE> &vec) const
{
    auto cross = Vec<DIM, DATATYPE>{};
    cross[0] = ALGORITHMS::ProdsDiff(this->data_[1], vec[2], this->data_[2], vec[1]);
    cross[1] = ALGORITHMS::ProdsDiff(this->data_[2], vec[0], this->data_[0], vec[2]);
    cross[2] = ALGORITHMS::ProdsDiff(this->data_[0], vec[1], this->data_[1], vec[0]);
    return cross;
}


template <int DIM, typename DATATYPE>
DATATYPE Vec<DIM, DATATYPE>::Distance2(const Vec<DIM, DATATYPE> &vec) const
{
    // Initialize squared distance of the vectors.
    DATATYPE sq_dist = 0;

    // Compute squared distance.
    for (int id = 0; id != DIM; ++id) {
        sq_dist += (this->data_[id] - vec.data_[id]) * (this->data_[id] - vec.data_[id]);
    }
    
    return sq_dist;

}


template <int DIM, typename DATATYPE>
DATATYPE Vec<DIM, DATATYPE>::Sum() const
{
    return std::accumulate(this->data_.begin(), this->data_.end(), static_cast<DATATYPE>(0));
}


template <int DIM, typename DATATYPE>
Vec<DIM, double> Vec<DIM, DATATYPE>::Pow(double exp)
{
    Vec<DIM, double> result;

    for (auto &el : this->data_) {
        auto i = std::distance(&this->data_[0], &el);

        result[i] = std::pow(el, exp);
    }
    
    return result;
}

template <int DIM, typename DATATYPE>
DATATYPE Vec<DIM, DATATYPE>::MaxCoeff() const
{
    return *std::max_element(this->data_.begin(), this->data_.end());
}


template <int DIM, typename DATATYPE>
DATATYPE Vec<DIM, DATATYPE>::MinCoeff() const
{
    return *std::min_element(this->data_.begin(), this->data_.end());
}


template <int DIM, typename DATATYPE>
void Vec<DIM, DATATYPE>::SetZero()
{
    this->data_.clear();
    this->data_.assign(this->dim_, static_cast<DATATYPE>(0));
}


template <int DIM, typename DATATYPE>
bool Vec<DIM, DATATYPE>::IsZero() const
{
    bool is_zero = true;
    for (const auto &val : this->data_) {
        // Find non zero value.
        if (val != static_cast<DATATYPE>(0)) { is_zero = false; break; }
    }

    return is_zero;
}


template <int DIM, typename DATATYPE>
inline Eigen::Matrix<DATATYPE,Eigen::Dynamic,1> Vec<DIM, DATATYPE>::CopyToEigen() const
{
    Eigen::Matrix<DATATYPE,Eigen::Dynamic,1> eigen_vec(DIM);
    
    std::size_t i = 0;
    for (const auto &val : this->data_) eigen_vec.coeffRef(i++) = val;

    return eigen_vec;
}


template <int DIM, typename DATATYPE>
DATATYPE & Vec<DIM, DATATYPE>::operator [] (std::size_t i)
{
    return this->data_[i];
}


template <int DIM, typename DATATYPE>
const DATATYPE & Vec<DIM, DATATYPE>::operator [] (std::size_t i) const
{
    return this->data_[i];
}


template <int DIM, typename DATATYPE>
DATATYPE & Vec<DIM, DATATYPE>::At(std::size_t i)
{
    return this->data_.at(i);
}


template <int DIM, typename DATATYPE>
const DATATYPE & Vec<DIM, DATATYPE>::At(std::size_t i) const
{
    return this->data_.at(i);
}


template <int DIM, typename DATATYPE>
bool Vec<DIM, DATATYPE>::operator == (const Vec<DIM, DATATYPE> &vec) const
{
    return this->data_ == vec.data_;
}


template <int DIM, typename DATATYPE>
bool Vec<DIM, DATATYPE>::operator != (const Vec<DIM, DATATYPE> &vec) const
{
    return !(*this == vec);
}


template <int DIM, typename DATATYPE>
bool Vec<DIM, DATATYPE>::operator < (const Vec<DIM, DATATYPE> &vec) const
{
    bool comp = true;
    for (const auto &val : this->data_) {
        auto id = &val - &this->data_[0];

        // Element-wise comparison.
        if (val >= vec.data_[id]) { comp = false; break; }
    }

    return comp;
}


template <int DIM, typename DATATYPE>
bool Vec<DIM, DATATYPE>::operator > (const Vec<DIM, DATATYPE> &vec) const
{
    return vec < *this;
}


template <int DIM, typename DATATYPE>
bool Vec<DIM, DATATYPE>::operator <= (const Vec<DIM, DATATYPE> &vec) const
{
    return *this < vec || *this == vec;
}



template <int DIM, typename DATATYPE>
bool Vec<DIM, DATATYPE>::operator >= (const Vec<DIM, DATATYPE> &vec) const
{
    return *this > vec || *this == vec;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> & Vec<DIM, DATATYPE>::operator = (const Vec<DIM, DATATYPE> &vec)
{
    if (*this != vec) {
        this->data_ = vec.data_;
    }
    return *this;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> & Vec<DIM, DATATYPE>::operator = (Vec<DIM, DATATYPE> &&vec) noexcept
{
    if (*this != vec) {
        this->data_ = std::move(vec.data_);
    }
    return *this;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> & Vec<DIM, DATATYPE>::operator += (const Vec<DIM, DATATYPE> &vec)
{
    std::transform(this->data_.begin(), this->data_.end(), vec.data_.begin(),
                   this->data_.begin(), std::plus<DATATYPE>());
    return *this;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> & Vec<DIM, DATATYPE>::operator -= (const Vec<DIM, DATATYPE> &vec)
{
    std::transform(this->data_.begin(), this->data_.end(), vec.data_.begin(),
                   this->data_.begin(), std::minus<DATATYPE>());
    return *this;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::operator - () const
{
    Vec<DIM, DATATYPE> negative = *this;
    std::transform(negative.data_.begin(), negative.data_.end(), negative.data_.begin(),
                   std::bind1st(std::multiplies<DATATYPE>(), static_cast<DATATYPE>(-1)));
    return negative;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::operator *= (const DATATYPE &scalar)
{
    std::for_each(this->data_.begin(), this->data_.end(), [scalar](DATATYPE &val){ val *= scalar; });
    return *this;
}


template <int DIM, typename DATATYPE>
Vec<DIM, DATATYPE> Vec<DIM, DATATYPE>::operator /= (const DATATYPE &scalar)
{
    // Check for zero-division.
    if (scalar == static_cast<DATATYPE>(0)) {
        throw std::overflow_error(Logger::Error("Division by zero not allowed for IMP::Vec<DIM, DATATYPE>").c_str());
    }

    std::for_each(this->data_.begin(), this->data_.end(), [scalar](DATATYPE &val){ val /= scalar; });
    return *this;
}


} // End of namespace IMP


#endif //IMP_ENGINE_VECTORS_VEC_TPP_
