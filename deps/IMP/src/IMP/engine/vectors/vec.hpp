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


/**
   \file vec.hpp
   \brief Vec class header file.
   \author Konstantinos A. Mountris
   \date 14/01/2018
*/


#ifndef IMP_ENGINE_VECTORS_VEC_HPP_
#define IMP_ENGINE_VECTORS_VEC_HPP_

#include "IMP/engine/utilities/algorithms.hpp"
#include "IMP/engine/utilities/logger.hpp"

#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <utility>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <exception>
#include <stdexcept>
#include <initializer_list>
#include <functional>


namespace IMP {

/** \addtogroup Vectors \{ */


/**
 * \class Vec
 * \brief Class implenting a N-dimensional vector.
 */

template <int DIM, typename DATATYPE>
class Vec {

private:

    std::vector<DATATYPE> data_;                                  /**< The data of the Vec vector */

    const std::size_t dim_ = static_cast<std::size_t>(DIM);       /**< The number of dimensions of the Vec vector */


public:

    /**
     * \brief Vec vector default constructor.
     */
    inline Vec();


    /**
     * \overload
     * \brief Vec copy consturctor.
     * \param [in] vec The Vec vector to be copied.
     */
    inline Vec(const Vec<DIM, DATATYPE> &vec);


    /**
     * \overload
     * \brief Vec move consturctor.
     * \param [in] vec The Vec vector to be moved.
     */
    inline Vec(Vec<DIM, DATATYPE> &&vec) noexcept;


    /**
     * \overload
     * \brief Vec consturctor from std::initializer_list.
     *
     * Uses IMP::Vec<DIM, DATATYPE>::Set(std::initializer_list<DATATYPE> val_list) internally.
     *
     * \param [in] val_list The data std::initializer_list.
     */
    inline Vec(std::initializer_list<DATATYPE> val_list);


    /**
     * \overload
     * \brief Vec consturctor from copying std::vector.
     *
     * Uses IMP::Vec<DIM, DATATYPE>::Set(const std::vector<DATATYPE> &values) internally.
     *
     * \param [in] values The data std::vector to be copied.
     */
    inline Vec(const std::vector<DATATYPE> &values);


    /**
     * \overload
     * \brief Vec consturctor from moving std::vector.
     *
     * Uses IMP::Vec<DIM, DATATYPE>::Set(std::vector<DATATYPE> &&values) internally.
     *
     * \param [in] values The data std::vector to be copied.
     */
    inline Vec(std::vector<DATATYPE> &&values) noexcept;


    /**
     * \overload
     * \brief Vec consturctor from dynamic memory array.
     *
     * Uses IMP::Vec<DIM, DATATYPE>::Set(DATATYPE *data_ptr, const std::size_t &ptr_len) internally.
     *
     * \param [in] data_ptr The dynamic memory array to set the vector's data.
     * \param [in] ptr_len The length of the dynamic memory array.
     */
    inline Vec(DATATYPE *data_ptr, std::size_t ptr_len);


    /**
     * \brief Vec vector default destructor.
     */
    inline virtual ~Vec();


    /**
     * \brief Set vectors data from std::initializer_list.
     * \param [in] val_list The data std::initializer_list.
     * \return [void]
     */
    inline void Set(std::initializer_list<DATATYPE> val_list);


    /**
     * \overload
     * \brief Set vectors data from std::vector by coping.
     * \param [in] values The data std::vector to be copied.
     * \return [void]
     */
    inline void Set(const std::vector<DATATYPE> &values);


    /**
     * \overload
     * \brief Set vectors data from std::vector by moving.
     * \param [in] values The data std::vector to be moved.
     * \return [void]
     */
    inline void Set(std::vector<DATATYPE> &&values);


    /**
     * \overload
     * \brief Set vector's data from dynamic memory array.
     * \param [in] data_ptr The dynamic memory array to set the vector's data.
     * \param [in] ptr_lent The length of the dynamic memory array.
     * \return [void]
     */
    inline void Set(DATATYPE *data_ptr, std::size_t ptr_len);


    /**
     * \brief Component-wise multiplication with another IMP::Vec<DIM, DATATYPE> vector.
     * \param [in] vec The vector to be multiplied with in component-wise manner.
     * \return [IMP::Vec<DIM, DATATYPE>] The component-wise multiplication product.
     */
    inline Vec<DIM, DATATYPE> CwiseMul(const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Component-wise division  with another IMP::Vec<DIM, DATATYPE> vector.
     * \param [in] vec The vector to be divided by in component-wise manner.
     * \return [IMP::Vec<DIM, DATATYPE>] The coefficient-wise division product.
     * \note It is up to the user to ensure that division by zero does not occur.
     */
    inline Vec<DIM, DATATYPE> CwiseDiv(const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Convert the components of the vector to absolute values.
     * \return [IMP::Vec<DIM, DATATYPE>] The component-wise absolute vector.
     */
    inline Vec<DIM, DATATYPE> CwiseAbs() const;


    /**
     * \brief Get the sign for each component of the vector.
     * \return [IMP::Vec<DIM, DATATYPE>] The vector of signs of the original vector's components.
     */
    inline Vec<DIM, DATATYPE> Sign() const;


    /**
     * \brief Compute the squared length of the vector.
     */
    inline DATATYPE Length2() const;


    /**
     * \brief Compute the Euclidean norm of the vector.
     */
    inline DATATYPE Norm() const;


    inline DATATYPE Dot(const Vec<DIM, DATATYPE> &vec) const;


    inline Vec<DIM, DATATYPE> Cross(const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Compute the squared distance with another IMP::Vec<DIM, DATATYPE> vector.
     * \param [in] vec The vector to compute its distance from this vector.
     */
    inline DATATYPE Distance2(const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Compute the sum of the stored data in the IMP::Vec<DIM, DATATYPE> vector.
     * \return DATATYPE The sum of the vector's stored data.
     */
    inline DATATYPE Sum() const;


    /**
     * \brief Compute the power of the stored data in the IMP::Vec<DIM, DATATYPE> vector.
     * \param [in] exp The exponent to raise the vector.
     * \return [IMP::Vec<DIM, double>] The initial vector raised in the given exponent.
     */
    inline Vec<DIM, double> Pow(double exp);


    /**
     * \brief Get the value of the maximum component of the vector.
     * \return [DATATYPE] The value of the maximum component of the vector.
     */
    inline DATATYPE MaxCoeff() const;


    /**
     * \brief Get the value of the minimum component of the vector.
     * \return [DATATYPE] The value of the minimum componennt of the vector.
     */
    inline DATATYPE MinCoeff() const;


    /**
     * \brief Set all the components of the vector to zero.
     * \return [void]
     */
    inline void SetZero();


    /**
     * \brief Check if all the components of vector are zero.
     * \return [bool] TRUE if all the components are zero, FALSE otherwise.
     */
    inline bool IsZero() const;


    /**
     * \brief Copy the contents of the Vec vector to a dynamic Eigen::VectorXd vector.
     * \return [Eigen::Matrix<DATATYPE,Eigen::Dynamic,1>] The dynamic Eigen::VectorXd vector with the copied contents of the Vec vector.
     */
    inline Eigen::Matrix<DATATYPE,Eigen::Dynamic,1> CopyToEigen() const;


    /**
     * \brief Get the number of the vector's dimensions.
     * \return [std::vector<DATATYPE>] The number of the vector's dimensions.
     */
    inline const std::size_t & Dim() const { return this->dim_; }


    /**
     * \brief Get the data of the Vec.
     * \return [std::vector<DATATYPE>&] The data of the Vec.
     */
    inline auto & Data() const { return this->data_; }


    /**
     * \brief Get write-permission iterator pointing to the first component of the vector's data.
     * \return [std::vector<DATATYPE>::iterator] The write-permission iterator pointing to the
     *                                           first component of the vector's data.
     */
    inline typename std::vector<DATATYPE>::iterator begin() { return this->data_.begin(); }


    /**
     * \brief Get read-only iterator pointing to the first component of the vector's data.
     * \return [std::vector<DATATYPE>::iterator] The read-only iterator pointing to the
     *                                           first component of the vector's data.
     */
    inline typename std::vector<DATATYPE>::const_iterator begin() const { return this->data_.begin(); }


    /**
     * \brief Get read-only iterator pointing to the first component of the vector's data.
     * \return [std::vector<DATATYPE>::iterator] The read-only iterator pointing to the
     *                                           first component of the vector's data.
     */
    inline typename std::vector<DATATYPE>::const_iterator cbegin() const { return this->data_.cbegin(); }


    /**
     * \brief Get write-permission iterator pointing to the last component of the vector's data.
     * \return [std::vector<DATATYPE>::iterator] The write-permission iterator pointing to the
     *                                           last component of the vector's data.
     */
    inline typename std::vector<DATATYPE>::iterator end() { return this->data_.end(); }


    /**
     * \brief Get read-only iterator pointing to the last component of the vector's data.
     * \return [std::vector<DATATYPE>::iterator] The read-only iterator pointing to the
     *                                           last component of the vector's data.
     */
    inline typename std::vector<DATATYPE>::const_iterator end() const { return this->data_.cend(); }


    /**
     * \brief Get read-only iterator pointing to the last component of the vector's data.
     * \return [std::vector<DATATYPE>::iterator] The read-only iterator pointing to the
     *                                           last component of the vector's data.
     */
    inline typename std::vector<DATATYPE>::const_iterator cend() const { return this->data_.cend(); }


    /**
     * \brief Get with write-access permision the i-th component of the vector. No range check is performed.
     * \param [in] i The index of the requested component of the vector.
     * \return [DATATYPE] The i-th component of the vector vector.
     */
    inline DATATYPE & operator [] (std::size_t i);


    /**
     * \brief Get with read-only permision the i-th component of the vector. No range check is performed.
     * \param [in] i The index of the requested component of the vector.
     * \return [DATATYPE] The i-th component of the vector vector.
     */
    inline const DATATYPE & operator [] (std::size_t i) const;


    /**
     * \brief Get with write-access permision the i-th component of the vector. Range check is performed.
     * \param [in] i The index of the requested component of the vector.
     * \return [DATATYPE] The i-th component of the vector vector.
     */
    inline DATATYPE & At(std::size_t i);


    /**
     * \brief Get with read-only permision the i-th component of the vector. Range check is performed.
     * \param [in] i The index of the requested component of the vector.
     * \return [DATATYPE] The i-th component of the vector vector.
     */
    inline const DATATYPE & At(std::size_t i) const;


    /**
     * \brief Equal operator.
     * \param [in] vec The vector to be compared for equality.
     * \return [bool] True if this and vec vectors are equal.
     */
    inline bool operator == (const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Not equal operator.
     * \param [in] vec The vector to be compared for not equality.
     * \return [bool] True if this and vec vectors are not equal.
     */
    inline bool operator != (const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Less than operator.
     * \param [in] vec The vector to be compared.
     * \return [bool] True if all the components of this are less than the components of vec vector.
     */
    inline bool operator < (const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Greater than operator.
     * \param [in] vec The vector to be compared.
     * \return [bool] True if all the components of this are greater than the components of vec vector.
     */
    inline bool operator > (const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Less or equal operator.
     * \param [in] vec The vector to be compared.
     * \return [bool] True True if all the components of this are less or equal than the components of vec vector.
     */
    inline bool operator <= (const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Greater or equal operator.
     * \param [in] vec The vector to be compared.
     * \return [bool] True True if all the components of this are greater or equal than the components of vec vector.
     */
    inline bool operator >= (const Vec<DIM, DATATYPE> &vec) const;


    /**
     * \brief Copy-assignment operator.
     * \param [in] vec The vector to be copy-assigned.
     * \return [IMP::Vec<DIM, DATATYPE>] This containing the data of vec vector.
     */
    inline Vec<DIM, DATATYPE> & operator = (const Vec<DIM, DATATYPE> &vec);


    /**
     * \brief Move-assignment operator.
     * \param [in] vec The vector to be move-assigned.
     * \return [IMP::Vec<DIM, DATATYPE>] This containing the data of vec vector.
     */
    inline Vec<DIM, DATATYPE> & operator = (Vec<DIM, DATATYPE> &&vec) noexcept;


    /**
     * \brief Addition-assignment operator.
     * \param [in] vec The vector to be added to this.
     * \return [IMP::Vec<DIM, DATATYPE>] This after the addition of the vec vector.
     */
    inline Vec<DIM, DATATYPE> & operator += (const Vec<DIM, DATATYPE> &vec);


    /**
     * \brief Substraction-assignment operator.
     * \param [in] vec The vector to be substructed from this.
     * \return [IMP::Vec<DIM, DATATYPE>] This after the substraction of the vec vector.
     */
    inline Vec<DIM, DATATYPE> & operator -= (const Vec<DIM, DATATYPE> &vec);


    /**
     * \brief Convert Vec to its negative.
     * \return [Vec<DIM, DATATYPE>] The negative vector.
     */
    inline Vec<DIM, DATATYPE> operator - () const;


    /**
     * \brief Scalar multiplication-assignment operator.
     * \param [in] scalar The scalar to multiply the vector coefficients.
     * \return [IMP::Vec<DIM, DATATYPE>] This after the scalar multiplication.
     */
    inline Vec<DIM, DATATYPE> operator *= (const DATATYPE & scalar);


    /**
     * \brief Scalar division-assignment operator.
     * \param [in] scalar The scalar to divide the vector coefficients.
     * \return [IMP::Vec<DIM, DATATYPE>] This after the scalar division.
     */
    inline Vec<DIM, DATATYPE> operator /= (const DATATYPE &scalar);


    /**
     * \brief Addition operator for IMP::Vec<DIM, DATATYPE> objects.
     * \param [in] vec1 First vector to add.
     * \param [in] vec1 Second vector to add.
     * \return[IMP::Vec<DIM, DATATYPE>] The sum of the vec1 and vec2 vectors.
     */
    inline friend Vec<DIM, DATATYPE> operator + (Vec<DIM, DATATYPE> vec1, const Vec<DIM, DATATYPE> &vec2)
    {
        return vec1 += vec2;
    }


    /**
     * \brief Subtraction operator for IMP::Vec<DIM, DATATYPE> objects.
     * \param [in] vec1 First vector to subtract.
     * \param [in] vec1 Second vector to subtract.
     * \return[IMP::Vec<DIM, DATATYPE>] The difference of the vec1 and vec2 vectors.
     */
    inline friend Vec<DIM, DATATYPE> operator - (Vec<DIM, DATATYPE> vec1, const Vec<DIM, DATATYPE> &vec2)
    {
        return vec1 -= vec2;
    }


    /**
     * \brief Scalar multiplication operator for IMP::Vec<DIM, DATATYPE> object and DATATYPE scalar.
     * \param [in] vec Vector to be multiplied.
     * \param [in] scalar Scalar number to multiply the vector.
     * \return[IMP::Vec<DIM, DATATYPE>] The product of the vector and scalar multiplication.
     */
    inline friend Vec<DIM, DATATYPE> operator * (Vec<DIM, DATATYPE> vec, const DATATYPE &scalar)
    {
        return vec *= scalar;
    }


    /**
     * \brief Scalar multiplication operator for DATATYPE scalar and IMP::Vec<DIM, DATATYPE> object.
     * \param [in] scalar Scalar number to multiply the vector.
     * \param [in] vec Vector to be multiplied.
     * \return[IMP::Vec<DIM, DATATYPE>] The product of the scalar and vector multiplication.
     */
    inline friend Vec<DIM, DATATYPE> operator * (const DATATYPE &scalar, Vec<DIM, DATATYPE> vec)
    {
        return vec *= scalar;
    }


    /**
     * \brief Scalar division operator for IMP::Vec<DIM, DATATYPE> object by DATATYPE scalar.
     * \param [in] vec Vector to be divided.
     * \param [in] scalar The divisor scalar number.
     * \return[IMP::Vec<DIM, DATATYPE>] The quotient of the vector division by scalar.
     */
    inline friend Vec<DIM, DATATYPE> operator / (Vec<DIM, DATATYPE> vec, const DATATYPE &scalar)
    {
        return vec /= scalar;
    }


    /**
       \brief Stream output for IMP::Vec<DIM, DATATYPE> vector.
       \param [out] out The output stream.
       \param [in] vec The vector to be outputted.
       \return [std::ostream&] The output of the vector's data components.
    */
    inline friend std::ostream & operator << (std::ostream &out, const Vec<DIM, DATATYPE> &vec)
    {
        for (const auto &val : vec) { out << val << " "; }
        return out;
    }


    /**
     * \brief Stream input for IMP::Vec<DIM, DATATYPE> vector.
     * \param [in] is The input stream.
     * \param [out] vec The vec to load the input stream.
     * \return [std::istream&] The loaded input stream.
     */
    inline friend std::istream & operator >> (std::istream &is, Vec<DIM, DATATYPE> &vec)
    {
        for (auto &val : vec) { is >> val; }
        return is;
    }

};


/** @} End of Doxygen Groups*/


} //namespace IMP


#include "IMP/engine/vectors/vec.tpp"

#endif //IMP_ENGINE_VECTORS_VEC_HPP_
