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
   \file predicates.hpp
   \brief Predicate functions header file.
   \author Konstantinos A. Mountris
   \date 10/12/2020
*/

#ifndef IMP_ENGINE_UTILITIES_PREDICATES_HPP_
#define IMP_ENGINE_UTILITIES_PREDICATES_HPP_

#include <cmath>		
#ifdef CPU86
#include <float.h>
#endif /* CPU86 */
#ifdef LINUX
#include <fpu_control.h>
#endif /* LINUX */

namespace IMP {


/**
 * \namespace PREDICATES
 * \brief Predicate functions for adaptive precision operations used in computational geometry programs \cite{shewchuk1997adaptive}.
 */
namespace PREDICATES {

/** \addtogroup Utilities \{ */

static double splitter;     /**< Used to split floats in half, splitter = 2^ceiling(p / 2) + 1 */

static double epsilon;      /**< Used to estimate roundoff errors, epsilon = 2^(-p) */

/* A set of coefficients used to calculate maximum roundoff errors.          */
static double resulterrbound;

static double ccwerrboundA, ccwerrboundB, ccwerrboundC;

static double o3derrboundA, o3derrboundB, o3derrboundC;

static double iccerrboundA, iccerrboundB, iccerrboundC;

static double isperrboundA, isperrboundB, isperrboundC;

// Options to choose types of geometric computtaions. 
// Added by H. Si, 2012-08-23.
static int  _use_inexact_arith; // -X option.
static int  _use_static_filter; // Default option, disable it by -X1

// Static filters for orient3d() and insphere(). 
// They are pre-calcualted and set in exactinit().
// Added by H. Si, 2012-08-23.
static double o3dstaticfilter;
static double ispstaticfilter;


/**
 * \fn FastTwoSumTail(double a, double b, double &x, double &y)
 * \brief 
 * \tparam double 
 * \param [in] a 
 * \param [in] b 
 * \param [in] x 
 * \param [in] y 
 */
inline void FastTwoSumTail(double a, double b, double &x, double &y) 
{
    double bvirt = x - a;
    y = b - bvirt;
}


/**
 * \fn FastTwoSum(double a, double b, double &x, double &y)
 * \brief
 * \tparam double 
 * \param [in] a 
 * \param [in] b 
 * \param [in] x 
 * \param [in] y
 * \return [void] 
 */
inline void FastTwoSum(double a, double b, double &x, double &y) 
{
    x = a + b;
    FastTwoSumTail(a, b, x, y);
}


inline void FastTwoDiffTail(double a, double b, double &x, double &y)
{
    double bvirt = a - x;
    y = bvirt - b;
}


inline void FastTwoDiff(double a, double b, double &x, double &y)
{
    x = a - b; 
    FastTwoDiffTail(a, b, x, y);
}


inline void TwoSumTail(double a, double b, double &x, double &y)
{
    double bvirt = x - a;
    double avirt = x - bvirt;
    double bround = b - bvirt; 
    double around = a - avirt;
    y = around + bround;
}


inline void TwoSum(double a, double b, double &x, double &y) 
{
    x = a + b; 
    TwoSumTail(a, b, x, y);
}


inline void TwoDiffTail(double a, double b, double &x, double &y)
{
    double bvirt = a - x;
    double avirt = x + bvirt;
    double bround = bvirt - b;
    double around = a - avirt;
    y = around + bround;
}


inline void TwoDiff(double a, double b, double &x, double &y)
{
    x = a - b;
    TwoDiffTail(a, b, x, y);
}


inline void Split(double a, double &ahi, double &alo) 
{
    double c = splitter * a;
    double abig = c - a;
    ahi = c - abig;
    alo = a - ahi;
}


inline void TwoProductTail(double a, double b, double &x, double &y)
{
    Split(a, ahi, alo);
    Split(b, bhi, blo);
    double err1 = x - (ahi * bhi);
    double err2 = err1 - (alo * bhi);
    double err3 = err2 - (ahi * blo);
    y = (alo * blo) - err3;
}


inline void TwoProduct(double a, double b, double &x, double &y)
{
    x = a * b;
    TwoProductTail(a, b, x, y);
}


/**
 * @brief Similar to TwoProduct where one of the inputs has already been split. Avoids redundant splitting.
 * 
 * @tparam double 
 * @param a 
 * @param b 
 * @param bhi 
 * @param blo 
 * @param x 
 * @param y 
 */
inline void TwoProductPresplit(double a, double b, double bhi, double blo, double &x, double &y)
{
    x = a * b;
    Split(a, ahi, alo);
    double err1 = x - (ahi * bhi);
    double err2 = err1 - (alo * bhi);
    double err3 = err2 - (ahi * blo);
    y = (alo * blo) - err3;
}


/* Two_Product_2Presplit() is Two_Product() where both of the inputs have    */
/*   already been split.  Avoids redundant splitting.                        */
inline void TwoProduct2Presplit(double a, double ahi, double alo, double b, double bhi, double blo, double &x, double &y)
{
    x = a * b;
    err1 = x - (ahi * bhi);
    err2 = err1 - (alo * bhi);
    err3 = err2 - (ahi * blo);
    y = (alo * blo) - err3;
}


/**
 * @brief 
 * @tparam double 
 * @param a 
 * @param x 
 * @param y 
 */
inline void SquareTail(double a, double &x, double &y)
{
    Split(a, ahi, alo);
    err1 = x - (ahi * ahi);
    err3 = err1 - ((ahi + ahi) * alo);
    y = (alo * alo) - err3;
}


/**
 * @brief Faster than TwoProduct. 
 * @tparam double 
 * @param a 
 * @param x 
 * @param y 
 */
inline void Square(double a, double &x, double &y)
{
    x = a * a;
    SquareTail(a, x, y);
}


/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

/**
 * @brief Unrolled version of ExpansionSum for fixed length of two one.
 * 
 * @tparam double 
 */
inline void TwoOneSum(double a1, double a0, double b, double &x2, double &x1, double &x0)
{
    TwoSum(a0, b , _i, x0);
    TwoSum(a1, _i, x2, x1);
}


inline void TwoOneDiff(double a1, double a0, double b, double &x2, double &x1, double &x0) 
{  
    TwoDiff(a0, b , _i, x0);
    TwoSum( a1, _i, x2, x1);
}


inline void Two_Two_Sum(double a1, double a0, double b1, double b0, double &x3, double &x2, double &x1, double &x0) 
{
    TwoOneSum(a1, a0, b0, _j, _0, x0);
    TwoOneSum(_j, _0, b1, x3, x2, x1);
}


inline void Two_Two_Diff(double a1, double a0, double b1, double b0, double &x3, double &x2, double &x1, double &x0) {
    TwoOneDiff(a1, a0, b0, _j, _0, x0);
    TwoOneDiff(_j, _0, b1, x3, x2, x1);
}


inline void Four_One_Sum(double a3, double a2, double a1, double a0, double b, double &x4, double &x3, double &x2, double &x1, double &x0) 
{
    TwoOneSum(a1, a0, b , _j, x1, x0);
    TwoOneSum(a3, a2, _j, x4, x3, x2)
}


inline void Four_Two_Sum(double a3, double a2, double a1, double a0, double b1, double b0, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0);
    Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1);
}


inline void Four_Four_Sum(double a3, double a2, double a1, double a0, double b4, double b3, double b1, double b0, double &x7, double &x6, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0);
    Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2);
}


inline void Eight_One_Sum(double a7, double a6, double a5, double a4, double a3, double a2, double a1, double a0, double b, double &x8, double &x7, double &x6, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Four_One_Sum(a3, a2, a1, a0, b , _j, x3, x2, x1, x0);
    Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4);
}


inline void Eight_Two_Sum(double a7, double a6, double a5, double a4, double a3, double a2, double a1, double a0, double b1, double b0, 
                          double &x9, double &x8, double &x7, double &x6, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, _1, _0, x0);
    Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, x3, x2, x1);
}


inline void Eight_Four_Sum(double a7, double a6, double a5, double a4, double a3, double a2, double a1, double a0, double b4, double b3, double b1, double b0, 
                           double &x11, double &x10, double &x9, double &x8, double &x7, double &x6, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, _2, _1, _0, x1, x0);
    Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2);
}


inline void Two_One_Product(double a1, double a0, double b, double &x3, double &x2, double &x1, double &x0)
{
    Split(b, bhi, blo);
    TwoProductPresplit(a0, b, bhi, blo, _i, x0);
    TwoProductPresplit(a1, b, bhi, blo, _j, _0);
    TwoSum(_i, _0, _k, x1);
    FastTwoSum(_j, _k, x3, x2);
}


inline void Four_One_Product(double a3, double a2, double a1, double a0, double b, double &x7, double &x6, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Split(b, bhi, blo);
    TwoProductPresplit(a0, b, bhi, blo, _i, x0);
    TwoProductPresplit(a1, b, bhi, blo, _j, _0);
    TwoSum(_i, _0, _k, x1);
    FastTwoSum(_j, _k, _i, x2);
    TwoProductPresplit(a2, b, bhi, blo, _j, _0);
    TwoSum(_i, _0, _k, x3);
    FastTwoSum(_j, _k, _i, x4);
    TwoProductPresplit(a3, b, bhi, blo, _j, _0);
    TwoSum(_i, _0, _k, x5);
    FastTwoSum(_j, _k, x7, x6);
}


inline void Two_Two_Product(double a1, double a0, double b1, double b0, double &x7, double &x6, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Split(a0, a0hi, a0lo);
    Split(b0, bhi, blo);
    TwoProduct2Presplit(a0, a0hi, a0lo, b0, bhi, blo, _i, x0);
    Split(a1, a1hi, a1lo);
    TwoProduct2Presplit(a1, a1hi, a1lo, b0, bhi, blo, _j, _0);
    TwoSum(_i, _0, _k, _1);
    FastTwoSum(_j, _k, _l, _2);
    Split(b1, bhi, blo);
    TwoProduct2Presplit(a0, a0hi, a0lo, b1, bhi, blo, _i, _0);
    TwoSum(_1, _0, _k, x1);
    TwoSum(_2, _k, _j, _1);
    TwoSum(_l, _j, _m, _2);
    TwoProduct2Presplit(a1, a1hi, a1lo, b1, bhi, blo, _j, _0);
    TwoSum(_i, _0, _n, _0);
    TwoSum(_1, _0, _i, x2);
    TwoSum(_2, _i, _k, _1);
    TwoSum(_m, _k, _l, _2);
    TwoSum(_j, _n, _k, _0);
    TwoSum(_1, _0, _j, x3);
    TwoSum(_2, _j, _i, _1);
    TwoSum(_l, _i, _m, _2);
    TwoSum(_1, _k, _i, x4);
    TwoSum(_2, _i, _k, x5);
    TwoSum(_m, _k, x7, x6);
}


/* An expansion of length two can be squared more quickly than finding the   */
/*   product of two different expansions of length two, and the result is    */
/*   guaranteed to have no more than six (rather than eight) components.     */
inline void Two_Square(double a1, double a0, double &x5, double &x4, double &x3, double &x2, double &x1, double &x0)
{
    Square(a0, _j, x0);
    _0 = a0 + a0;
    TwoProduct(a1, _0, _k, _1);
    TwoOneSum(_k, _1, _j, _l, _2, x1);
    Square(a1, _j, _1);
    Two_Two_Sum(_j, _1, _l, _2, x5, x4, x3, x2);
}


// The following codes were part of "IEEE 754 floating-point test software"
//          http://www.math.utah.edu/~beebe/software/ieee/
// The original program was "fpinfo2.c".
static double fppow2(int n)
{
    double x, power;
    x = (n < 0) ? (static_cast<double>(1.0)/static_cast<double>(2.0)) : static_cast<double>(2.0);
    n = (n < 0) ? -n : n;
    power = static_cast<double>(1.0);
    while (n-- > 0) power *= x;
    return (power);
}


static double dstore(double x)
{
    return (x);
}

static int test_double()
{
    double x = 1.0;
    int pass = 1;

    while (dstore(1.0 + x/2.0) != 1.0) x /= 2.0;
    if (x != static_cast<double>(fppow2(-52))) pass = 0;

    x = 1.0;
    while (dstore(x / 2.0) != 0.0) x /= 2.0;
    if (x != static_cast<double>(fppow2(-1074)) && x != static_cast<double>(fppow2(-1022))) pass = 0;

    return pass;
}


/*****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-point arithmetic is too much like real arithmetic.              */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/*****************************************************************************/

void exactinit(int noexact, int nofilter, double maxx, double maxy, double maxz)
{
    #ifdef CPU86
    _control87(_PC_53, _MCW_PC); /* Set FPU control word for double precision. */
    #endif // CPU86

    #ifdef LINUX
    // set FPU control word for double precision.
    int cword = 4722;
    _FPU_SETCW(cword);
    #endif // LINUX
    
    test_double();

    int every_other = 1;
    double half = 0.5;
    epsilon = 1.0;
    splitter = 1.0;
    double check = 1.0;
    double lastcheck;
  
    /* Repeatedly divide `epsilon' by two until it is too small to add to    */
    /*   one without causing roundoff.  (Also check if the sum is equal to   */
    /*   the previous sum, for machines that round up instead of using exact */
    /*   rounding.  Not that this library will work on such machines anyway. */
    do {
        lastcheck = check;
        epsilon *= half;
        if (every_other) splitter *= 2.0;
        
        every_other = !every_other;
        check = 1.0 + epsilon;
    } while ((check != 1.0) && (check != lastcheck));
    splitter += 1.0;

    /* Error bounds for orientation and incircle tests. */
    resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
    ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
    ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
    ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
    o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
    o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
    o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
    iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
    iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
    iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
    isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
    isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
    isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;

    // Set TetGen options.  Added by H. Si, 2012-08-23.
    _use_inexact_arith = noexact;
    _use_static_filter = !nofilter;

    // Calculate the two static filters for orient3d() and insphere() tests.
    // Added by H. Si, 2012-08-23.
    // Sort maxx < maxy < maxz. Re-use 'half' for swapping.
    if (maxx > maxz) {
      half = maxx; maxx = maxz; maxz = half;
    }
    if (maxy > maxz) {
      half = maxy; maxy = maxz; maxz = half;
    }
    else if (maxy < maxx) {
      half = maxy; maxy = maxx; maxx = half;
    }

    o3dstaticfilter = 5.1107127829973299e-15 * maxx * maxy * maxz;
    ispstaticfilter = 1.2466136531027298e-13 * maxx * maxy * maxz * (maxz * maxz);

}


/**
 * @brief Add a scalar to an expansion
 * Sets h = e + b.  e and h can be the same.
 *  Maintains the nonoverlapping property.  If round-to-even is used (as 
 *  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent
 *  properties as well.  (That is, if e has one of these properties, so  
 *  will h.) 
 * @param elen 
 * @param e 
 * @param b 
 * @param h 
 * @return int 
 */
int grow_expansion(int elen, double *e, double b, double *h)
{
    double Q = b;
    double Qnew, enow;
    int eindex;
    for (eindex = 0; eindex != elen; ++eindex) {
        enow = e[eindex];
        TwoSum(Q, enow, Qnew, h[eindex]);
        Q = Qnew;
    }
    
    h[eindex] = Q;
    return eindex + 1;
}


/**
 * @brief Add a scalar to an expansion, eliminating zero components from the output expansion.
 *  Sets h = e + b. e and h can be the same.
 *                                                                       
 *  Maintains the nonoverlapping property.  If round-to-even is used (as 
 *  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent
 *  properties as well.  (That is, if e has one of these properties, so  
 *  will h.)                                                  
 * @param elen 
 * @param e 
 * @param b 
 * @param h 
 * @return int 
 */
int grow_expansion_zeroelim(int elen, double *e, double b, double *h)
{
    int hindex = 0;
    double Q = b;

    double hh, Qnew, enow;
    for (int eindex = 0; eindex != elen; ++eindex) {
        enow = e[eindex];
        TwoSum(Q, enow, Qnew, hh);
        Q = Qnew;
        if (hh != 0.0) {
          h[hindex++] = hh;
        }
    }
    if ((Q != 0.0) || (hindex == 0)) {
      h[hindex++] = Q;
    }
    return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  expansion_sum()   Sum two expansions.                                    */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
/*  strongly nonoverlapping property.                                        */
/*                                                                           */
/*****************************************************************************/

/**
 * @brief Sum two expansions. 
 * Sets h = e + f.  e and h can be the same, but f and h cannot.
 * Maintains the nonoverlapping property.  If round-to-even is used (as   
 *  with IEEE 754), maintains the nonadjacent property as well.  (That is, 
 *  if e has one of these properties, so will h.)  Does NOT maintain the   
 *  strongly nonoverlapping property.
 * @param elen 
 * @param e 
 * @param flen 
 * @param f 
 * @param h 
 * @return int 
 */
int expansion_sum(int elen, double *e, int flen, double *f, double *h)
{   
    double Q = f[0];

    double Qnew;
    int hindex;
    double hnow;
    for (hindex = 0; hindex != elen; ++hindex) {
        hnow = e[hindex];
        TwoSum(Q, hnow, Qnew, h[hindex]);
        Q = Qnew;
    }
    h[hindex] = Q;

    int findex, hlast = hindex;
    for (findex = 1; findex != flen; ++findex) {
        Q = f[findex];
        for (hindex = findex; hindex <= hlast; ++hindex) {
            hnow = h[hindex];
            TwoSum(Q, hnow, Qnew, h[hindex]);
            Q = Qnew;
        }
        h[++hlast] = Q;
    }
    return hlast + 1;
}

/*****************************************************************************/
/*                                                                           */
/*  expansion_sum_zeroelim1()   Sum two expansions, eliminating zero         */
/*                              components from the output expansion.        */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
/*  strongly nonoverlapping property.                                        */
/*                                                                           */
/*****************************************************************************/

int expansion_sum_zeroelim1(int elen, double *e, int flen, double *f, double *h)
/* e and h can be the same, but f and h cannot. */
{
  double Q;
   double Qnew;
  int index, findex, hindex, hlast;
  double hnow;
   double bvirt;
  double avirt, bround, around;

  Q = f[0];
  for (hindex = 0; hindex < elen; hindex++) {
    hnow = e[hindex];
    TwoSum(Q, hnow, Qnew, h[hindex]);
    Q = Qnew;
  }
  h[hindex] = Q;
  hlast = hindex;
  for (findex = 1; findex < flen; findex++) {
    Q = f[findex];
    for (hindex = findex; hindex <= hlast; hindex++) {
      hnow = h[hindex];
      TwoSum(Q, hnow, Qnew, h[hindex]);
      Q = Qnew;
    }
    h[++hlast] = Q;
  }
  hindex = -1;
  for (index = 0; index <= hlast; index++) {
    hnow = h[index];
    if (hnow != 0.0) {
      h[++hindex] = hnow;
    }
  }
  if (hindex == -1) {
    return 1;
  } else {
    return hindex + 1;
  }
}

/*****************************************************************************/
/*                                                                           */
/*  expansion_sum_zeroelim2()   Sum two expansions, eliminating zero         */
/*                              components from the output expansion.        */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
/*  strongly nonoverlapping property.                                        */
/*                                                                           */
/*****************************************************************************/

int expansion_sum_zeroelim2(int elen, double *e, int flen, double *f, double *h)
/* e and h can be the same, but f and h cannot. */
{
  double Q, hh;
   double Qnew;
  int eindex, findex, hindex, hlast;
  double enow;
   double bvirt;
  double avirt, bround, around;

  hindex = 0;
  Q = f[0];
  for (eindex = 0; eindex < elen; eindex++) {
    enow = e[eindex];
    TwoSum(Q, enow, Qnew, hh);
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  h[hindex] = Q;
  hlast = hindex;
  for (findex = 1; findex < flen; findex++) {
    hindex = 0;
    Q = f[findex];
    for (eindex = 0; eindex <= hlast; eindex++) {
      enow = h[eindex];
      TwoSum(Q, enow, Qnew, hh);
      Q = Qnew;
      if (hh != 0) {
        h[hindex++] = hh;
      }
    }
    h[hindex] = Q;
    hlast = hindex;
  }
  return hlast + 1;
}

/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum()   Sum two expansions.                               */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/

static int fast_expansion_sum(int elen, double *e, int flen, double *f, double *h)
/* h cannot be e or f. */
{
  double Q;
   double Qnew;
   double bvirt;
  double avirt, bround, around;
  int eindex, findex, hindex;
  double enow, fnow;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;
    enow = e[++eindex];
  } else {
    Q = fnow;
    fnow = f[++findex];
  }
  hindex = 0;
  if ((eindex < elen) && (findex < flen)) {
    if ((fnow > enow) == (fnow > -enow)) {
      FastTwoSum(enow, Q, Qnew, h[0]);
      enow = e[++eindex];
    } else {
      FastTwoSum(fnow, Q, Qnew, h[0]);
      fnow = f[++findex];
    }
    Q = Qnew;
    hindex = 1;
    while ((eindex < elen) && (findex < flen)) {
      if ((fnow > enow) == (fnow > -enow)) {
        TwoSum(Q, enow, Qnew, h[hindex]);
        enow = e[++eindex];
      } else {
        TwoSum(Q, fnow, Qnew, h[hindex]);
        fnow = f[++findex];
      }
      Q = Qnew;
      hindex++;
    }
  }
  while (eindex < elen) {
    TwoSum(Q, enow, Qnew, h[hindex]);
    enow = e[++eindex];
    Q = Qnew;
    hindex++;
  }
  while (findex < flen) {
    TwoSum(Q, fnow, Qnew, h[hindex]);
    fnow = f[++findex];
    Q = Qnew;
    hindex++;
  }
  h[hindex] = Q;
  return hindex + 1;
}

/*****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/*****************************************************************************/
static
int fast_expansion_sum_zeroelim(int elen, double *e, int flen, double *f, double *h)
/* h cannot be e or f. */
{
  double Q;
   double Qnew;
   double hh;
   double bvirt;
  double avirt, bround, around;
  int eindex, findex, hindex;
  double enow, fnow;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    Q = enow;
    enow = e[++eindex];
  } else {
    Q = fnow;
    fnow = f[++findex];
  }
  hindex = 0;
  if ((eindex < elen) && (findex < flen)) {
    if ((fnow > enow) == (fnow > -enow)) {
      FastTwoSum(enow, Q, Qnew, hh);
      enow = e[++eindex];
    } else {
      FastTwoSum(fnow, Q, Qnew, hh);
      fnow = f[++findex];
    }
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
    while ((eindex < elen) && (findex < flen)) {
      if ((fnow > enow) == (fnow > -enow)) {
        TwoSum(Q, enow, Qnew, hh);
        enow = e[++eindex];
      } else {
        TwoSum(Q, fnow, Qnew, hh);
        fnow = f[++findex];
      }
      Q = Qnew;
      if (hh != 0.0) {
        h[hindex++] = hh;
      }
    }
  }
  while (eindex < elen) {
    TwoSum(Q, enow, Qnew, hh);
    enow = e[++eindex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  while (findex < flen) {
    TwoSum(Q, fnow, Qnew, hh);
    fnow = f[++findex];
    Q = Qnew;
    if (hh != 0.0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  linear_expansion_sum()   Sum two expansions.                             */
/*                                                                           */
/*  Sets h = e + f.  See either version of my paper for details.             */
/*                                                                           */
/*  Maintains the nonoverlapping property.  (That is, if e is                */
/*  nonoverlapping, h will be also.)                                         */
/*                                                                           */
/*****************************************************************************/

int linear_expansion_sum(int elen, double *e, int flen, double *f, double *h)
/* h cannot be e or f. */
{
  double Q, q;
   double Qnew;
   double R;
   double bvirt;
  double avirt, bround, around;
  int eindex, findex, hindex;
  double enow, fnow;
  double g0;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    g0 = enow;
    enow = e[++eindex];
  } else {
    g0 = fnow;
    fnow = f[++findex];
  }
  if ((eindex < elen) && ((findex >= flen)
                          || ((fnow > enow) == (fnow > -enow)))) {
    FastTwoSum(enow, g0, Qnew, q);
    enow = e[++eindex];
  } else {
    FastTwoSum(fnow, g0, Qnew, q);
    fnow = f[++findex];
  }
  Q = Qnew;
  for (hindex = 0; hindex < elen + flen - 2; hindex++) {
    if ((eindex < elen) && ((findex >= flen)
                            || ((fnow > enow) == (fnow > -enow)))) {
      FastTwoSum(enow, q, R, h[hindex]);
      enow = e[++eindex];
    } else {
      FastTwoSum(fnow, q, R, h[hindex]);
      fnow = f[++findex];
    }
    TwoSum(Q, R, Qnew, q);
    Q = Qnew;
  }
  h[hindex] = q;
  h[hindex + 1] = Q;
  return hindex + 2;
}

/*****************************************************************************/
/*                                                                           */
/*  linear_expansion_sum_zeroelim()   Sum two expansions, eliminating zero   */
/*                                    components from the output expansion.  */
/*                                                                           */
/*  Sets h = e + f.  See either version of my paper for details.             */
/*                                                                           */
/*  Maintains the nonoverlapping property.  (That is, if e is                */
/*  nonoverlapping, h will be also.)                                         */
/*                                                                           */
/*****************************************************************************/

int linear_expansion_sum_zeroelim(int elen, double *e, int flen, double *f,
                                  double *h)
/* h cannot be e or f. */
{
  double Q, q, hh;
   double Qnew;
   double R;
   double bvirt;
  double avirt, bround, around;
  int eindex, findex, hindex;
  int count;
  double enow, fnow;
  double g0;

  enow = e[0];
  fnow = f[0];
  eindex = findex = 0;
  hindex = 0;
  if ((fnow > enow) == (fnow > -enow)) {
    g0 = enow;
    enow = e[++eindex];
  } else {
    g0 = fnow;
    fnow = f[++findex];
  }
  if ((eindex < elen) && ((findex >= flen)
                          || ((fnow > enow) == (fnow > -enow)))) {
    FastTwoSum(enow, g0, Qnew, q);
    enow = e[++eindex];
  } else {
    FastTwoSum(fnow, g0, Qnew, q);
    fnow = f[++findex];
  }
  Q = Qnew;
  for (count = 2; count < elen + flen; count++) {
    if ((eindex < elen) && ((findex >= flen)
                            || ((fnow > enow) == (fnow > -enow)))) {
      FastTwoSum(enow, q, R, hh);
      enow = e[++eindex];
    } else {
      FastTwoSum(fnow, q, R, hh);
      fnow = f[++findex];
    }
    TwoSum(Q, R, Qnew, q);
    Q = Qnew;
    if (hh != 0) {
      h[hindex++] = hh;
    }
  }
  if (q != 0) {
    h[hindex++] = q;
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  scale_expansion()   Multiply an expansion by a scalar.                   */
/*                                                                           */
/*  Sets h = be.  See either version of my paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/
static
int scale_expansion(int elen, double *e, double b, double *h)
/* e and h cannot be the same. */
{
   double Q;
   double sum;
   double product1;
  double product0;
  int eindex, hindex;
  double enow;
   double bvirt;
  double avirt, bround, around;
   double c;
   double abig;
  double ahi, alo, bhi, blo;
  double err1, err2, err3;

  Split(b, bhi, blo);
  TwoProductPresplit(e[0], b, bhi, blo, Q, h[0]);
  hindex = 1;
  for (eindex = 1; eindex < elen; eindex++) {
    enow = e[eindex];
    TwoProductPresplit(enow, b, bhi, blo, product1, product0);
    TwoSum(Q, product0, sum, h[hindex]);
    hindex++;
    TwoSum(product1, sum, Q, h[hindex]);
    hindex++;
  }
  h[hindex] = Q;
  return elen + elen;
}

/*****************************************************************************/
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.  See either version of my paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/*****************************************************************************/
static
int scale_expansion_zeroelim(int elen, double *e, double b, double *h)
/* e and h cannot be the same. */
{
   double Q, sum;
  double hh;
   double product1;
  double product0;
  int eindex, hindex;
  double enow;
   double bvirt;
  double avirt, bround, around;
   double c;
   double abig;
  double ahi, alo, bhi, blo;
  double err1, err2, err3;

  Split(b, bhi, blo);
  TwoProductPresplit(e[0], b, bhi, blo, Q, hh);
  hindex = 0;
  if (hh != 0) {
    h[hindex++] = hh;
  }
  for (eindex = 1; eindex < elen; eindex++) {
    enow = e[eindex];
    TwoProductPresplit(enow, b, bhi, blo, product1, product0);
    TwoSum(Q, product0, sum, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
    FastTwoSum(product1, sum, Q, hh);
    if (hh != 0) {
      h[hindex++] = hh;
    }
  }
  if ((Q != 0.0) || (hindex == 0)) {
    h[hindex++] = Q;
  }
  return hindex;
}

/*****************************************************************************/
/*                                                                           */
/*  compress()   Compress an expansion.                                      */
/*                                                                           */
/*  See the long version of my paper for details.                            */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), then any nonoverlapping expansion is converted to a      */
/*  nonadjacent expansion.                                                   */
/*                                                                           */
/*****************************************************************************/
static
int compress(int elen, double *e, double *h)
/* e and h may be the same. */
{
  double Q, q;
   double Qnew;
  int eindex, hindex;
   double bvirt;
  double enow, hnow;
  int top, bottom;

  bottom = elen - 1;
  Q = e[bottom];
  for (eindex = elen - 2; eindex >= 0; eindex--) {
    enow = e[eindex];
    FastTwoSum(Q, enow, Qnew, q);
    if (q != 0) {
      h[bottom--] = Qnew;
      Q = q;
    } else {
      Q = Qnew;
    }
  }
  top = 0;
  for (hindex = bottom + 1; hindex < elen; hindex++) {
    hnow = h[hindex];
    FastTwoSum(hnow, Q, Qnew, q);
    if (q != 0) {
      h[top++] = q;
    }
    Q = Qnew;
  }
  h[top] = Q;
  return top + 1;
}

/*****************************************************************************/
/*                                                                           */
/*  estimate()   Produce a one-word estimate of an expansion's value.        */
/*                                                                           */
/*  See either version of my paper for details.                              */
/*                                                                           */
/*****************************************************************************/
static
double estimate(int elen, double *e)
{
  double Q;
  int eindex;

  Q = e[0];
  for (eindex = 1; eindex < elen; eindex++) {
    Q += e[eindex];
  }
  return Q;
}


/*****************************************************************************/
/*                                                                           */
/*  orient2dfast()   Approximate 2D orientation test.  Nonrobust.            */
/*  orient2dslow()   Another exact 2D orientation test.  Robust.             */
/*                                                                           */
/*               Return a positive value if the points pa, pb, and pc occur  */
/*               in counterclockwise order; a negative value if they occur   */
/*               in clockwise order; and zero if they are collinear.  The    */
/*               result is also a rough approximation of twice the signed    */
/*               area of the triangle defined by the three points.           */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In orient2d() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, orient2d() is usually quite */
/*  fast, but will run more slowly when the input points are collinear or    */
/*  nearly so.                                                               */
/*                                                                           */
/*****************************************************************************/
static
double orient2dfast(double *pa, double *pb, double *pc)
{
  double acx, bcx, acy, bcy;

  acx = pa[0] - pc[0];
  bcx = pb[0] - pc[0];
  acy = pa[1] - pc[1];
  bcy = pb[1] - pc[1];
  return acx * bcy - acy * bcx;
}


// orient2dexact()   Exact 2D orientation test.  Robust.   
double orient2dexact(double *pa, double *pb, double *pc)
{
    double axby1, axcy1, bxcy1, bxay1, cxay1, cxby1;
    double axby0, axcy0, bxcy0, bxay0, cxay0, cxby0;
    double aterms[4], bterms[4], cterms[4];
    double aterms3, bterms3, cterms3;
    double v[8], w[12];
    int vlength, wlength;

    double bvirt;
    double avirt, bround, around;
    double c;
    double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
    double _i, _j;
    double _0;

    TwoProduct(pa[0], pb[1], axby1, axby0);
    TwoProduct(pa[0], pc[1], axcy1, axcy0);
    Two_Two_Diff(axby1, axby0, axcy1, axcy0,
                 aterms3, aterms[2], aterms[1], aterms[0]);
    aterms[3] = aterms3;

    TwoProduct(pb[0], pc[1], bxcy1, bxcy0);
    TwoProduct(pb[0], pa[1], bxay1, bxay0);
    Two_Two_Diff(bxcy1, bxcy0, bxay1, bxay0,
                 bterms3, bterms[2], bterms[1], bterms[0]);
    bterms[3] = bterms3;

    TwoProduct(pc[0], pa[1], cxay1, cxay0);
    TwoProduct(pc[0], pb[1], cxby1, cxby0);
    Two_Two_Diff(cxay1, cxay0, cxby1, cxby0,
                 cterms3, cterms[2], cterms[1], cterms[0]);
    cterms[3] = cterms3;

    vlength = fast_expansion_sum_zeroelim(4, aterms, 4, bterms, v);
    wlength = fast_expansion_sum_zeroelim(vlength, v, 4, cterms, w);

    return w[wlength - 1];
}


double orient2dslow(double *pa, double *pb, double *pc)
{
    double acx, acy, bcx, bcy;
    double acxtail, acytail;
    double bcxtail, bcytail;
    double negate, negatetail;
    double axby[8], bxay[8];
    double axby7, bxay7;
    double deter[16];
    int deterlen;

    double bvirt;
    double avirt, bround, around;
    double c;
    double abig;
    double a0hi, a0lo, a1hi, a1lo, bhi, blo;
    double err1, err2, err3;
    double _i, _j, _k, _l, _m, _n;
    double _0, _1, _2;

    TwoDiff(pa[0], pc[0], acx, acxtail);
    TwoDiff(pa[1], pc[1], acy, acytail);
    TwoDiff(pb[0], pc[0], bcx, bcxtail);
    TwoDiff(pb[1], pc[1], bcy, bcytail);

    Two_Two_Product(acx, acxtail, bcy, bcytail,
                    axby7, axby[6], axby[5], axby[4],
                    axby[3], axby[2], axby[1], axby[0]);
    axby[7] = axby7;
    negate = -acy;
    negatetail = -acytail;
    Two_Two_Product(bcx, bcxtail, negate, negatetail,
                    bxay7, bxay[6], bxay[5], bxay[4],
                    bxay[3], bxay[2], bxay[1], bxay[0]);
    bxay[7] = bxay7;

    deterlen = fast_expansion_sum_zeroelim(8, axby, 8, bxay, deter);

    return deter[deterlen - 1];
}


double orient2dadapt(double *pa, double *pb, double *pc, double detsum)
{
    double acx, acy, bcx, bcy;
    double acxtail, acytail, bcxtail, bcytail;
    double detleft, detright;
    double detlefttail, detrighttail;
    double det, errbound;
    double B[4], C1[8], C2[12], D[16];
    double B3;
    int C1length, C2length, Dlength;
    double u[4];
    double u3;
    double s1, t1;
    double s0, t0;

    double bvirt;
    double avirt, bround, around;
    double c;
    double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
    double _i, _j;
    double _0;

    acx = (double) (pa[0] - pc[0]);
    bcx = (double) (pb[0] - pc[0]);
    acy = (double) (pa[1] - pc[1]);
    bcy = (double) (pb[1] - pc[1]);

    TwoProduct(acx, bcy, detleft, detlefttail);
    TwoProduct(acy, bcx, detright, detrighttail);

    Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
                 B3, B[2], B[1], B[0]);
    B[3] = B3;

    det = estimate(4, B);
    errbound = ccwerrboundB * detsum;
    if ((det >= errbound) || (-det >= errbound)) { return det; }

    TwoDiffTail(pa[0], pc[0], acx, acxtail);
    TwoDiffTail(pb[0], pc[0], bcx, bcxtail);
    TwoDiffTail(pa[1], pc[1], acy, acytail);
    TwoDiffTail(pb[1], pc[1], bcy, bcytail);

    if ((acxtail == 0.0) && (acytail == 0.0) && 
        (bcxtail == 0.0) && (bcytail == 0.0)) {
        return det;
    }

    errbound = ccwerrboundC * detsum + resulterrbound * std::fabs(det);
    det += (acx * bcytail + bcy * acxtail) - (acy * bcxtail + bcx * acytail);
    if ((det >= errbound) || (-det >= errbound)) { return det; }

    TwoProduct(acxtail, bcy, s1, s0);
    TwoProduct(acytail, bcx, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    C1length = fast_expansion_sum_zeroelim(4, B, 4, u, C1);

    TwoProduct(acx, bcytail, s1, s0);
    TwoProduct(acy, bcxtail, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    C2length = fast_expansion_sum_zeroelim(C1length, C1, 4, u, C2);

    TwoProduct(acxtail, bcytail, s1, s0);
    TwoProduct(acytail, bcxtail, t1, t0);
    Two_Two_Diff(s1, s0, t1, t0, u3, u[2], u[1], u[0]);
    u[3] = u3;
    Dlength = fast_expansion_sum_zeroelim(C2length, C2, 4, u, D);

    return(D[Dlength - 1]);
}


/*  orient2d()   Adaptive exact 2D orientation test.  Robust.                */
double orient2d(double *pa, double *pb, double *pc)
{
    double detleft, detright, det;
    double detsum, errbound;

    detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    det = detleft - detright;

    if (detleft > 0.0) {
        if (detright <= 0.0) { return det; } 
        else { detsum = detleft + detright; }
    } else if (detleft < 0.0) {
        if (detright >= 0.0) { return det; } 
        else { detsum = -detleft - detright; }
    } else {
        return det;
    }

    errbound = ccwerrboundA * detsum;
    if ((det >= errbound) || (-det >= errbound)) { return det; }

    return orient2dadapt(pa, pb, pc, detsum);
}


/*****************************************************************************/
/*                                                                           */

/*                                                                           */
/*               Return a positive value if the point pd lies below the      */
/*               plane passing through pa, pb, and pc; "below" is defined so */
/*               that pa, pb, and pc appear in counterclockwise order when   */
/*               viewed from above the plane.  Returns a negative value if   */
/*               pd lies above the plane.  Returns zero if the points are    */
/*               coplanar.  The result is also a rough approximation of six  */
/*               times the signed volume of the tetrahedron defined by the   */
/*               four points.                                                */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In orient3d() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, orient3d() is usually quite */
/*  fast, but will run more slowly when the input points are coplanar or     */
/*  nearly so.                                                               */
/*                                                                           */
/*****************************************************************************/

/*  orient3dfast()   Approximate 3D orientation test.  Nonrobust.  */
double orient3dfast(double *pa, double *pb, double *pc, double *pd)
{
    double adx = pa[0] - pd[0]; 
    double ady = pa[1] - pd[1]; 
    double adz = pa[2] - pd[2];

    double bdx = pb[0] - pd[0];    
    double bdy = pb[1] - pd[1];
    double bdz = pb[2] - pd[2];

    double cdx = pc[0] - pd[0];
    double cdy = pc[1] - pd[1]; 
    double cdz = pc[2] - pd[2];

    return adx * (bdy * cdz - bdz * cdy) +
           bdx * (cdy * adz - cdz * ady) + 
           cdx * (ady * bdz - adz * bdy);
}


/*  orient3dexact()   Exact 3D orientation test.  Robust.                    */
double orient3dexact(double *pa, double *pb, double *pc, double *pd)
{
    double axby1, bxcy1, cxdy1, dxay1, axcy1, bxdy1;
    double bxay1, cxby1, dxcy1, axdy1, cxay1, dxby1;
    double axby0, bxcy0, cxdy0, dxay0, axcy0, bxdy0;
    double bxay0, cxby0, dxcy0, axdy0, cxay0, dxby0;
    double ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
    double temp8[8];
    int templen;
    double abc[12], bcd[12], cda[12], dab[12];
    int abclen, bcdlen, cdalen, dablen;
    double adet[24], bdet[24], cdet[24], ddet[24];
    int alen, blen, clen, dlen;
    double abdet[48], cddet[48];
    int ablen, cdlen;
    double deter[96];
    int deterlen;
    int i;

    double bvirt;
    double avirt, bround, around;
    double c;
    double abig;
    double ahi, alo, bhi, blo;
    double err1, err2, err3;
    double _i, _j;
    double _0;

    TwoProduct(pa[0], pb[1], axby1, axby0);
    TwoProduct(pb[0], pa[1], bxay1, bxay0);
    Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

    TwoProduct(pb[0], pc[1], bxcy1, bxcy0);
    TwoProduct(pc[0], pb[1], cxby1, cxby0);
    Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

    TwoProduct(pc[0], pd[1], cxdy1, cxdy0);
    TwoProduct(pd[0], pc[1], dxcy1, dxcy0);
    Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

    TwoProduct(pd[0], pa[1], dxay1, dxay0);
    TwoProduct(pa[0], pd[1], axdy1, axdy0);
    Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

    TwoProduct(pa[0], pc[1], axcy1, axcy0);
    TwoProduct(pc[0], pa[1], cxay1, cxay0);
    Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

    TwoProduct(pb[0], pd[1], bxdy1, bxdy0);
    TwoProduct(pd[0], pb[1], dxby1, dxby0);
    Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

    templen = fast_expansion_sum_zeroelim(4, cd, 4, da, temp8);
    cdalen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, cda);
    templen = fast_expansion_sum_zeroelim(4, da, 4, ab, temp8);
    dablen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, dab);
    for (i = 0; i < 4; i++) {
        bd[i] = -bd[i];
        ac[i] = -ac[i];
    }
    templen = fast_expansion_sum_zeroelim(4, ab, 4, bc, temp8);
    abclen = fast_expansion_sum_zeroelim(templen, temp8, 4, ac, abc);
    templen = fast_expansion_sum_zeroelim(4, bc, 4, cd, temp8);
    bcdlen = fast_expansion_sum_zeroelim(templen, temp8, 4, bd, bcd);

    alen = scale_expansion_zeroelim(bcdlen, bcd, pa[2], adet);
    blen = scale_expansion_zeroelim(cdalen, cda, -pb[2], bdet);
    clen = scale_expansion_zeroelim(dablen, dab, pc[2], cdet);
    dlen = scale_expansion_zeroelim(abclen, abc, -pd[2], ddet);

    ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
    deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, deter);

    return deter[deterlen - 1];
}


/*  orient3dslow()   Another exact 3D orientation test.  Robust.             */
double orient3dslow(double *pa, double *pb, double *pc, double *pd)
{
    double adx, ady, adz, bdx, bdy, bdz, cdx, cdy, cdz;
    double adxtail, adytail, adztail;
    double bdxtail, bdytail, bdztail;
    double cdxtail, cdytail, cdztail;
    double negate, negatetail;
    double axby7, bxcy7, axcy7, bxay7, cxby7, cxay7;
    double axby[8], bxcy[8], axcy[8], bxay[8], cxby[8], cxay[8];
    double temp16[16], temp32[32], temp32t[32];
    int temp16len, temp32len, temp32tlen;
    double adet[64], bdet[64], cdet[64];
    int alen, blen, clen;
    double abdet[128];
    int ablen;
    double deter[192];
    int deterlen;

    double bvirt;
    double avirt, bround, around;
    double c;
    double abig;
    double a0hi, a0lo, a1hi, a1lo, bhi, blo;
    double err1, err2, err3;
    double _i, _j, _k, _l, _m, _n;
    double _0, _1, _2;

    TwoDiff(pa[0], pd[0], adx, adxtail);
    TwoDiff(pa[1], pd[1], ady, adytail);
    TwoDiff(pa[2], pd[2], adz, adztail);
    TwoDiff(pb[0], pd[0], bdx, bdxtail);
    TwoDiff(pb[1], pd[1], bdy, bdytail);
    TwoDiff(pb[2], pd[2], bdz, bdztail);
    TwoDiff(pc[0], pd[0], cdx, cdxtail);
    TwoDiff(pc[1], pd[1], cdy, cdytail);
    TwoDiff(pc[2], pd[2], cdz, cdztail);

    Two_Two_Product(adx, adxtail, bdy, bdytail,
                    axby7, axby[6], axby[5], axby[4],
                    axby[3], axby[2], axby[1], axby[0]);
    axby[7] = axby7;
    negate = -ady;
    negatetail = -adytail;
    Two_Two_Product(bdx, bdxtail, negate, negatetail,
                    bxay7, bxay[6], bxay[5], bxay[4],
                    bxay[3], bxay[2], bxay[1], bxay[0]);
    bxay[7] = bxay7;
    Two_Two_Product(bdx, bdxtail, cdy, cdytail,
                    bxcy7, bxcy[6], bxcy[5], bxcy[4],
                    bxcy[3], bxcy[2], bxcy[1], bxcy[0]);
    bxcy[7] = bxcy7;
    negate = -bdy;
    negatetail = -bdytail;
    Two_Two_Product(cdx, cdxtail, negate, negatetail,
                    cxby7, cxby[6], cxby[5], cxby[4],
                    cxby[3], cxby[2], cxby[1], cxby[0]);
    cxby[7] = cxby7;
    Two_Two_Product(cdx, cdxtail, ady, adytail,
                    cxay7, cxay[6], cxay[5], cxay[4],
                    cxay[3], cxay[2], cxay[1], cxay[0]);
    cxay[7] = cxay7;
    negate = -cdy;
    negatetail = -cdytail;
    Two_Two_Product(adx, adxtail, negate, negatetail,
                    axcy7, axcy[6], axcy[5], axcy[4],
                    axcy[3], axcy[2], axcy[1], axcy[0]);
    axcy[7] = axcy7;

    temp16len = fast_expansion_sum_zeroelim(8, bxcy, 8, cxby, temp16);
    temp32len = scale_expansion_zeroelim(temp16len, temp16, adz, temp32);
    temp32tlen = scale_expansion_zeroelim(temp16len, temp16, adztail, temp32t);
    alen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, adet);

    temp16len = fast_expansion_sum_zeroelim(8, cxay, 8, axcy, temp16);
    temp32len = scale_expansion_zeroelim(temp16len, temp16, bdz, temp32);
    temp32tlen = scale_expansion_zeroelim(temp16len, temp16, bdztail, temp32t);
    blen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, bdet);

    temp16len = fast_expansion_sum_zeroelim(8, axby, 8, bxay, temp16);
    temp32len = scale_expansion_zeroelim(temp16len, temp16, cdz, temp32);
    temp32tlen = scale_expansion_zeroelim(temp16len, temp16, cdztail, temp32t);
    clen = fast_expansion_sum_zeroelim(temp32len, temp32, temp32tlen, temp32t, cdet);

    ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
    deterlen = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, deter);

    return deter[deterlen - 1];
}

double orient3dadapt(double *pa, double *pb, double *pc, double *pd, double permanent)
{
    double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
    double det, errbound;

   double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
  double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
  double bc[4], ca[4], ab[4];
   double bc3, ca3, ab3;
  double adet[8], bdet[8], cdet[8];
  int alen, blen, clen;
  double abdet[16];
  int ablen;
  double *finnow, *finother, *finswap;
  double fin1[192], fin2[192];
  int finlength;


  double adxtail, bdxtail, cdxtail;
  double adytail, bdytail, cdytail;
  double adztail, bdztail, cdztail;
   double at_blarge, at_clarge;
   double bt_clarge, bt_alarge;
   double ct_alarge, ct_blarge;
  double at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
  int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
   double bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
   double adxt_cdy1, adxt_bdy1, bdxt_ady1;
  double bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
  double adxt_cdy0, adxt_bdy0, bdxt_ady0;
   double bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
   double adyt_cdx1, adyt_bdx1, bdyt_adx1;
  double bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
  double adyt_cdx0, adyt_bdx0, bdyt_adx0;
  double bct[8], cat[8], abt[8];
  int bctlen, catlen, abtlen;
   double bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
   double adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
  double bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
  double adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
  double u[4], v[12], w[16];
   double u3;
  int vlength, wlength;
  double negate;

   double bvirt;
  double avirt, bround, around;
   double c;
   double abig;
  double ahi, alo, bhi, blo;
  double err1, err2, err3;
   double _i, _j, _k;
  double _0;


  adx = (double) (pa[0] - pd[0]);
  bdx = (double) (pb[0] - pd[0]);
  cdx = (double) (pc[0] - pd[0]);
  ady = (double) (pa[1] - pd[1]);
  bdy = (double) (pb[1] - pd[1]);
  cdy = (double) (pc[1] - pd[1]);
  adz = (double) (pa[2] - pd[2]);
  bdz = (double) (pb[2] - pd[2]);
  cdz = (double) (pc[2] - pd[2]);

  TwoProduct(bdx, cdy, bdxcdy1, bdxcdy0);
  TwoProduct(cdx, bdy, cdxbdy1, cdxbdy0);
  Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
  bc[3] = bc3;
  alen = scale_expansion_zeroelim(4, bc, adz, adet);

  TwoProduct(cdx, ady, cdxady1, cdxady0);
  TwoProduct(adx, cdy, adxcdy1, adxcdy0);
  Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
  ca[3] = ca3;
  blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

  TwoProduct(adx, bdy, adxbdy1, adxbdy0);
  TwoProduct(bdx, ady, bdxady1, bdxady0);
  Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
  ab[3] = ab3;
  clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

  ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
  finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

  det = estimate(finlength, fin1);
  errbound = o3derrboundB * permanent;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  TwoDiffTail(pa[0], pd[0], adx, adxtail);
  TwoDiffTail(pb[0], pd[0], bdx, bdxtail);
  TwoDiffTail(pc[0], pd[0], cdx, cdxtail);
  TwoDiffTail(pa[1], pd[1], ady, adytail);
  TwoDiffTail(pb[1], pd[1], bdy, bdytail);
  TwoDiffTail(pc[1], pd[1], cdy, cdytail);
  TwoDiffTail(pa[2], pd[2], adz, adztail);
  TwoDiffTail(pb[2], pd[2], bdz, bdztail);
  TwoDiffTail(pc[2], pd[2], cdz, cdztail);

  if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
      && (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
      && (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
    return det;
  }

  errbound = o3derrboundC * permanent + resulterrbound * std::fabs(det);
  det += (adz * ((bdx * cdytail + cdy * bdxtail)
                 - (bdy * cdxtail + cdx * bdytail))
          + adztail * (bdx * cdy - bdy * cdx))
       + (bdz * ((cdx * adytail + ady * cdxtail)
                 - (cdy * adxtail + adx * cdytail))
          + bdztail * (cdx * ady - cdy * adx))
       + (cdz * ((adx * bdytail + bdy * adxtail)
                 - (ady * bdxtail + bdx * adytail))
          + cdztail * (adx * bdy - ady * bdx));
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }

  finnow = fin1;
  finother = fin2;

  if (adxtail == 0.0) {
    if (adytail == 0.0) {
      at_b[0] = 0.0;
      at_blen = 1;
      at_c[0] = 0.0;
      at_clen = 1;
    } else {
      negate = -adytail;
      TwoProduct(negate, bdx, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      TwoProduct(adytail, cdx, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    }
  } else {
    if (adytail == 0.0) {
      TwoProduct(adxtail, bdy, at_blarge, at_b[0]);
      at_b[1] = at_blarge;
      at_blen = 2;
      negate = -adxtail;
      TwoProduct(negate, cdy, at_clarge, at_c[0]);
      at_c[1] = at_clarge;
      at_clen = 2;
    } else {
      TwoProduct(adxtail, bdy, adxt_bdy1, adxt_bdy0);
      TwoProduct(adytail, bdx, adyt_bdx1, adyt_bdx0);
      Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
                   at_blarge, at_b[2], at_b[1], at_b[0]);
      at_b[3] = at_blarge;
      at_blen = 4;
      TwoProduct(adytail, cdx, adyt_cdx1, adyt_cdx0);
      TwoProduct(adxtail, cdy, adxt_cdy1, adxt_cdy0);
      Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
                   at_clarge, at_c[2], at_c[1], at_c[0]);
      at_c[3] = at_clarge;
      at_clen = 4;
    }
  }
  if (bdxtail == 0.0) {
    if (bdytail == 0.0) {
      bt_c[0] = 0.0;
      bt_clen = 1;
      bt_a[0] = 0.0;
      bt_alen = 1;
    } else {
      negate = -bdytail;
      TwoProduct(negate, cdx, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      TwoProduct(bdytail, adx, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    }
  } else {
    if (bdytail == 0.0) {
      TwoProduct(bdxtail, cdy, bt_clarge, bt_c[0]);
      bt_c[1] = bt_clarge;
      bt_clen = 2;
      negate = -bdxtail;
      TwoProduct(negate, ady, bt_alarge, bt_a[0]);
      bt_a[1] = bt_alarge;
      bt_alen = 2;
    } else {
      TwoProduct(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
      TwoProduct(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
      Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
                   bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
      bt_c[3] = bt_clarge;
      bt_clen = 4;
      TwoProduct(bdytail, adx, bdyt_adx1, bdyt_adx0);
      TwoProduct(bdxtail, ady, bdxt_ady1, bdxt_ady0);
      Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
                  bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
      bt_a[3] = bt_alarge;
      bt_alen = 4;
    }
  }
  if (cdxtail == 0.0) {
    if (cdytail == 0.0) {
      ct_a[0] = 0.0;
      ct_alen = 1;
      ct_b[0] = 0.0;
      ct_blen = 1;
    } else {
      negate = -cdytail;
      TwoProduct(negate, adx, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      TwoProduct(cdytail, bdx, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    }
  } else {
    if (cdytail == 0.0) {
      TwoProduct(cdxtail, ady, ct_alarge, ct_a[0]);
      ct_a[1] = ct_alarge;
      ct_alen = 2;
      negate = -cdxtail;
      TwoProduct(negate, bdy, ct_blarge, ct_b[0]);
      ct_b[1] = ct_blarge;
      ct_blen = 2;
    } else {
      TwoProduct(cdxtail, ady, cdxt_ady1, cdxt_ady0);
      TwoProduct(cdytail, adx, cdyt_adx1, cdyt_adx0);
      Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
                   ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
      ct_a[3] = ct_alarge;
      ct_alen = 4;
      TwoProduct(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
      TwoProduct(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
      Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
                   ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
      ct_b[3] = ct_blarge;
      ct_blen = 4;
    }
  }

  bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
  wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
  wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
  wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
  finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                          finother);
  finswap = finnow; finnow = finother; finother = finswap;

  if (adztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, bc, adztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  if (adxtail != 0.0) {
    if (bdytail != 0.0) {
      TwoProduct(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
      Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (cdytail != 0.0) {
      negate = -adxtail;
      TwoProduct(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
      Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (bdxtail != 0.0) {
    if (cdytail != 0.0) {
      TwoProduct(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
      Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (adytail != 0.0) {
      negate = -bdxtail;
      TwoProduct(negate, adytail, bdxt_adyt1, bdxt_adyt0);
      Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (cdztail != 0.0) {
        Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }
  if (cdxtail != 0.0) {
    if (adytail != 0.0) {
      TwoProduct(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
      Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (bdztail != 0.0) {
        Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
    if (bdytail != 0.0) {
      negate = -cdxtail;
      TwoProduct(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
      Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
      u[3] = u3;
      finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                              finother);
      finswap = finnow; finnow = finother; finother = finswap;
      if (adztail != 0.0) {
        Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
        u[3] = u3;
        finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
                                                finother);
        finswap = finnow; finnow = finother; finother = finswap;
      }
    }
  }

  if (adztail != 0.0) {
    wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (bdztail != 0.0) {
    wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }
  if (cdztail != 0.0) {
    wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
    finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
                                            finother);
    finswap = finnow; finnow = finother; finother = finswap;
  }

  return finnow[finlength - 1];
}

/*  orient3d()   Adaptive exact 3D orientation test.  Robust.                */
double orient3d(double *pa, double *pb, double *pc, double *pd)
{
    double adx = pa[0] - pd[0];
    double ady = pa[1] - pd[1];
    double adz = pa[2] - pd[2];

    double bdx = pb[0] - pd[0]; 
    double bdy = pb[1] - pd[1]; 
    double bdz = pb[2] - pd[2];

    double cdx = pc[0] - pd[0];
    double cdy = pc[1] - pd[1];
    double cdz = pc[2] - pd[2];
  
    double bdxcdy = bdx * cdy;
    double cdxbdy = cdx * bdy;

    double cdxady = cdx * ady;
    double adxcdy = adx * cdy;

    double adxbdy = adx * bdy;
    double bdxady = bdx * ady;

    double det = adz * (bdxcdy - cdxbdy) +
                 bdz * (cdxady - adxcdy) + 
                 cdz * (adxbdy - bdxady);

    if (_use_inexact_arith)  return det;

    if (_use_static_filter) {
      if (det > o3dstaticfilter) return det;
      if (det < -o3dstaticfilter) return det;
    }

    double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * std::fabs(adz) + 
                       (std::fabs(cdxady) + std::fabs(adxcdy)) * std::fabs(bdz) + 
                       (std::fabs(adxbdy) + std::fabs(bdxady)) * std::fabs(cdz);
    double errbound = o3derrboundA * permanent;
    
    if ((det > errbound) || (-det > errbound)) return det;
    return orient3dadapt(pa, pb, pc, pd, permanent);
}



/** \} End of Doxygen Groups*/
} //end of namespace PREDICATES
} //end of namespace IMP


#endif //IMP_ENGINE_UTILITIES_PREDICATES_HPP_