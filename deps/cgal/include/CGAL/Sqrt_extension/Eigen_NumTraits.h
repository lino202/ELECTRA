// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-5.0.1/Number_types/include/CGAL/Sqrt_extension/Eigen_NumTraits.h $
// $Id: Eigen_NumTraits.h 52164b1 2019-10-19T15:34:59+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>


#ifndef CGAL_SQRT_EXTENSION_EIGEN_NUMTRAITS_H
#define CGAL_SQRT_EXTENSION_EIGEN_NUMTRAITS_H

namespace Eigen {
  template<class> struct NumTraits;
  template <class NT,class ROOT, class ACDE_TAG, class FP_TAG>
  struct NumTraits<CGAL::Sqrt_extension<NT, ROOT, ACDE_TAG, FP_TAG> >
  {
    typedef CGAL::Sqrt_extension<NT, ROOT, ACDE_TAG, FP_TAG> Real;
    typedef Real NonInteger;
    typedef Real Nested;
    typedef Real Literal;

    static inline Real epsilon() { return NumTraits<NT>::epsilon(); }

    enum {
      IsInteger = 0,
      IsSigned = 1,
      IsComplex = 0,
      RequireInitialization = 1,
      ReadCost = 2*NumTraits<NT>::ReadCost+NumTraits<ROOT>::ReadCost,
      AddCost = 2*NumTraits<NT>::AddCost+NumTraits<ROOT>::ReadCost,
      MulCost = 5*NumTraits<NT>::MulCost+2*NumTraits<NT>::AddCost
    };
  };
}

#endif
