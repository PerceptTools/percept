// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef LocalQuadraticRationalSplineFit_hpp
#define LocalQuadraticRationalSplineFit_hpp

#include <percept/mesh/geometry/stk_geom/BSplineFit.hpp>

  namespace geom {

    class LocalQuadraticRationalSplineFit : public BSplineFit
    {
    public:
      // computes m_CV, m_U
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T);
    };

  }

#endif
