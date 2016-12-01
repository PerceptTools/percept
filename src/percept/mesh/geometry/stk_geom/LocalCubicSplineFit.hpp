// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef LocalCubicSplineFit_hpp
#define LocalCubicSplineFit_hpp

#include <percept/mesh/geometry/stk_geom/BSplineFit.hpp>

  namespace geom {

    /// see P&T The Nurbs Book, 2nd Edition, section 9.3.4
    class LocalCubicSplineFit : public BSplineFit
    {
    public:
      LocalCubicSplineFit(Option option = ThreePoint) : BSplineFit(option) {}

      /// create an OpenNURBS curve that fits the given input points
      virtual void fit_internal(int n, Vectors2D& Q, Vectors2D& T);
    };

  }

#endif
