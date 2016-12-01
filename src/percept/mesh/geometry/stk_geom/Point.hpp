// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Point_hpp
#define Point_hpp

#include <opennurbs.h>
#include <vector>

  namespace geom {

    typedef ON_2dVector Vector2D;
    typedef std::vector<Vector2D> Vectors2D;

    typedef ON_3dVector Vector3D;
    typedef ON_3dPoint Point3D;
    typedef std::vector<Vector3D> Vectors3D;
    typedef std::vector<Point3D> Points3D;

  }
#endif
