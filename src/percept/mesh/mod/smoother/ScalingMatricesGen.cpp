// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/ScalingMatricesGen.hpp>

namespace percept {
  ScalingMatrices ScalingMatrices::s_scalingMatrices;
  Indices Indices::s_indices;
}

#endif
