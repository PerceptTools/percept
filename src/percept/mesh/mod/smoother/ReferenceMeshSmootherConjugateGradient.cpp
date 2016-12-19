// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <percept/math/DenseMatrix.hpp>
#include <percept/math/Math.hpp>

#include <percept/mesh/geometry/kernel/GeometryKernel.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <percept/MeshType.hpp>
//#include <percept/MeshTypeSpec.hpp>

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradientDef.hpp>

#include <stdio.h>
#include <limits>

#include "mpi.h"
#include <cstdio>


namespace percept {



  template class ReferenceMeshSmootherConjugateGradientImpl<STKMesh>;
  template class ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid>;

}



#endif

