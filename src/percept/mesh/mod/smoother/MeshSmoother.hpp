// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef MeshSmoother_hpp
#define MeshSmoother_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/PerceptMesh.hpp>

#if defined(STK_PERCEPT_HAS_GEOMETRY)
#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#endif

#undef USE_CALLGRIND_MESH_SMOOTHER
//#define USE_CALLGRIND_MESH_SMOOTHER
#ifdef USE_CALLGRIND_MESH_SMOOTHER
#include "/usr/netpub/valgrind-3.8.1/include/valgrind/callgrind.h"
#endif

  namespace percept {

#if !defined(STK_PERCEPT_HAS_GEOMETRY)
    class MeshGeometry {};
#endif

    const double PS_DEF_UNT_BETA = 1e-8;
    const double PS_DEF_SUC_EPS = 1e-4;


    enum NodeClassifyType {
      MS_VERTEX,
      MS_CURVE,
      MS_SURFACE,
      MS_VOLUME,
      MS_ON_BOUNDARY,
      MS_NOT_ON_BOUNDARY
    };

    /// Abstract base class smoother
    class MeshSmoother
    {

    protected:
      PerceptMesh *m_eMesh;
      int innerIter;
      double gradNorm;
      int parallelIterations;
      stk::mesh::Selector *m_boundarySelector;
    public:
      MeshGeometry *m_meshGeometry;

    public:

      MeshSmoother(PerceptMesh *eMesh,
                   stk::mesh::Selector *boundary_selector=0,
                   MeshGeometry *meshGeometry=0,
                   int innerIter=100, double gradNorm = 1.e-8, int parallelIterations=20) :
        m_eMesh(eMesh), innerIter(innerIter), gradNorm(gradNorm), parallelIterations(parallelIterations),
        m_boundarySelector(boundary_selector), m_meshGeometry(meshGeometry)
      {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
        if (m_meshGeometry)
          m_meshGeometry->m_cache_classify_bucket_is_active = true;
#endif
      }

      /// @deprecated
      static size_t count_invalid_elements();

      static size_t parallel_count_invalid_elements(PerceptMesh *eMesh);

      void run( bool always_smooth=true, int debug=0);
      virtual void run_algorithm() = 0;

      static bool select_bucket(stk::mesh::Bucket& bucket, PerceptMesh *eMesh);
      std::pair<bool,int> get_fixed_flag(stk::mesh::Entity node_ptr);
      int classify_node(stk::mesh::Entity node, size_t& curveOrSurfaceEvaluator) const;
      void project_delta_to_tangent_plane(stk::mesh::Entity node, double *delta, double *norm=0);
      void enforce_tangent_plane(stk::mesh::Entity node, double rhs[3], double lhs[3][3], double *norm=0);
      /// if reset is true, don't actually modify the node's coordinates and only
      /// return the snapped position in @param coordinate
      void snap_to(stk::mesh::Entity node_ptr,  double *coordinate, bool reset=false) const;


    };

  }

#endif
#endif
