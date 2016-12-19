// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmootherLocalNewton_hpp
#define ReferenceMeshSmootherLocalNewton_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>

  namespace percept {

    class ReferenceMeshSmootherLocalNewton : public ReferenceMeshSmootherConjugateGradientImpl<STKMesh> {

    public:

      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length
      ReferenceMeshSmootherLocalNewton(PerceptMesh *eMesh,
                            stk::mesh::Selector *boundary_selector=0,
                            MeshGeometry *meshGeometry=0,
                            int inner_iterations = 100,
                            double grad_norm =1.e-8,
                                       int parallel_iterations = 20);

    protected:

      virtual void get_step();
      virtual double run_one_iteration();
      double local_line_search(stk::mesh::Entity node, stk::mesh::FieldBase *cg_g_field);

    };

  }

#endif
#endif
