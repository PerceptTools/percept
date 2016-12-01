// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmoother1_hpp
#define ReferenceMeshSmoother1_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/ReferenceMeshSmoother.hpp>

  namespace percept {

    /// A Jacobian based optimization smoother, e.g. 1/A - 1/W, W/A - I, etc. (A = local current Jacobian, W is for original mesh)
    /// Conjugate-gradient version, element-based metrics
    class ReferenceMeshSmoother1 : public ReferenceMeshSmoother {

    public:

      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length
      ReferenceMeshSmoother1(PerceptMesh *eMesh,
                            stk::mesh::Selector *boundary_selector=0,
                            MeshGeometry *meshGeometry=0,
                            int inner_iterations = 100,
                            double grad_norm =1.e-8,
                            int parallel_iterations = 20)
        : ReferenceMeshSmoother(eMesh, boundary_selector, meshGeometry, inner_iterations, grad_norm, parallel_iterations)
        ,m_max_edge_length_factor(1.0)
      {}

    protected:
      double m_max_edge_length_factor;

      virtual void get_gradient();
      //virtual void get_scale();
      void debug_print(double alpha);

      typedef long double Double;
      virtual double run_one_iteration();

      virtual Double total_metric( Double alpha, double multiplicative_edge_scaling, bool& valid, size_t *num_invalid=0);
      virtual Double metric(stk::mesh::Entity entity, bool& valid);
      virtual void update_node_positions( Double alpha);
      virtual bool check_convergence();

      Double line_search(bool& restarted, double mfac_mult=1.0);
      Double get_alpha_0();

      double nodal_metric(stk::mesh::Entity node, double alpha, double *coord_current, double *cg_d,  bool& valid );
      void nodal_gradient(stk::mesh::Entity node, double alpha, double *coord_current, double *cg_d,  bool& valid, double *ng);
      Double nodal_edge_length_ave(stk::mesh::Entity node);
      void get_edge_lengths(PerceptMesh * eMesh);
      void get_surface_normals(PerceptMesh * eMesh);
      void snap_nodes();

    };

  }

#endif
#endif
