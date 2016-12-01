// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmoother3_hpp
#define ReferenceMeshSmoother3_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT) && ENABLE_SMOOTHER3

#include <percept/mesh/mod/smoother/ReferenceMeshSmoother1.hpp>
#include <percept/mesh/mod/smoother/TpetraLinearSolver.hpp>
#include <TpetraBaseLinearSystem.h>


  namespace percept {

    typedef tftk::linsys::TpetraCrsLinearSystemImpl<tftk::linsys::MeshManagerNonFieldBased> TpetraBaseLinearSystem;
    class TpetraLinearSystem : public TpetraBaseLinearSystem {
    public:
      TpetraLinearSystem(
                         const unsigned numDof,
                         const std::string & name,
                         stk::mesh::BulkData & bulk_data,
                         std::shared_ptr<tftk::linsys::MeshManagerNonFieldBased> & mesh_manager,
                         const tftk::linsys::GlobalIdField * idField,
                         stk::topology::rank_t algorithm_topology_rank = stk::topology::NODE_RANK)
        : TpetraBaseLinearSystem (numDof, name, bulk_data, idField, {}, nullptr, mesh_manager, algorithm_topology_rank) {}

      void applyDirichletBCs(percept::PerceptMesh *eMesh,
                             stk::mesh::FieldBase * rhsField,
                             stk::mesh::FieldBase * solutionField,
                             const unsigned beginPos,
                             const unsigned endPos);


    };

    class ReferenceMeshSmoother3 : public ReferenceMeshSmoother1 {

    public:

      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length
      ReferenceMeshSmoother3(PerceptMesh *eMesh,
                            stk::mesh::Selector *boundary_selector=0,
                            MeshGeometry *meshGeometry=0,
                            int inner_iterations = 100,
                            double grad_norm =1.e-8,
                            int parallel_iterations = 20)
        : ReferenceMeshSmoother1(eMesh, boundary_selector, meshGeometry, inner_iterations, grad_norm, parallel_iterations),
          m_isSetup(false)
      {}

      ~ReferenceMeshSmoother3();

    protected:
      bool m_isSetup;

      virtual void get_step();
      virtual double run_one_iteration();
      virtual void setup_linsys();
      double local_line_search(stk::mesh::Entity node, stk::mesh::FieldBase *cg_g_field);

      typedef tftk::linsys::MeshManagerNonFieldBased TpetraMeshManager;
      TpetraLinearSystem * m_linearSystem;
      TpetraLinearSolver * m_linearSolver;
      std::shared_ptr<tftk::linsys::MeshManagerNonFieldBased> m_meshManager;
    };

  }

#endif
#endif
