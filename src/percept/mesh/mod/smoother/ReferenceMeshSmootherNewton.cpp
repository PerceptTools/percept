// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT) && ENABLE_SMOOTHER3


#include <percept/mesh/mod/smoother/ReferenceMeshSmootherNewton.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/SmootherMetricUntangleGen.hpp>
#include <percept/mesh/mod/smoother/SmootherMetricShapeB1Gen.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <percept/math/DenseMatrix.hpp>
#include <percept/math/Math.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stdio.h>
#include <limits>

#include "mpi.h"
#include <cstdio>

#define DEBUG_PRINT 0
#define PRINT(a) do { if (DEBUG_PRINT && !m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_1(a) do { if (!m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_2(a) do {  std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)

namespace percept {

  ReferenceMeshSmootherNewton::~ReferenceMeshSmootherNewton()
  {
    delete m_linearSystem;
    delete m_linearSolver;
  }

  void ReferenceMeshSmootherNewton::setup_linsys()
  {
    if (m_isSetup)
      return;

    m_isSetup = true;

    PerceptMesh *eMesh = m_eMesh;

    tftk::linsys::GlobalIdField * cg_gid_field =
      eMesh->get_fem_meta_data()->get_field<tftk::linsys::GlobalIdField>(eMesh->node_rank(), "cg_gid");
    VERIFY_OP_ON(cg_gid_field, !=, 0, "cg_gid_field");

    int spatialDim = eMesh->get_spatial_dim();
    int extraLambdaDof = 0;
    if (m_eMesh->get_smooth_surfaces() && m_eMesh->getProperty("ReferenceMeshSmootherNewton.use_lambda") == "true")
      extraLambdaDof = 1;

    m_meshManager = std::make_shared<TpetraMeshManager>(unsigned(spatialDim + extraLambdaDof), *m_eMesh->get_bulk_data(), *cg_gid_field);

    // FIXME - only need extra dof on the boundaries, but currently tftk/linsys doesn't support this
    m_linearSystem = new TpetraLinearSystem(unsigned(spatialDim + extraLambdaDof), "test", *m_eMesh->get_bulk_data(), m_meshManager, cg_gid_field);
    m_linearSolver = new TpetraLinearSolver(*m_linearSystem, m_meshManager);

    stk::mesh::Selector selu = m_eMesh->get_fem_meta_data()->universal_part();

    {
      const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          const unsigned num_nodes_in_bucket = bucket.size();
          for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
            {
              stk::mesh::Entity node = bucket[i_node];
              stk::mesh::EntityId *entityId = stk::mesh::field_data(*cg_gid_field, node);
              *entityId = eMesh->id(node);
            }
        }
    }

    m_meshManager->computeMeshLocalIds(*cg_gid_field, selu);
    m_linearSystem->beginLinearSystemConstruction();

    if (1)
      {
        // element loop: compute deltas
        const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (MeshSmoother::select_bucket(**k, m_eMesh))
              {
                stk::mesh::Bucket & bucket = **k ;

                const unsigned num_elements_in_bucket = bucket.size();
                if (!bucket.owned())
                  continue;

                static std::vector<stk::mesh::Entity> nodes;

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];
                    if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
                      continue;

                    stk::mesh::Entity const* elem_nodes = m_eMesh->get_bulk_data()->begin_nodes(element);
                    unsigned num_node = m_eMesh->get_bulk_data()->num_nodes(element);

                    nodes.resize(num_node);

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        stk::mesh::Entity node = elem_nodes[ inode ];
                        nodes[inode] = node;
                      }
                    m_meshManager->addSymmetricConnections(nodes);
                  }
              }
          }
      }
    m_linearSystem->finalizeLinearSystem();
    bool doPrecond = true;
    if (m_stage == 0)
      doPrecond = false;
    m_linearSolver->setupLinearSolver(doPrecond);
  }

  void
  TpetraLinearSystem::applyDirichletBCs(percept::PerceptMesh *eMesh,
                                        stk::mesh::FieldBase * rhsField,
                                        stk::mesh::FieldBase * solutionField,
                                        const unsigned beginPos,
                                        const unsigned endPos)
  {
    std::ofstream *file = 0;
    if (0)
      {
        file = new std::ofstream("mat.mm");
        *file << "{\n";
      }
    stk::mesh::MetaData& stkmesh_meta_data = *eMesh->get_fem_meta_data();
    stk::mesh::BulkData& stkmesh_bulk_data = *eMesh->get_bulk_data();
    unsigned numDof_ = unsigned(eMesh->get_spatial_dim());

    const stk::mesh::Selector selector = stkmesh_meta_data.locally_owned_part();

    stk::mesh::BucketVector const& buckets =
      stkmesh_bulk_data.get_buckets( stk::topology::NODE_RANK, selector );

    int nentity = 0;
    for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
          ib != buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;

      const unsigned fieldSize = field_bytes_per_entity(*solutionField, b) / sizeof(double);
      ThrowRequire(fieldSize == numDof_);

      if (!b.owned())
        continue;

      const stk::mesh::Bucket::size_type length   = b.size();
      const double * solution = (double*)stk::mesh::field_data(*solutionField, b);
      const double * rhs_data = (double *)stk::mesh::field_data(*rhsField, b);
      const double bcValues[3] = {0.0, 0.0, 0.0};

      typedef tftk::linsys::LocalOrdinal LocalOrdinal;
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const double> values;
      std::vector<LocalOrdinal> diagonalIds(1);
      std::vector<double> new_values;
      std::vector<double> diagonal_values(1);

      for (stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const stk::mesh::Entity entity = b[k];

        const LocalOrdinal localIdOffset = stkmesh_bulk_data.local_id(entity)*numDof_;
        ThrowAssert(!(localIdOffset < 0));

        ThrowAssert(localIdOffset < this->numOwnedRows_);
        double rn=0.0;
        for(unsigned d = beginPos; d < endPos; ++d)
          {
            double r = rhs_data[k*fieldSize + d];
            rn += r*r;
          }
        if (1)
          {
            for(unsigned d  =beginPos; d < endPos; ++d, ++nentity)
              {
                const LocalOrdinal localId = localIdOffset + d;

                diagonal_values[0] = 1.0;
                diagonalIds[0] = localId;

                this->ownedMatrix_->getLocalRowView(localId, indices, values);
                const size_t rowLength = values.size();
                new_values.resize(rowLength);
                double vn = 0.0, vs = 0.0;
                for(size_t i=0; i < rowLength; ++i) {
                  new_values[i] = 0;
                  vn += std::fabs(values[i]);
                  vs += values[i];
                }

                if (rn == 0.0)
                  {
                    this->ownedMatrix_->replaceLocalValues(localId, indices, new_values);
                    this->ownedMatrix_->replaceLocalValues(localId, diagonalIds, diagonal_values);

                    const double bc_residual = (bcValues[d] - 0.0*solution[k*fieldSize + d]);
                    this->ownedRhs_->replaceLocalValue(localId, bc_residual);
                  }

                if (1)
                  {
                    this->ownedMatrix_->getLocalRowView(localId, indices, values);
                    const size_t RowLength = values.size();
                    //double vn = 0.0, vs = 0.0;
                    if (file) *file << (nentity ? ",{" : "{");
                    for(size_t i=0; i < RowLength; ++i) {
                      if (file) *file << Util::convert_to_mm(values[i],20) << (i == RowLength-1 ? " " : ",");
                      new_values[i] = 0;
                      //vn += std::fabs(values[i]);
                      //vs += values[i];
                    }
                    if (file) *file << "}" << "(* rn= " << Util::convert_to_mm(rn, 20) << "*)\n";
                  }

              }
          }
      }
    }
    if (file) *file << "}\n";
    if (file) delete file;
  }

  void ReferenceMeshSmootherNewton::get_step()
  {
    double t0 = stk::cpu_time();

    PerceptMesh *eMesh = m_eMesh;
    stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
    stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    stk::mesh::FieldBase *cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
    stk::mesh::FieldBase *cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");
    stk::mesh::FieldBase *cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    stk::mesh::FieldBase *cg_lambda_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_lambda");
    //stk::mesh::FieldBase *cg_normal_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_normal");

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    m_scale = 1.e-10;

    // g=0
    eMesh->nodal_field_set_value(cg_g_field, 0.0);
    eMesh->nodal_field_set_value(cg_r_field, 0.0);

    SmootherMetricUntangle smu(eMesh);

    m_linearSystem->zeroSystem();

    if (1)
      {
        SmootherMetricElementFunction smf(m_eMesh, m_metric);
        typedef FDGradient<SmootherMetricElementFunction> Gradient;
#define USE_FD 1
#define USE_FDHessianFromGrad 0
#define USE_ANALYTIC !USE_FD
        typedef FDHessian<SmootherMetricElementFunction> Hessian;
        typedef FDHessianFromGrad<SmootherMetricElementFunction> HessianFromGrad;

#if USE_FD
        Hessian hessian;
#elif USE_FDHessianFromGrad
        HessianFromGrad hessian(smf);
#elif USE_ANALYTIC
        //Hessian hessian;
#endif
        Gradient grad;
        (void)grad;

        // grad.test();
        // hessian.test();

        int extraLambdaDof = 0;

        if (m_eMesh->get_smooth_surfaces() && m_eMesh->getProperty("ReferenceMeshSmootherNewton.use_lambda") == "true")
          {
            VERIFY_OP_ON(cg_lambda_field, !=, 0, "bad lambda field");
            extraLambdaDof = 1;
          }

        // element loop: compute deltas
        const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (MeshSmoother::select_bucket(**k, m_eMesh))
              {
                stk::mesh::Bucket & bucket = **k ;

                const unsigned num_elements_in_bucket = bucket.size();
                m_metric->m_topology_data = m_eMesh->get_cell_topology(bucket);

                shards::CellTopology topo(m_metric->m_topology_data);
                int nnodes = m_metric->m_topology_data->node_count;
                int nc = eMesh->get_spatial_dim();
                int ndof = extraLambdaDof + nc;
                int N = nnodes*ndof;
                std::vector<double> G(N,0.0), H(N*N,0.0), U(N,0.0);

                static std::vector<stk::mesh::Entity> nodes;

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];
                    if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
                      continue;

                    smf.set_element(element);
                    VERIFY_OP_ON(N, ==, smf.get_n(), "bad N");

                    const MyPairIterRelation elem_nodes(*m_eMesh, element, stk::topology::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    nodes.resize(num_node);
                    std::vector<double> num_elem(num_node,0.0);
                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        stk::mesh::Entity node = elem_nodes[ inode ].entity();

                        num_elem[inode] = (double)m_eMesh->get_bulk_data()->num_elements(node);

                        nodes[inode] = node;

                        double *coord = m_eMesh->field_data(coord_field, node);
                        for (int jc=0; jc < nc; ++jc)
                          U[inode*ndof + jc] = coord[jc];
                        if (extraLambdaDof)
                          {
                            double *lambda = m_eMesh->field_data(cg_lambda_field, node);
                            VERIFY_OP_ON(lambda, !=, 0, "bad lambda field");
                            U[inode*ndof + nc] = lambda[0];
                          }
                      }

#if USE_ANALYTIC
                    bool valid = true;
                    typedef double GradType[8][4];
                    typedef double HessType[8][4][8][4];
                    double lgrad[8][4], lhess[8][4][8][4];
                    m_metric->grad_and_hessian(element, valid, lgrad, lhess, true, true);
                    for (unsigned inode = 0; inode < num_node; inode++)
                      {
                        for (int jc = 0; jc < ndof; ++jc)
                          {
                            G[(inode*ndof + jc)] = lgrad[inode][jc];

                            for (unsigned jnode = 0; jnode < num_node; jnode++)
                              {
                                for (int kc = 0; kc < ndof; ++kc)
                                  {
                                    H[(inode*ndof + jc)*N + (jnode*ndof + kc)] = lhess[inode][jc][jnode][kc];
                                  }
                              }
                          }
                      }
#else
                    if (USE_FD)
                      grad(smf, &U[0], &G[0]);
                    else
                      smf.analytic_gradient(&U[0], &G[0]);
                    if (1)
                      {
                        hessian(smf, &U[0], &H[0]);
                      }
                    else
                      {
                        for (unsigned inode=0; inode < num_node; inode++)
                          {
                            for (int jc=0; jc < ndof; ++jc)
                              {
                                int ii = inode*ndof + jc;
                                for (int jj = 0; jj < N; ++jj)
                                  H[ii*N + jj] = 0.0;
                                H[ii*N + ii] = 1.0/num_elem[inode];
                              }
                          }
                      }
#endif

                    // Dirichlet
                    if (1)
                      {
                        for (unsigned inode=0; inode < num_node; inode++)
                          {
                            stk::mesh::Entity node = elem_nodes[ inode ].entity();

                            std::pair<bool,int> fixed = this->get_fixed_flag(node);

                            if (fixed.first || eMesh->aura(node))
                              {
                                for (int jc=0; jc < ndof; ++jc)
                                  {
                                    int ii = inode*ndof + jc;
                                    G[ii] = 0.0;
                                    if (1)
                                      {
                                        for (int jj = 0; jj < N; ++jj)
                                          H[ii*N + jj] = 0.0;
                                        H[ii*N + ii] = 1.0/num_elem[inode];
                                      }
                                  }
                              }

                            else if (!fixed.first && (fixed.second == MS_SURFACE || fixed.second == MS_ON_BOUNDARY))
                              {
                                double lhs[3][3], rhs[3];
                                for (int jc=0; jc < nc; ++jc)
                                  {
                                    int ii = inode*ndof + jc;
                                    rhs[jc] = G[ii];
                                    for (int kc=0; kc < nc; ++kc)
                                      {
                                        int jj = inode*ndof + kc;
                                        lhs[jc][kc] = H[ii*N + jj];
                                      }
                                  }
                                // FIXME
                                if (1)
                                  {
                                    //double norm[3];
                                    enforce_tangent_plane(node, rhs, lhs);
                                  }
                                for (int jc=0; jc < nc; ++jc)
                                  {
                                    int ii = inode*ndof + jc;
                                    G[ii] = rhs[jc];
                                    for (int kc=0; kc < nc; ++kc)
                                      {
                                        int jj = inode*ndof + kc;
                                        H[ii*N + jj] = lhs[jc][kc];
                                      }
                                  }
                              }

                          }
                      }

                    if (0)
                      {
                        Teuchos::Array<int> dim(2);
                        dim[0] = N; dim[1] = N;
                        MDArray HM(dim, &H[0], true, true);
                        //std::cout << "HM=\n" << printForMathematica(HM) << std::endl;
                        std::cout << "HM=\n" << HM << std::endl;
                      }



                    m_linearSystem->sumInto(nodes, G, H);
                  }
              }
          }
      }

    double t_1_form = stk::cpu_time();

    m_linearSystem->loadComplete();

    if (m_stage <= 1)
      {
        m_meshManager->copy_tpetra_to_stk(m_linearSystem->vector(), cg_g_field);

        /// r = -g
        eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);

        //m_linearSystem->applyDirichletBCs(m_eMesh, cg_r_field, cg_r_field, 0, unsigned(m_eMesh->get_spatial_dim()));
      }

    double t_2_load = stk::cpu_time();

    //if (m_stage == 0) m_linearSystem->writeToFile("test");

    m_linearSolver->solve(cg_g_field);
    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(cg_g_field);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
      stk::mesh::copy_owned_to_shared(*eMesh->get_bulk_data(), fields);
    }

    double t_3_solve = stk::cpu_time();

    if (m_stage <= 1)
      {
        double dn = eMesh->nodal_field_dot(cg_r_field, cg_g_field);
        if (m_eMesh->get_rank() == 0 && dn > 0) std::cout << "WARNING, not a descent direction, should be < 0, is= " << dn << std::endl;
      }

    // project deltas to surface
    if (m_eMesh->get_smooth_surfaces())
      {
        m_eMesh->copy_field(cg_s_field, cg_g_field);
        check_project_all_delta_to_tangent_plane(cg_g_field);
        if (0) project_all_delta_to_tangent_plane(cg_g_field);
      }

    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(cg_g_field);
      fields.push_back(cg_r_field);
      fields.push_back(cg_d_field);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
      stk::mesh::copy_owned_to_shared(*eMesh->get_bulk_data(), fields);
    }

    if (m_eMesh->get_smooth_surfaces()) // FIXME
      {
        double dn = eMesh->nodal_field_dot(cg_r_field, cg_g_field);
        double dnNew = 0.0;
        if (dn > 0)
          {
            // copy field that hasn't been projected
            m_eMesh->copy_field(cg_g_field, cg_s_field);
            dnNew = eMesh->nodal_field_dot(cg_r_field, cg_g_field);
            if (m_eMesh->get_rank() == 0 && dn > 0) std::cout << "WARNING, not a descent direction, after project_delta_to_tangent_plane: should be < 0, is= "
                                                              << dn << " dnNew= " << dnNew << std::endl;
          }
        VERIFY_OP_ON(dnNew, <=, 0.0, "bad dnNew");
      }

    double t_4_comm = stk::cpu_time();
    if (m_eMesh->get_rank() == 0)
      {
        std::cout << "\ntimes= form= " << (t_1_form - t0) << " loadComplete= " << (t_2_load - t_1_form)
                  << " solve= " << (t_3_solve - t_2_load) << " comm= " << (t_4_comm - t_3_solve) << std::endl;
      }
  }

  double ReferenceMeshSmootherNewton::run_one_iteration()
  {
    PerceptMesh *eMesh = m_eMesh;

    stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    stk::mesh::FieldBase *cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
    stk::mesh::FieldBase *cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");
    stk::mesh::FieldBase *cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    stk::mesh::FieldBase *coordinates_0 = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_0");
    stk::mesh::FieldBase *cg_lambda_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_lambda");

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    bool total_valid=false;

    if (m_iter == 0)
      {
        get_edge_lengths(m_eMesh);
        eMesh->nodal_field_set_value(cg_g_field, 0.0);
        eMesh->nodal_field_set_value(cg_r_field, 0.0);
        eMesh->nodal_field_set_value(cg_d_field, 0.0);
        eMesh->nodal_field_set_value(cg_s_field, 0.0);
        if (cg_lambda_field) eMesh->nodal_field_set_value(cg_lambda_field, 0.0);
        setup_linsys();
      }

    // this field is used for smoothing surfaces
    if (eMesh->get_smooth_surfaces())
      {
        eMesh->copy_field(coordinates_0, m_coord_field_original);
        get_surface_normals(eMesh);
      }

    get_step();

    eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);
    eMesh->copy_field(cg_s_field, cg_r_field);

    m_dnew = eMesh->nodal_field_dot(cg_s_field, cg_s_field);
    if (m_iter == 0)
      {
        m_d0 = m_dnew;
      }

    double metric_check = total_metric( 0.0, 1.0, total_valid);
    m_total_metric = metric_check;
    //PRINT_1( "INFO: tmp srk m_iter= " << m_iter << " m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " metric_check= " << metric_check );
    if (check_convergence() || metric_check == 0.0)
      {
        PRINT_1( "INFO: tmp srk already converged m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " metric_check= " << metric_check );
        //update_node_positions
        return total_metric(0.0,1.0, total_valid);
      }
    if (m_stage)
      VERIFY_OP_ON(total_valid, ==, true, "bad mesh before line search");

#if 1
    bool restarted = false;
    double alpha = line_search(restarted, 0.0);
#else
    double alpha = 0.1;
#endif
    double snorm = eMesh->nodal_field_dot(cg_s_field, cg_s_field);
    m_grad_norm_scaled = m_alpha_0*std::sqrt(snorm)/double(m_num_nodes);

    /// x = x + alpha*d
    m_alpha = alpha;
    update_node_positions(alpha);

    // check if re-snapped geometry is acceptable
    if (m_eMesh->get_smooth_surfaces())
      {
        // FIXME
        snap_nodes();
        if (m_stage != 0)
          {
            bool total_valid_0=true;
            total_metric( 0.0, 1.0, total_valid_0);
            VERIFY_OP_ON(total_valid_0, ==, true, "bad mesh after snap_node_positions...");
          }
      }

    double tm = total_metric(0.0,1.0, total_valid);
    //PRINT_1( "tmp srk iter= "<< m_iter << " dmax= " << m_dmax << " alpha= " << alpha << " global metric= " << tm);

    return tm;
  }


}


#endif
