// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/ReferenceMeshSmoother.hpp>
#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>

//#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stdio.h>

#include "mpi.h"

  namespace percept {

    void ReferenceMeshSmoother::sync_fields(int iter)
    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(m_eMesh->get_coordinates_field());
      fields.push_back(m_coord_field_projected);
      fields.push_back(m_coord_field_lagged);
      fields.push_back(m_coord_field_current);
      stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
    }

    void ReferenceMeshSmoother::project_all_delta_to_tangent_plane(stk::mesh::FieldBase *field)
    {
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );

      // node loop
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  std::pair<bool,int> fixed = this->get_fixed_flag(node);
                  if (fixed.first)
                    {
                      continue;
                    }

                  double *cg_g = m_eMesh->field_data(field, node);

                  if (fixed.second == MS_SURFACE)
                    {
                      project_delta_to_tangent_plane(node, cg_g);
                    }
                }
            }
        }
    }

    void ReferenceMeshSmoother::check_project_all_delta_to_tangent_plane(stk::mesh::FieldBase *field)
    {
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );

      int nc = m_eMesh->get_spatial_dim();

      // node loop
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  std::pair<bool,int> fixed = this->get_fixed_flag(node);
                  if (fixed.first)
                    {
                      continue;
                    }

                  double dproj[3];
                  double *cg_g = m_eMesh->field_data(field, node);
                  for (int ii=0; ii < nc; ++ii)
                    dproj[ii] = cg_g[ii];

                  if (fixed.second == MS_SURFACE)
                    {
                      double normal[3];
                      project_delta_to_tangent_plane(node, dproj, normal);
                      double dot=0.0, sc=0.0;
                      for (int ii=0; ii < nc; ++ii)
                        {
                          dot += dproj[ii]*normal[ii];
                          sc += dproj[ii]*dproj[ii];
                        }
                      sc = std::sqrt(sc);
                      if (sc > 1.e-10)
                        VERIFY_OP_ON(std::fabs(dot), <=, 1.e-6*sc, "bad norm");
                    }
                }
            }
        }
    }

    bool ReferenceMeshSmoother::check_convergence()
    {
      throw std::runtime_error("not implemented");
#if 0
      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & m_dmax ) );
      bool cond = (m_num_invalid == 0 && m_dmax < gradNorm);
      return cond;
#endif
      return false;
    }

    double ReferenceMeshSmoother::run_one_iteration()
    {
      throw std::runtime_error("not implemented");
      return 0.0;
    }

    static void print_comm_list( const stk::mesh::BulkData & mesh , bool doit )
    {
      using namespace stk::mesh;
      if ( doit ) {
        std::ostringstream msg ;

        msg << std::endl ;

        const std::vector<stk::mesh::Ghosting *> &ghostings = mesh.ghostings();
        std::vector<int> commProcs;
        stk::mesh::EntityVector entities;
        for(stk::topology::rank_t rank=stk::topology::NODE_RANK; rank<stk::topology::END_RANK; ++rank) {
          stk::mesh::get_entities(mesh, rank, entities);
          for(size_t e=0; e<entities.size(); ++e) {
            for(size_t g=0; g<ghostings.size(); ++g) {
              Entity entity = entities[e];
              EntityKey key = mesh.entity_key(entity);
              msg << "P" << mesh.parallel_rank() << ": " ;

              print_entity_key( msg , MetaData::get(mesh) , key );

              msg << " owner(" << mesh.parallel_owner_rank(entity) << ")" ;

              if ( Modified == mesh.state(entity) ) { msg << " mod" ; }
              else if ( Deleted == mesh.state(entity) ) { msg << " del" ; }
              else { msg << "    " ; }

              mesh.comm_procs(*ghostings[g], key, commProcs);
              for(size_t k=0; k<commProcs.size(); ++k) {
                msg << " gid, proc (" << g << "," << commProcs[k] << ")" ;
              }
              msg << std::endl ;
            }
          }
        }

        std::cout << msg.str();
      }
    }

    void ReferenceMeshSmoother::run_algorithm()
    {
      //std::cout << "\nP[" << m_eMesh->get_rank() << "] tmp srk ReferenceMeshSmoother innerIter= " << innerIter << " parallelIterations= " << parallelIterations << std::endl;

      //if (!m_eMesh->get_rank())
      //std::cout << "\nP[" << m_eMesh->get_rank() << "] tmp srk ReferenceMeshSmoother: running shape improver... \n" << std::endl;

      PerceptMesh *eMesh = m_eMesh;
      m_num_nodes = m_eMesh->get_number_nodes();

      print_comm_list(*eMesh->get_bulk_data(), false);

      stk::mesh::FieldBase *coord_field           = eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field;
      stk::mesh::FieldBase *coord_field_projected = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_N");
      stk::mesh::FieldBase *coord_field_original  = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_NM1");
      stk::mesh::FieldBase *coord_field_lagged    = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");

      m_coord_field_original  = coord_field_original;
      m_coord_field_projected = coord_field_projected;
      m_coord_field_lagged    = coord_field_lagged;
      m_coord_field_current   = coord_field_current;

      eMesh->copy_field(coord_field_lagged, coord_field_original);

      static int anim_step = 0;

      // untangle
      SmootherMetricUntangle untangle_metric(eMesh);
      SmootherMetricUntangleGen untangle_metric_gen(eMesh, this);

      // shape-size-orient smooth
      SmootherMetricShapeSizeOrient shape_size_orient_metric(eMesh);

      // shape
      SmootherMetricShapeB1 shape_b1_metric(eMesh, m_use_ref_mesh);
      SmootherMetricShapeB1Gen shape_b1_metric_gen(eMesh, this);

      SmootherMetricShapeC1 shape_c1_metric(eMesh);

      SmootherMetricShapeSizeOrientB1 shape_size_orient_b1_metric(eMesh);


      SmootherMetricShapeMeanRatio shape_mean_metric(eMesh, m_use_ref_mesh);

      // size
      SmootherMetricSizeB1 size_b1_metric(eMesh, m_use_ref_mesh);

      // laplace
      SmootherMetricLaplace laplace_metric(eMesh);

      // laplace weighted
      SmootherMetricLaplaceInverseVolumeWeighted laplace_weighted_metric(eMesh);

      // scaled jacobian
      SmootherMetricScaledJacobianElemental scaled_jac_metric(eMesh);

      // scaled jacobian - nodal
      SmootherMetricScaledJacobianNodal scaled_jac_metric_nodal(eMesh);

      //double omegas[] = {0.0, 0.001, 0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0};
      //double omegas[] = {0.001, 1.0};
      //double omegas[] = { 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.45,0.46,0.47,0.48,0.49,0.5,0.52,0.54,0.56,0.59, 0.6, 0.8, 1.0};
      double omegas[] = { 1.0};
      //double omegas[] = { 0.001, 0.01, 0.1, 0.2, 0.4, 1.0};
      //double omegas[] = {0.0, 0.001, 0.01, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.4, 0.6, 0.8, 1.0};
      int nomega = sizeof(omegas)/sizeof(omegas[0]);

      for (int outer = 0; outer < nomega; outer++)
        {
          double omega = (outer < nomega ? omegas[outer] : 1.0);
          m_omega = omega;
          m_omega_prev = omega;
          if (outer > 0) m_omega_prev = omegas[outer-1];

          // set current state and evaluate mesh validity (current = omega*project + (1-omega)*original)
          eMesh->nodal_field_axpbypgz(omega, coord_field_projected, (1.0-omega), coord_field_original, 0.0, coord_field_current);

          size_t num_invalid = MeshSmoother::parallel_count_invalid_elements(m_eMesh);

          m_num_invalid = num_invalid;
          m_untangled = (m_num_invalid == 0);

          if (0) eMesh->save_as("outer_iter_"+toString(outer)+"_m1_mesh.e");
          int iter_all=0;

          int do_anim = m_do_animation; // = frequency of anim writes
          if (do_anim)
            {
              std::ostringstream fileid_ss;
              fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);

              std::string oname = "anim_all.e";
              if (anim_step > 0) oname += "-s" + fileid_ss.str();
              eMesh->save_as(oname);
              ++anim_step;
            }

          int nstage = get_num_stages();
          //std::cout << "nstage= " << nstage << std::endl;
          for (int stage = 0; stage < nstage; stage++)
            {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
              if (m_meshGeometry && m_meshGeometry->geomKernel)
                {
                  m_meshGeometry->geomKernel->m_useFoundFaceMap = true;
                  //m_meshGeometry->geomKernel->m_nodeToFoundFaceMap.clear();
                }
#endif
              m_stage = stage;
              if (stage==0)
                {
                  if (eMesh->getProperty("smoother_metric_type") == "generated")
                    {
                      m_metric = &untangle_metric_gen;
                      if (!eMesh->get_rank()) std::cout << "using generated for untangle" << std::endl;
                    }
                  else
                    m_metric = &untangle_metric;
                  //m_metric = &scaled_jac_metric_nodal;
                }
              else
                {
                  int num_invalid_1 = MeshSmoother::parallel_count_invalid_elements(m_eMesh);
                  VERIFY_OP_ON(num_invalid_1, ==, 0, "Invalid elements exist for start of stage 2, quiting");

                  //m_metric = &shape_size_orient_metric;
                  ////m_metric = &shape_size_orient_b1_metric;
                  if (eMesh->getProperty("smoother_metric_type") == "generated")
                    {
                      m_metric = &shape_b1_metric_gen;
                      if (!eMesh->get_rank()) std::cout << "using generated for B1" << std::endl;
                    }
                  else
                    {
                      m_metric = &shape_b1_metric;
                      if (eMesh->getProperty("smoother_metric_type") == "size")
                        {
                          if (eMesh->get_rank() == 0) std::cout << "INFO: using SmootherMetricSizeB1 metric" << std::endl;
                          m_metric = &size_b1_metric;
                        }
                      if (eMesh->getProperty("smoother_metric_type") == "mean_ratio")
                        {
                          if (eMesh->get_rank() == 0) std::cout << "INFO: using SmootherMetricShapeMeanRatio metric" << std::endl;
                          m_metric = &shape_mean_metric;
                        }

                      //m_metric = &shape_c1_metric;
                      //m_metric = &laplace_metric;
                      //m_metric = &laplace_weighted_metric;
                      //m_metric = &scaled_jac_metric;
                    }
                }

              for (int iter = 0; iter < innerIter; ++iter, ++iter_all)
                {
                  m_iter = iter;
                  m_num_invalid = MeshSmoother::parallel_count_invalid_elements(m_eMesh);
                  int num_invalid_0 = m_num_invalid;
                  if (!m_untangled && m_num_invalid == 0)
                    {
                      m_untangled = true;
                    }

                  //               if (!m_eMesh->get_rank() && num_invalid_0)
                  //                 std::cout << "\ntmp srk ReferenceMeshSmoother num_invalid current= " << num_invalid_0
                  //                           << (num_invalid ? " WARNING: invalid elements exist before smoothing" : "OK")
                  //                           << std::endl;

                  m_global_metric = run_one_iteration();

                  sync_fields(iter);
                  //num_invalid_0 = MeshSmoother::parallel_count_invalid_elements(m_eMesh);
                  //m_num_invalid = num_invalid_0;
                  bool conv = check_convergence();
                  if (!m_eMesh->get_rank())
                  {
                    std::cout << "P[" << m_eMesh->get_rank() << "] " << "INFO: iter= " << iter << " dmax= " << m_dmax << " dmax_rel= " << m_dmax_relative << " grad/g0= " << std::sqrt(m_dnew/m_d0)
                              << " gradScaled= " << m_grad_norm_scaled
                              << " m_alpha= " << m_alpha  << " m_alpha_0 = " << m_alpha_0
                              << " num_invalid= " << num_invalid_0
                              << " m_global_metric= " << m_global_metric << " stage= " << stage << " m_untangled= " << m_untangled
                              << std::endl;
                  }

                  if (do_anim)
                    {
                      //eMesh->save_as("iter_"+toString(outer)+"_"+toString(stage)+"."+toString(iter+1)+".e");
                       if (iter_all % do_anim == 0)
                         {
                           std::ostringstream fileid_ss;
                           fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);

                           std::string oname = "anim_all.e";
                           if (anim_step > 0) oname += "-s" + fileid_ss.str();
                           eMesh->save_as(oname);
                           ++anim_step;
                         }
                    }

                  if (!m_untangled && m_num_invalid == 0)
                    {
                      m_untangled = true;
                    }
                  if (conv && m_untangled)
                    {
                      //std::cout << "P[" << m_eMesh->get_rank() << "] tmp srk a1 converged, stage= " << m_stage << std::endl;
                      break;
                    }
                }

              if (0) eMesh->save_as("outer_iter_"+toString(outer)+"_"+toString(stage)+"_mesh.e");
            }

          eMesh->copy_field(coord_field_lagged, coord_field);

        }

      //if (!m_eMesh->get_rank())

      MPI_Barrier( MPI_COMM_WORLD );

      if (!m_eMesh->get_rank())
        std::cout << "\nP[" << m_eMesh->get_rank() << "] tmp srk ReferenceMeshSmoother: running shape improver... done \n" << std::endl;

    }

  }

#endif
