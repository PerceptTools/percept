// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/geometry/kernel/GeometryKernel.hpp>

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherBase.hpp>
#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stdio.h>

#include "mpi.h"

  namespace percept {

  template <typename MeshType>
  ReferenceMeshSmootherBaseImpl<MeshType>::
  ReferenceMeshSmootherBaseImpl(PerceptMesh *eMesh,
                            typename MeshType::MTSelector *boundary_selector,
                            typename MeshType::MTMeshGeometry *meshGeometry,
                            int inner_iterations,
                            double grad_norm,
                            int parallel_iterations)
    : MeshSmootherImpl<MeshType>(eMesh, boundary_selector, meshGeometry, inner_iterations, grad_norm, parallel_iterations),
      m_scale(0),
      m_dmax(0),
      m_dmax_relative(0),
      m_dnew(0), m_dold(0), m_d0(1), m_dmid(0), m_dd(0), m_alpha(0), m_alpha_0(0), m_grad_norm(0), m_grad_norm_scaled(0),
      m_total_metric(0),
      m_stage(0),
      m_omega(0),
      m_omega_prev(0),
      m_iter(0),
      m_num_invalid(0), m_global_metric(std::numeric_limits<double>::max()), m_untangled(false), m_num_nodes(0),
      m_use_ref_mesh(true), m_do_animation(0)
  {
    if (eMesh->getProperty("ReferenceMeshSmootherBase_do_anim") != "")
      {
        m_do_animation = toInt(eMesh->getProperty("ReferenceMeshSmootherBase_do_anim"));
        std::cout << "m_do_animation= " << m_do_animation << std::endl;
      }
  }

  template<>
    void ReferenceMeshSmootherBaseImpl<STKMesh>::sync_fields(int iter)
    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(m_eMesh->get_coordinates_field());
      fields.push_back(m_coord_field_projected);
      fields.push_back(m_coord_field_lagged);
      fields.push_back(m_coord_field_current);
      stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
    }

    template<>
    void ReferenceMeshSmootherBaseImpl<StructuredGrid>::sync_fields(int iter)
    {}

    template<>
    void ReferenceMeshSmootherBaseImpl<STKMesh>::project_all_delta_to_tangent_plane(stk::mesh::FieldBase *field)
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

    template<>
    void ReferenceMeshSmootherBaseImpl<STKMesh>::check_project_all_delta_to_tangent_plane(stk::mesh::FieldBase *field)
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

    template<typename MeshType>
    bool ReferenceMeshSmootherBaseImpl<MeshType>::check_convergence()
    {
      throw std::runtime_error("not implemented");
#if 0
      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & m_dmax ) );
      bool cond = (m_num_invalid == 0 && m_dmax < gradNorm);
      return cond;
#endif
      return false;
    }

    template<typename MeshType>
    double ReferenceMeshSmootherBaseImpl<MeshType>::run_one_iteration()
    {
      throw std::runtime_error("not implemented");
      return 0.0;
    }

    template<typename MeshType>
    struct GA_run_algorithm
    {
      using This = GA_run_algorithm<MeshType>;
      using RefMeshSmootherBase = ReferenceMeshSmootherBaseImpl<MeshType>;
      RefMeshSmootherBase *m_rms;
      PerceptMesh *m_eMesh;
      typename MeshType::MTField *coord_field;
      typename MeshType::MTField *coord_field_current;
      typename MeshType::MTField *coord_field_projected;
      typename MeshType::MTField *coord_field_original;
      typename MeshType::MTField *coord_field_lagged;

      GA_run_algorithm(RefMeshSmootherBase *rms, PerceptMesh *eMesh);

    };

    template<>
    GA_run_algorithm<STKMesh>::GA_run_algorithm(RefMeshSmootherBase *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
    {
      coord_field           = eMesh->get_coordinates_field();
      coord_field_current   = coord_field;
      coord_field_projected = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_N");
      coord_field_original  = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_NM1");
      coord_field_lagged    = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");

      m_rms->m_coord_field_original  = coord_field_original;
      m_rms->m_coord_field_projected = coord_field_projected;
      m_rms->m_coord_field_lagged    = coord_field_lagged;
      m_rms->m_coord_field_current   = coord_field_current;
    }

    template<>
    GA_run_algorithm<StructuredGrid>::GA_run_algorithm(RefMeshSmootherBase *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
    {
      std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
      coord_field                              = bsg->m_fields["coordinates"].get();
      coord_field_current                      = coord_field;
      coord_field_projected                    = bsg->m_fields["coordinates_N"].get();
      coord_field_original                     = bsg->m_fields["coordinates_NM1"].get();
      coord_field_lagged                       = bsg->m_fields["coordinates_lagged"].get();

      m_rms->m_coord_field_original  = coord_field_original;
      m_rms->m_coord_field_projected = coord_field_projected;
      m_rms->m_coord_field_lagged    = coord_field_lagged;
      m_rms->m_coord_field_current   = coord_field_current;
    }

    template<typename MeshType>
    void anim(PerceptMesh *eMesh, int& anim_step)
    {
    }

    template<>
    void anim<STKMesh>(PerceptMesh *eMesh, int& anim_step)
    {
      std::ostringstream fileid_ss;
      fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);

      std::string oname = "anim_all.e";
      if (anim_step > 0) oname += "-s" + fileid_ss.str();
      eMesh->save_as(oname);
      ++anim_step;
    }

    template<>
    void anim<StructuredGrid>(PerceptMesh *eMesh, int& anim_step)
    {
      std::ostringstream fileid_ss;
      fileid_ss << "."  << anim_step;

      std::string oname = "anim_all"+fileid_ss.str();
      eMesh->get_block_structured_grid()->dump_vtk(oname);
      ++anim_step;
    }

    template<typename MeshType>
    void set_geometry_info(typename MeshType::MTMeshGeometry *meshGeom)
    {
    }

    template<>
    void set_geometry_info<STKMesh>(typename STKMesh::MTMeshGeometry *meshGeometry)
    {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
      if (meshGeometry && meshGeometry->geomKernel)
        {
          meshGeometry->geomKernel->m_useFoundFaceMap = true;
        }
#endif
    }

    template<typename MeshType>
    void ReferenceMeshSmootherBaseImpl<MeshType>::run_algorithm()
    {
      GA_run_algorithm<MeshType> ga(this, Base::m_eMesh);

      typename MeshType::MTField *coord_field           = ga.coord_field;
      typename MeshType::MTField *coord_field_current   = ga.coord_field_current;
      typename MeshType::MTField *coord_field_projected = ga.coord_field_projected;
      typename MeshType::MTField *coord_field_original  = ga.coord_field_original;
      typename MeshType::MTField *coord_field_lagged    = ga.coord_field_lagged;

      PerceptMesh *eMesh = Base::m_eMesh;
      m_num_nodes = eMesh->get_number_nodes();

      eMesh->copy_field(coord_field_lagged, coord_field_original);

      static int anim_step = 0;

      // untangle
      SmootherMetricUntangleImpl<MeshType> untangle_metric(eMesh);

      // shape
      SmootherMetricShapeB1Impl<MeshType> shape_b1_metric(eMesh, m_use_ref_mesh);

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

          size_t num_invalid = MeshSmootherImpl<MeshType>::parallel_count_invalid_elements(eMesh);

          m_num_invalid = num_invalid;
          m_untangled = (m_num_invalid == 0);

          int iter_all=0;

          int do_anim = m_do_animation; // = frequency of anim writes
          if (do_anim)
            {
              anim<MeshType>(eMesh, anim_step);
            }

          int nstage = get_num_stages();
          for (int stage = 0; stage < nstage; stage++)
            {
              set_geometry_info<MeshType>(Base::m_meshGeometry);

              m_stage = stage;
              if (stage==0)
                {
                  m_metric = &untangle_metric;
                }
              else
                {
                  int num_invalid_1 = MeshSmootherImpl<MeshType>::parallel_count_invalid_elements(eMesh);
                  VERIFY_OP_ON(num_invalid_1, ==, 0, "Invalid elements exist for start of stage 2, quiting");

                  m_metric = &shape_b1_metric;
                }

              for (int iter = 0; iter < Base::innerIter; ++iter, ++iter_all)
                {
                  m_iter = iter;
                  m_num_invalid = MeshSmootherImpl<MeshType>::parallel_count_invalid_elements(eMesh);
                  int num_invalid_0 = m_num_invalid;
                  if (!m_untangled && m_num_invalid == 0)
                    {
                      m_untangled = true;
                    }

                  m_global_metric = run_one_iteration();

                  sync_fields(iter);
                  bool conv = check_convergence();
                  if (!eMesh->get_rank())
                    {
                      std::cout << " iter= " << iter << " dmax= " << m_dmax << " dmax_rel= " << m_dmax_relative << " grad/g0= " << std::sqrt(m_dnew/m_d0)
                                << " gradScaled= " << m_grad_norm_scaled
                                << " m_alpha= " << m_alpha  << " m_alpha_0 = " << m_alpha_0
                                << " num_invalid= " << num_invalid_0
                                << " m_global_metric= " << m_global_metric << " stage= " << stage << " m_untangled= " << m_untangled
                                << std::endl;
                    }

                  if (do_anim)
                    {
                      if (iter_all % do_anim == 0)
                        {
                          anim<MeshType>(eMesh, anim_step);
                        }
                    }

                  if (!m_untangled && m_num_invalid == 0)
                    {
                      m_untangled = true;
                    }
                  if (conv && m_untangled)
                    {
                      break;
                    }
                }
            }

          eMesh->copy_field(coord_field_lagged, coord_field);
        }

      if (!eMesh->get_rank())
        std::cout << "\nP[" << eMesh->get_rank() << "] tmp srk ReferenceMeshSmootherBase: running shape improver... done \n" << std::endl;

    }

    template class ReferenceMeshSmootherBaseImpl<STKMesh>;
    template class ReferenceMeshSmootherBaseImpl<StructuredGrid>;
  }

#endif
