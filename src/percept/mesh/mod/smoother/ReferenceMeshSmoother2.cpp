// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)


#include <percept/mesh/mod/smoother/ReferenceMeshSmoother2.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <percept/math/DenseMatrix.hpp>
#include <percept/math/Math.hpp>

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

  //static double macheps = std::numeric_limits<double>::epsilon();
  //static double sqrt_eps = std::sqrt(macheps);
  // static double cbrt_eps = std::pow(macheps, 1./3.);

  double ReferenceMeshSmoother2::local_line_search(stk::mesh::Entity node, stk::mesh::FieldBase *cg_g_field)
  {
    SmootherMetricFunction smf(m_eMesh, m_metric);
    double *cg_g = m_eMesh->field_data(cg_g_field, node);
    VERIFY_OP_ON(cg_g, !=, 0, "cg_g");
    double *coord = m_eMesh->field_data(m_eMesh->get_coordinates_field(), node);
    smf.set_node(node);

    double u[3];
    for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
      u[jc] = coord[jc];
    double alpha = 1.0;
    Double tau = 0.5;

    double sm0 = smf(u);
    Double min_alpha_factor = 1.e-12;

    bool converged = false;
    int liter = 0, niter = 1000;
    while (!converged && liter < niter)
      {
        for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
          {
            u[jc] = coord[jc] - alpha*cg_g[jc];
          }
        double sm1 = smf(u);
        if (sm1 < sm0)
          converged = true;

        if (!converged)
          alpha *= tau;
        if (alpha < min_alpha_factor*m_alpha_0)
          break;
        ++liter;
      }
    if (!converged)
      {
        //PRINT_1("not converged");
        alpha = 0.0;
      }
    return alpha;
  }

  // fills cg_s with locally solved Newton at nodes
  void ReferenceMeshSmoother2::get_step()
  {
    PerceptMesh *eMesh = m_eMesh;
    //stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
    //stk::mesh::FieldBase *coord_field_current   = coord_field;
    stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    stk::mesh::FieldBase *cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    //int spatialDim = eMesh->get_spatial_dim();

    m_scale = 1.e-10;

    // r=0
    eMesh->nodal_field_set_value(cg_g_field, 0.0);

    if (1)
      {
        // node loop
        SmootherMetricFunction smf(m_eMesh, m_metric);
        typedef FDGradient<SmootherMetricFunction> Gradient;
        typedef FDHessian<SmootherMetricFunction> Hessian;
        Gradient grad;
        Hessian hessian;
        DenseMatrix<3,1> G;
        DenseMatrix<3,3> H, HI;
        H.identity();
        HI.zero();
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
                    VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");

                    double *cg_edge_length = m_eMesh->field_data(cg_edge_length_field, node);
                    VERIFY_OP_ON(cg_edge_length, !=, 0, "cg_edge_length");
                    VERIFY_OP_ON(cg_edge_length_field, !=, 0, "cg_edge_length_field");

                    double edge_length_ave = cg_edge_length[0];


                    double *cg_g = m_eMesh->field_data(cg_g_field, node);
                    VERIFY_OP_ON(cg_g, !=, 0, "cg_g");
                    double *coord = m_eMesh->field_data(m_eMesh->get_coordinates_field(), node);
                    smf.set_node(node);

                    double u[3];
                    for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
                      u[jc] = coord[jc];

                    int nsiter = 10;
                    for (int siter = 0; siter < 1; ++siter)
                      {
                        grad(smf, u, G.data());
                        hessian(smf, u, H.data());
                        double Hdet = det(H);
                        if (Hdet == 0.0)
                          {
                            for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
                              {
                                cg_g[jc] = 0.0;
                              }
                            //break;
                          }
                        else
                          {
                            inverse(H, HI);
                            DenseMatrix<3,1> delta;
                            delta.zero();
                            //delta = (Hdet < 0 ? -1.0 : 1.0) * HI*gg;
                            delta = HI*G;
                            double dd = Frobenius(delta);
                            double fac = 10.;
                            if (0 && dd > edge_length_ave*fac)
                              {
                                delta *= edge_length_ave*fac/dd;
                              }
                            double ddnew = Frobenius(delta);
                            if (DEBUG_PRINT) std::cout << "here 6 dd= " << dd << " edge_length_ave= " << edge_length_ave << " new= " << ddnew << "\n";
                            for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
                              {
                                cg_g[jc] = delta(jc,0);
                              }

                            double alpha = local_line_search(node, cg_g_field);
                            //double alpha = 1.0;

                            //PRINT_1("node= " << node << " alpha= " << alpha);
                            for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
                              {
                                cg_g[jc] *= alpha;
                              }

                            double alp = 1.0;
                            double dm = 0.0;
                            for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
                              {
                                dm += cg_g[jc]*cg_g[jc];
                                //coord[jc] = coord[jc] - alp*cg_g[jc];
                                u[jc] = u[jc] - alp*cg_g[jc];
                              }
                            dm = std::sqrt(dm);
                            //std::cout << "siter= " << siter << " dm = " << dm << std::endl;
                            if (dm < edge_length_ave*1.e-5 || siter == nsiter-1)
                              {
                                for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
                                  {
                                    cg_g[jc] = -(u[jc] - coord[jc]);
                                  }
                                break;
                              }
                          }
                      }
                  }
              }
          }
      }

    //     CoordinatesFieldType *cg_g_field_v = static_cast<CoordinatesFieldType *>(cg_g_field);
    //     std::vector<stk::mesh::FieldBase*> sum_fields(1, cg_g_field_v);
    //     stk::mesh::parallel_sum(*m_eMesh->get_bulk_data(), sum_fields);

    // project deltas to surface
    if (m_eMesh->get_smooth_surfaces())
      {
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

                    double *cg_g = m_eMesh->field_data(cg_g_field, node);

                    if (fixed.second == MS_SURFACE)
                      {
                        project_delta_to_tangent_plane(node, cg_g);
                      }
                  }
              }
          }
      }

    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(cg_g_field);

      // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
      stk::mesh::copy_owned_to_shared(*eMesh->get_bulk_data(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*m_eMesh->get_bulk_data()->ghostings()[0], fields);
    }



  }

  double ReferenceMeshSmoother2::run_one_iteration()
  {
    PerceptMesh *eMesh = m_eMesh;

    stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    stk::mesh::FieldBase *cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
    stk::mesh::FieldBase *cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");
    stk::mesh::FieldBase *cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");

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
      }

    get_step();

    eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);
    eMesh->copy_field(cg_s_field, cg_r_field);
    //eMesh->copy_field(cg_s_field, cg_g_field);

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
