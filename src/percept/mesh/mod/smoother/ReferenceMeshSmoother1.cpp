// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)


#include <percept/mesh/mod/smoother/ReferenceMeshSmoother1.hpp>
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

  static double macheps = std::numeric_limits<double>::epsilon();
  static double sqrt_eps = std::sqrt(macheps);
  static double cbrt_eps = std::pow(macheps, 1./3.);

  //const bool do_tot_test = false;
  bool do_print_elem_val = false;

  static bool s_do_snap = false;

  ReferenceMeshSmoother1::Double ReferenceMeshSmoother1::nodal_edge_length_ave(stk::mesh::Entity node)
  {
    //int spatialDim = m_eMesh->get_spatial_dim();
    Double nm=0.0;

    double min=std::numeric_limits<double>::max();
    const MyPairIterRelation node_elems(*m_eMesh, node, m_eMesh->element_rank() );
    Double nele = 0.0;
    for (unsigned i_elem=0; i_elem < node_elems.size(); i_elem++)
      {
        stk::mesh::Entity element = node_elems[i_elem].entity();
        if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
          continue;
        double lmin=0,lmax=0;
        double elem_edge_len = m_eMesh->edge_length_ave(element, m_coord_field_original, &lmin, &lmax);
        if (lmin < min) min=lmin;
        nm += elem_edge_len;
        nele += 1.0;
      }
    nm /= nele;
    if (0)
      return min;
    return nm;
  }

  /// gets a global scale factor so that local gradient*scale is approximately the size of the local mesh edges
  /// also uses the reference mesh to compute a local scaled gradient norm for convergence checks
  ReferenceMeshSmoother1::Double ReferenceMeshSmoother1::get_alpha_0()
  {
    PerceptMesh *eMesh = m_eMesh;
    stk::mesh::FieldBase *cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    stk::mesh::FieldBase *cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    int spatialDim = eMesh->get_spatial_dim();

    Double alpha = std::numeric_limits<double>::max();
    bool alpha_set = false;

    {
      const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // FIXME
          if (MeshSmoother::select_bucket(**k, m_eMesh) && (on_locally_owned_part(**k) || on_globally_shared_part(**k)))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");

                  double *cg_edge_length = m_eMesh->field_data(cg_edge_length_field, node);
                  Double edge_length_ave = cg_edge_length[0];

                  bool isGhostNode = !(on_locally_owned_part(m_eMesh->bucket(node)) || on_globally_shared_part(m_eMesh->bucket(node)));
                  VERIFY_OP_ON(isGhostNode, ==, false, "hmmmm");
                  bool fixed = this->get_fixed_flag(node).first;
                  if (fixed || isGhostNode)
                    continue;

                  double *cg_s = m_eMesh->field_data(cg_s_field, node);
                  Double sn = 0.0;
                  for (int idim=0; idim < spatialDim; idim++)
                    {
                      sn += cg_s[idim]*cg_s[idim];
                    }
                  sn = std::sqrt(sn);
                  if (sn > 0.0)
                    {
                      Double alpha_new = edge_length_ave / sn;
                      if (!alpha_set)
                        {
                          alpha_set = true;
                          alpha = alpha_new;
                        }
                      else if (alpha_new < alpha)
                        {
                          alpha = alpha_new;
                        }
                    }
                }
            }
        }
    }

    stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & alpha_set ) );
    if (!alpha_set)
      alpha = 1.0;

    stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMin<1>( & alpha ) );
    VERIFY_OP_ON(alpha, > , 0.0, "bad alpha");

    return alpha;
  }

  bool ReferenceMeshSmoother1::check_convergence()
  {
    Double grad_check = gradNorm;
    bool retval=false;
    //std::cout << "P[" << m_eMesh->get_rank() << "] tmp srk a1.0 m_num_invalid= " << m_num_invalid << " m_dnew= " << m_dnew << " m_dmax = " << m_dmax << " m_grad_norm_scaled= " << m_grad_norm_scaled << std::endl;
    int type = -1;
    if (m_stage == 0 && (m_dnew == 0.0 || m_total_metric == 0.0))
      {
        //std::cout << "P[" << m_eMesh->get_rank() << "] tmp srk a1.1 untangle m_dnew= " << m_dnew << " m_total_metric = " << m_total_metric << std::endl;
        retval = true; // for untangle
        type=1;
      }
    else if (m_stage == 0 && m_num_invalid == 0 && (m_dmax_relative < grad_check))
      {
        //std::cout << "P[" << m_eMesh->get_rank() << "] tmp srk a1.2 untangle m_num_invalid= " << m_num_invalid << " m_dmax = " << m_dmax << " m_grad_norm_scaled= " << m_grad_norm_scaled << std::endl;
        retval = true;
        type=2;
      }
    //grad_check = 1.e-8;
    else if (m_num_invalid == 0 && (m_iter > 10 && (m_dmax_relative < grad_check && (m_dnew < grad_check*grad_check*m_d0 || m_grad_norm_scaled < grad_check))))
      {
        if (0 && m_eMesh->get_rank() == 0)
          {
            std::cout << "P[" << m_eMesh->get_rank() << "] tmp srk a1.3 untangle m_dnew(scaled nnode) check= "
                      << " m_grad_norm_scaled check= " << (m_grad_norm_scaled < grad_check)
                      << " m_dmax check= " << (m_dmax < grad_check)
                      << " m_dnew check= " << (m_dnew < grad_check*grad_check*m_d0)
                      << " m_total_metric = " << m_total_metric << std::endl;
          }
        type=3;
        retval = true;
      }
    if (0)
      {
        if (m_eMesh->get_rank() == 0)
          std::cout << "P[" << m_eMesh->get_rank() << "] tmp srk a1 retval= " << retval << " type= " << type << " m_num_invalid= " << m_num_invalid << " m_total_metric= " << m_total_metric
                    << " m_dnew= " << m_dnew << " m_dmax = " << m_dmax
                    << " m_grad_norm_scaled= " << m_grad_norm_scaled << " m_iter= " << m_iter
                    << " m_d0= " << m_d0 << " m_num_nodes= " << m_num_nodes << std::endl;
        //std::cout << "P[" << m_eMesh->get_rank() << "] tmp srk a1 retval= " << retval << " m_dnew= " << m_dnew << " m_total_metric = " << m_total_metric << std::endl;
      }
    return retval;
  }

  ReferenceMeshSmoother1::Double ReferenceMeshSmoother1::metric(stk::mesh::Entity element, bool& valid)
  {
    return m_metric->metric(element,valid);
  }

  /// fills cg_g_field with f'(x)
  void ReferenceMeshSmoother1::get_gradient()
  {
    PerceptMesh *eMesh = m_eMesh;
    stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
    stk::mesh::FieldBase *coord_field_current   = coord_field;
    stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    stk::mesh::FieldBase *cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
    //stk::mesh::FieldBase *cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");
    stk::mesh::FieldBase *cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    int spatialDim = eMesh->get_spatial_dim();

    m_scale = 1.e-10;

    // g=0
    eMesh->nodal_field_set_value(cg_g_field, 0.0);

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
                m_metric->m_topology_data = m_eMesh->get_cell_topology(bucket);

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];
                    if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
                      continue;
                    const MyPairIterRelation elem_nodes(*m_eMesh, element, stk::topology::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    Double edge_length_ave = 0; // m_eMesh->edge_length_ave(element, m_coord_field_original);

                    const bool use_analytic_grad = true;
                    double analytic_grad[8][4];
                    const bool test_analytic_grad = false;
                    if ((test_analytic_grad || use_analytic_grad) && m_stage >= 0 && m_metric->has_gradient())
                      {
                        bool gmvalid = true;
                        Double gm = m_metric->grad_metric(element, gmvalid, analytic_grad);
                        (void)gm;
                        if ((gmvalid || m_stage == 0) && !test_analytic_grad)
                          {
                            for (unsigned inode=0; inode < num_node; inode++)
                              {
                                stk::mesh::Entity node = elem_nodes[ inode ].entity();

                                bool isGhostNode = !(on_locally_owned_part(m_eMesh->bucket(node)) || on_globally_shared_part(m_eMesh->bucket(node)));
                                bool node_locally_owned = (eMesh->get_rank() == eMesh->owner_rank(node));
                                bool fixed = this->get_fixed_flag(node).first;
                                if (fixed || isGhostNode)
                                  continue;

                                double *cg_g = m_eMesh->field_data(cg_g_field, node);
                                for (int jdim=0; jdim < spatialDim; jdim++)
                                  {
                                    if (node_locally_owned)
                                      cg_g[jdim] += analytic_grad[inode][jdim];
                                    else
                                      cg_g[jdim] = 0.0;
                                  }
                              }
                          }
                        if (!test_analytic_grad)
                          continue;
                      }

                      {
                        // finite-different grad
                        for (unsigned inode=0; inode < num_node; inode++)
                          {
                            stk::mesh::Entity node = elem_nodes[ inode ].entity();

                            bool isGhostNode = !(on_locally_owned_part(m_eMesh->bucket(node)) || on_globally_shared_part(m_eMesh->bucket(node)));
                            bool node_locally_owned = (eMesh->get_rank() == eMesh->owner_rank(node));
                            bool fixed = this->get_fixed_flag(node).first;
                            if (fixed || isGhostNode)
                              continue;

                            // FIXME
                            double *cg_edge_length = m_eMesh->field_data(cg_edge_length_field, node);

                            edge_length_ave = cg_edge_length[0];

                            m_metric->set_node(node);
                            double *coord_current = m_eMesh->field_data(coord_field_current, node);
                            double *cg_g = m_eMesh->field_data(cg_g_field, node);

                            //Double eps1 = cbrt_eps*edge_length_ave;
                            Double eps1 = sqrt_eps*edge_length_ave;

                            double gsav[3]={0,0,0};
                            for (int idim=0; idim < spatialDim; idim++)
                              {
                                Double cc = coord_current[idim];
                                coord_current[idim] += eps1;
                                bool pvalid=false, mvalid=false;
                                Double mp = metric(element, pvalid);
                                const bool second_order = true;
                                if (second_order)
                                  coord_current[idim] -= 2.0*eps1;
                                else
                                  coord_current[idim] = cc;
                                Double mm = metric(element, mvalid);
                                coord_current[idim] = cc;
                                Double dd = 0.0;
                                if (second_order)
                                  {
                                    dd = (mp - mm)/(2*eps1);
                                  }
                                else
                                  {
                                    dd = (mp - mm)/(eps1);
                                  }
                                gsav[idim] = dd;

                                //if (std::fabs(mp) > 1.e-10) std::cout << "tmp srk mp = " << mp << " mm= " << mm << " dd= " << dd << std::endl;

                                if (node_locally_owned)
                                  {
                                    cg_g[idim] += dd;
                                  }
                                else
                                  {
                                    cg_g[idim] = 0.0;
                                  }
                              }
                          if (test_analytic_grad)
                            {
                              double fd = std::max(std::fabs( gsav[0]),std::fabs( gsav[1]));
                              if (spatialDim==3) fd = std::max(fd, std::fabs( gsav[2]));
                              double ag = std::max(std::fabs( analytic_grad[inode][0]),std::fabs( analytic_grad[inode][1] ));
                              if (spatialDim==3) ag = std::max(ag, std::fabs( analytic_grad[inode][2] ));
                              double diff = std::fabs(ag-fd);
                              //if (ag > 1.e-3 && diff > 1.e-3*ag)
                              //if (diff > 1.e-3/edge_length_ave)
                              double comp = (ag+fd)/2.0;
                              if (comp > 1.e-6 && diff > 1.e-8 && diff > 1.e-3*comp)
                                {
                                  std::cout << "analytic_grad= " << analytic_grad[inode][0] << " " << analytic_grad[inode][1] << " "
                                            << " fd_grad= " << gsav[0] << " " << gsav[1]
                                            << " diff = " << analytic_grad[inode][0]-gsav[0] << " " << analytic_grad[inode][1] -gsav[1] 
                                            << " ag= " << ag << " fd= " << fd
                                            << std::endl;
                                  //exit(1);
                                }

                            }
                          } // inode...
                      } // use_analytic_grad
                  } // i_element
              } // on_locally_owned_part...
          } // buckets
      }

    CoordinatesFieldType *cg_g_field_v = static_cast<CoordinatesFieldType *>(cg_g_field);
    std::vector<stk::mesh::FieldBase*> sum_fields(1, cg_g_field_v);
    stk::mesh::parallel_sum(*m_eMesh->get_bulk_data(), sum_fields);

    m_eMesh->copy_field(cg_r_field, cg_g_field);

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
      fields.push_back(cg_r_field);

      // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
      stk::mesh::copy_owned_to_shared(*eMesh->get_bulk_data(), fields);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*m_eMesh->get_bulk_data()->ghostings()[0], fields);
    }

  }

  void ReferenceMeshSmoother1::get_edge_lengths(PerceptMesh * eMesh)
  {
    stk::mesh::FieldBase *cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");
    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                stk::mesh::Entity node = bucket[i_node];
                double *cg_edge_length = m_eMesh->field_data(cg_edge_length_field, node);

                //if (on_locally_owned_part(node) || on_globally_shared_part(node))
                {
                  Double edge_length_ave = nodal_edge_length_ave(node);
                  cg_edge_length[0] = edge_length_ave;
                }
              }
          }
      }
  }

  void ReferenceMeshSmoother1::get_surface_normals(PerceptMesh * eMesh)
  {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    stk::mesh::FieldBase *cg_normal_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_normal");
    if (!cg_normal_field) return;

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    std::vector<double> norm(3,0.0);
    const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                stk::mesh::Entity node = bucket[i_node];
                double *cg_normal = m_eMesh->field_data(cg_normal_field, node);

                std::pair<bool,int> fixed = this->get_fixed_flag(node);
                if (!fixed.first && (fixed.second == MS_SURFACE || fixed.second == MS_ON_BOUNDARY))
                  {
                    m_meshGeometry->normal_at(eMesh, node, norm);
                    for (int ii=0; ii < eMesh->get_spatial_dim(); ++ii)
                      {
                        cg_normal[ii] = norm[ii];
                      }
                  }
              }
          }
      }
#else
    VERIFY_MSG("you must configure Percept with geometry (STK_PERCEPT_HAS_GEOMETRY=1) to use get_surface_normals");
#endif
  }


  double ReferenceMeshSmoother1::run_one_iteration()
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

    /// f'(x)
    get_gradient();

    /// r = -g
    eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);

    m_dold = eMesh->nodal_field_dot(cg_d_field, cg_d_field);
    m_dmid = eMesh->nodal_field_dot(cg_r_field, cg_d_field);
    m_dnew = eMesh->nodal_field_dot(cg_r_field, cg_r_field);
    if (m_iter == 0)
      {
        m_d0 = m_dnew;
      }

    Double metric_check = total_metric( 0.0, 1.0, total_valid);
    m_total_metric = metric_check;

    if (check_convergence() || metric_check == 0.0)
      {
        PRINT_1( "INFO: already converged m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " metric_check= " << metric_check );
        //update_node_positions
        return total_metric(0.0,1.0, total_valid);
      }

    Double cg_beta = 0.0;
    if (m_dold == 0.0)
      cg_beta = 0.0;
    else if (m_iter > 0)
      cg_beta = (m_dnew - m_dmid) / m_dold;

    PRINT("tmp srk beta = " << cg_beta);

    // FIXME
    size_t N = m_num_nodes;
    if (m_iter % N == 0 || cg_beta <= 0.0)
      {
        /// s = r
        eMesh->copy_field(cg_s_field, cg_r_field);
      }
    else
      {
        /// s = r + beta * s
        eMesh->nodal_field_axpby(1.0, cg_r_field, cg_beta, cg_s_field);
      }

    eMesh->copy_field(cg_d_field, cg_r_field);

    bool restarted = false;
    Double alpha = line_search(restarted);
    Double snorm = eMesh->nodal_field_dot(cg_s_field, cg_s_field);
    m_grad_norm_scaled = m_alpha_0*std::sqrt(snorm)/Double(m_num_nodes);

    /// x = x + alpha*d
    m_alpha = alpha;
    update_node_positions(alpha);
    //PRINT_1( "tmp srk iter= "<< m_iter << " dmax= " << m_dmax << " alpha= " << alpha);

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

    Double tm = total_metric(0.0,1.0, total_valid);

    return tm;
  }

  ReferenceMeshSmoother1::Double ReferenceMeshSmoother1::line_search(bool& restarted, double mfac_mult)
  {
    PerceptMesh *eMesh = m_eMesh;
    restarted = false;
    bool extra_print = false;
    if (m_eMesh->getProperty("ReferenceMeshSmoother1::line_search.extra_print") == "true")
      extra_print = true;
    stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    stk::mesh::FieldBase *cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    stk::mesh::FieldBase *cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");

    m_alpha_0 = get_alpha_0();
    const Double alpha_fac = 10.0;
    Double alpha = alpha_fac*m_alpha_0;

    bool total_valid = false;
    size_t n_invalid = 0;
    size_t* n_invalid_p = &n_invalid;
    Double metric_0 = total_metric( 0.0, 1.0, total_valid, n_invalid_p);
    Double metric = 0.0;
    Double tau = 0.5;
    Double c0 = 1.e-1;
    Double min_alpha_factor = 1.e-12;

    //PRINT_1("metric_0= " << metric_0 << " m_stage= " << m_stage << " m_iter= " << m_iter);

    Double sDotGrad = eMesh->nodal_field_dot(cg_s_field, cg_g_field);
    if (sDotGrad >= 0.0)
      {
        Double sDotGradOld = sDotGrad;
        eMesh->copy_field(cg_s_field, cg_r_field);
        m_alpha_0 = get_alpha_0();
        alpha = alpha_fac*m_alpha_0;
        sDotGrad = eMesh->nodal_field_dot(cg_s_field, cg_g_field);
        PRINT_1("sDotGradOld= " << sDotGradOld << " sDotGrad= " << sDotGrad << " m_stage= " << m_stage << " m_iter= " << m_iter);
        restarted = true;
      }
    VERIFY_OP_ON(sDotGrad, <, 0.0, "bad sDotGrad");

    Double armijo_offset_factor = c0*sDotGrad;
    bool converged = false;
    total_valid = false;
    int liter = 0, niter = 1000;
    while (!converged && liter < niter)
      {
        metric = total_metric(alpha, 1.0, total_valid, n_invalid_p);

        Double mfac = alpha*armijo_offset_factor * mfac_mult;
        converged = (metric < metric_0 + mfac);
        if (m_untangled) converged = converged && total_valid;
        if (extra_print) PRINT_1(  "tmp srk alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac)
                                   << " sDotGrad= " << sDotGrad << " mfac= " << mfac << " n_invalid= " << n_invalid << " m_untangled = " << m_untangled
                                   << " total_valid= " << total_valid << " converged= " << converged);
        if (!converged)
          alpha *= tau;
        if (alpha < min_alpha_factor*m_alpha_0)
          break;
        ++liter;
      }

    if (!converged)
      {
        restarted = true;
        eMesh->copy_field(cg_s_field, cg_r_field);
        m_alpha_0 = get_alpha_0();
        alpha = alpha_fac*m_alpha_0;
        sDotGrad = eMesh->nodal_field_dot(cg_s_field, cg_g_field);
        PRINT_1("not converged, trying restart, sDotGrad new= " << sDotGrad << " m_stage= " << m_stage << " m_iter= " << m_iter);
        VERIFY_OP_ON(sDotGrad, <, 0.0, "bad sDotGrad 2nd time");
      }

    liter = 0;
    while (!converged && liter < niter)
      {
        metric = total_metric(alpha, 1.0, total_valid);

        Double mfac = alpha*armijo_offset_factor;
        converged = (metric < metric_0 + mfac);
        if (m_untangled) converged = converged && total_valid;
        PRINT_1(  "alpha 2nd time= " << alpha << " alpha_0= " << m_alpha_0 << " sDotGrad= " << sDotGrad << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac)
                  << " m_untangled = " << m_untangled << " m_stage= " << m_stage
                  << " total_valid= " << total_valid );
        if (!converged)
          alpha *= tau;
        if (alpha < min_alpha_factor*m_alpha_0)
          break;
        ++liter;
      }

    if (!converged)
      {
        PRINT_1( "WARNING: can't reduce metric 2nd time = " << metric << " metric_0 + armijo_offset " << metric_0+alpha*armijo_offset_factor << " sDotGrad = " << sDotGrad << " alpha_0= " << m_alpha_0 << " alpha= " << alpha);
        throw std::runtime_error("can't reduce metric");
      }
    else
      {
        Double a1 = alpha/2.;
        Double a2 = alpha;
        Double f0 = metric_0, f1 = total_metric( a1, 1.0, total_valid), f2 = total_metric( a2, 1.0, total_valid);
        Double den = 2.*(a2*(-f0 + f1) + a1*(f0 - f2));
        Double num = a2*a2*(f1-f0)+a1*a1*(f0-f2);
        if (std::fabs(den) > 1.e-10)
          {
            Double alpha_quadratic = num/den;
            if (alpha_quadratic > 1.e-10 && alpha_quadratic < 2*alpha)
              {
                Double fm=total_metric( alpha_quadratic, 1.0, total_valid);
                if (fm < f2 && (m_stage==0 || total_valid))
                  {
                    alpha = alpha_quadratic;
                    if (extra_print) PRINT_1( "\ntmp srk alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 << " f0= " << f0 << " f2= " << f2 << " fq= " << fm << "\n");
                  }
                if (fm < f2 && (m_stage!=0 && !total_valid))
                  {
                    PRINT_1( "WARNING: !total_valid alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 );
                  }
              }
          }
      }
    return alpha;
  }

  void ReferenceMeshSmoother1::debug_print(double alpha)
  {
  }

  void ReferenceMeshSmoother1::update_node_positions( Double alpha)
  {
    PerceptMesh *eMesh = m_eMesh;
    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    int spatialDim = eMesh->get_spatial_dim();
    stk::mesh::FieldBase *cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    stk::mesh::FieldBase *cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    m_dmax = 0.0;
    m_dmax_relative = 0.0;

    // node loop: update node positions
    {
      const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // update local and globally shared
          if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  bool fixed = this->get_fixed_flag(node).first;
                  bool isGhostNode = !(on_locally_owned_part(m_eMesh->bucket(node)) || on_globally_shared_part(m_eMesh->bucket(node)));
                  if (fixed || isGhostNode)
                    {
                      continue;
                    }

                  double *coord_current = m_eMesh->field_data(m_coord_field_current, node);
                  double *cg_s = m_eMesh->field_data(cg_s_field, node);
                  double *cg_edge_length = m_eMesh->field_data(cg_edge_length_field, node);

                  for (int i=0; i < spatialDim; i++)
                    {
                      Double dt = alpha*cg_s[i];
                      m_dmax = std::max(std::fabs(dt), m_dmax);
                      m_dmax_relative = std::max(std::fabs(dt)/cg_edge_length[0], m_dmax_relative);
                      coord_current[i] += dt;
                    }
                }
            }
        }
    }

    stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & m_dmax ) );
    stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & m_dmax_relative ) );

    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(m_eMesh->get_coordinates_field());
      //fields.push_back(cg_g_field);

      // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
      stk::mesh::copy_owned_to_shared(*eMesh->get_bulk_data(), fields);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*m_eMesh->get_bulk_data()->ghostings()[0], fields);
    }


  }

  ReferenceMeshSmoother1::Double ReferenceMeshSmoother1::total_metric( Double alpha, double multiplicative_edge_scaling, bool& valid, size_t* num_invalid)
  {
    stk::mesh::FieldBase *coord_field = m_eMesh->get_coordinates_field();
    stk::mesh::FieldBase *coord_field_current   = coord_field;
    stk::mesh::FieldBase *coord_field_lagged  = m_eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");

    stk::mesh::FieldBase *cg_s_field    = m_eMesh->get_field(stk::topology::NODE_RANK, "cg_s");

    stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
    int spatialDim = m_eMesh->get_spatial_dim();

    Double mtot = 0.0;
    size_t n_invalid=0;

    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(cg_s_field);
      fields.push_back(coord_field);
      stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
    }

    // cache coordinates
    m_eMesh->copy_field(coord_field_lagged, coord_field_current);

    // node loop - update positions
    {
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // update local and globally shared
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

                  double *coord_current = m_eMesh->field_data(coord_field_current, node);
                  double *cg_s = m_eMesh->field_data(cg_s_field, node);

                  // shouldn't be necessary
                  if (0 && fixed.second == MS_SURFACE)
                    {
                      project_delta_to_tangent_plane(node, cg_s);
                    }

                  double coord_project[3] = {0,0,0};
                  for (int i=0; i < spatialDim; i++)
                    {
                      Double dt = alpha * cg_s[i];
                      coord_current[i] += dt;
                      coord_project[i] = coord_current[i];
                    }
                  if (fixed.second == MS_SURFACE && (m_stage == 0 || s_do_snap))
                    {
                      if (0) snap_to(node, coord_project, false);
                    }
                }
            }
        }
    }

    valid = true;

    {
      // element loop
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->element_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (MeshSmoother::select_bucket(**k, m_eMesh) && on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              m_metric->m_topology_data = m_eMesh->get_cell_topology(bucket);

              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                {
                  stk::mesh::Entity element = bucket[i_element];
                  if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
                    continue;
                  bool local_valid=true;
                  Double mm = metric(element, local_valid);
                  if (!local_valid) n_invalid++;

                  valid = valid && local_valid;
                  if (do_print_elem_val) PRINT( "element= " << m_eMesh->identifier(element) << " metric= " << mm );
                  if (do_print_elem_val && m_eMesh->identifier(element) == 13) { std::cout << m_eMesh->identifier(element) << " iter= " << m_iter << " element= " << m_eMesh->identifier(element) << " metric= " << mm << std::endl;}
                  mtot += mm;
                }
            }
        }
    }

    // reset coordinates
    m_eMesh->copy_field(coord_field_current, coord_field_lagged);

    stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceSum<1>( & mtot ) );
    stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMin<1>( & valid ) );

    if (num_invalid)
      {
        *num_invalid = n_invalid;
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceSum<1>( num_invalid ) );
      }

    return mtot;

  }

  void ReferenceMeshSmoother1::snap_nodes()
  {
    static int anim_step = 0;
    bool save_anim = m_eMesh->getProperty("ReferenceMeshSmoother1.snap_nodes.save_anim") == "true";
    if (save_anim)
      {
        std::ostringstream fileid_ss;
        fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);
        std::string oname = "snap.e";
        if (anim_step > 0) oname += "-s" + fileid_ss.str();
        m_eMesh->save_as(oname);
        ++anim_step;
      }

    std::string prevSetting;
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    if (m_meshGeometry && m_meshGeometry->geomKernel)
      {
        prevSetting = m_meshGeometry->geomKernel->get_property("GKGP:use_unprojected_coordinates");
        m_meshGeometry->geomKernel->set_property("GKGP:use_unprojected_coordinates", "false");
      }
#endif
    stk::mesh::FieldBase *coord_field = m_eMesh->get_coordinates_field();
    stk::mesh::FieldBase *coord_field_current   = coord_field;

    stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
    int spatialDim = m_eMesh->get_spatial_dim();

    double dmax=0.0;

    // node loop - snap...
    {
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // update local and globally shared
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

                  double *coord_current = m_eMesh->field_data(coord_field_current, node);

                  double coord_project[3] = {0,0,0};
                  double coord_old[3] = {0,0,0};
                  for (int i=0; i < spatialDim; i++)
                    {
                      coord_project[i] = coord_current[i];
                      coord_old[i] = coord_current[i];
                    }

                  if (fixed.second == MS_SURFACE)
                    {
                      bool keep_node_unchanged = false;
                      //snap_to(node, coord_project, keep_node_unchanged);
                      snap_to(node, coord_current, keep_node_unchanged);
                    }

                  double dm = 0.0;
                  for (int i=0; i < spatialDim; i++)
                    {
                      dm += (coord_old[i] - coord_project[i])*(coord_old[i] - coord_project[i]);
                    }
                  dm = std::sqrt(dm);
                  dmax = std::max(dmax, dm);
                }
            }
        }
    }

    if (save_anim)
      {
        std::ostringstream fileid_ss;
        fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);
        std::string oname = "snap.e";
        if (anim_step > 0) oname += "-s" + fileid_ss.str();
        m_eMesh->save_as(oname);
        ++anim_step;
      }

    // FIXME - add coordinate comm
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    m_meshGeometry->geomKernel->set_property("GKGP:use_unprojected_coordinates", prevSetting);
#endif
    if (0)
      {
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & dmax ) );
        if (m_eMesh->get_rank()==0) std::cout << "tmp srk snap_nodes dmax= " << dmax << std::endl;
      }
  }


}


#endif
