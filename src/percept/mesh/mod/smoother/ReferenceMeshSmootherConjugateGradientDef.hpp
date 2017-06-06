// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmootherConjugateGradientDef_hpp
#define ReferenceMeshSmootherConjugateGradientDef_hpp

#include <percept/PerceptUtils.hpp>

#include <percept/mesh/mod/smoother/GenericAlgorithm_total_element_metric.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_update_coordinates.hpp>

namespace percept {

template <typename MeshType>
ReferenceMeshSmootherConjugateGradientImpl<MeshType>::
ReferenceMeshSmootherConjugateGradientImpl(PerceptMesh *eMesh,
                                       typename MeshType::MTSelector *boundary_selector,
                                       typename MeshType::MTMeshGeometry *meshGeometry,
                                       int inner_iterations,
                                       double grad_norm,
                                       int parallel_iterations)
  :  Base(eMesh, boundary_selector, meshGeometry, inner_iterations, grad_norm, parallel_iterations), m_max_edge_length_factor(1.0)
{}

template<>
  ReferenceMeshSmootherConjugateGradientImpl< STKMesh >::Double
  ReferenceMeshSmootherConjugateGradientImpl< STKMesh >::
  nodal_edge_length_ave(stk::mesh::Entity node)
  {
    //int spatialDim = m_eMesh->get_spatial_dim();
    Double nm=0.0;

    double min=std::numeric_limits<double>::max();
    const MyPairIterRelation node_elems(*m_eMesh, node, m_eMesh->element_rank() );
    RMSCG_PRINT("tmp srk1 node_elems.size= " << node_elems.size());
    Double nele = 0.0;
    for (unsigned i_elem=0; i_elem < node_elems.size(); i_elem++)
      {
        stk::mesh::Entity element = node_elems[i_elem].entity();
        if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
          continue;
        double lmin=0,lmax=0;
        double elem_edge_len = m_eMesh->edge_length_ave(element, m_coord_field_original, &lmin, &lmax);
        RMSCG_PRINT("tmp srk1 elem_edge_len= " << elem_edge_len);
        if (lmin < min) min=lmin;
        nm += elem_edge_len;
        nele += 1.0;
      }
    nm /= nele;
    if (0)
      return min;
    return nm;
  }

  void find_connected_cells(PerceptMesh *eMesh, typename StructuredGrid::MTNode node, std::vector<StructuredCellIndex>& cells)
  {
    unsigned iblock = node[3];
    std::shared_ptr<StructuredBlock> sgrid = eMesh->get_block_structured_grid()->m_sblocks[iblock];
    // const unsigned A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];
    // const unsigned L0 = sgrid->m_loop_ordering[0], L1 = sgrid->m_loop_ordering[1], L2 = sgrid->m_loop_ordering[2];

    const int sizes[3] = {static_cast<int>(sgrid->m_sizes.node_size[0]),
                          static_cast<int>(sgrid->m_sizes.node_size[1]),
                          static_cast<int>(sgrid->m_sizes.node_size[2])};
    //unsigned Asizes[3] = {sgrid->m_sizes.node_size[A0], sgrid->m_sizes.node_size[A1], sgrid->m_sizes.node_size[A2]};

    //unsigned indx[3]{0,0,0};

    for (int i0 = node[0]-1; i0 <= static_cast<int>(node[0]); ++i0)
      {
        if (i0 < 0 || i0 > sizes[0]-2) continue;
        for (int i1 = node[1]-1; i1 <= static_cast<int>(node[1]); ++i1)
          {
            if (i1 < 0 || i1 > sizes[1]-2) continue;
            for (int i2 = node[2]-1; i2 <= static_cast<int>(node[2]); ++i2)
              {
                if (i2 < 0 || i2 > sizes[2]-2) continue;
                cells.emplace_back( StructuredCellIndex{{static_cast<unsigned>(i0),static_cast<unsigned>(i1),static_cast<unsigned>(i2),iblock}} );
              }
          }
      }
  }

  template<>
  ReferenceMeshSmootherConjugateGradientImpl< StructuredGrid >::Double
  ReferenceMeshSmootherConjugateGradientImpl< StructuredGrid >::
  nodal_edge_length_ave(typename StructuredGrid::MTNode node)
  {
    Double nm=0.0;

    double min=std::numeric_limits<double>::max();
    Double nele = 0.0;
    std::vector<StructuredCellIndex> node_elems;
    find_connected_cells(m_eMesh, node, node_elems);
    RMSCG_PRINT("tmp srk1 node_elems.size= " << node_elems.size());
    for (unsigned ii=0; ii < node_elems.size(); ++ii)
      {
        StructuredCellIndex element = node_elems[ii];
        double lmin=0,lmax=0;
        double elem_edge_len = m_eMesh->edge_length_ave(element, m_coord_field_original, &lmin, &lmax);
        RMSCG_PRINT("tmp srk1 elem_edge_len= " << elem_edge_len);
        if (lmin < min) min=lmin;
        nm += elem_edge_len;
        nele += 1.0;
      }
    nm /= nele;
    return nm;
  }


  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  debug_print(double alpha)
  {
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  snap_nodes
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename MeshType>
  struct GenericAlgorithm_snap_nodes
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using Double = typename RefMeshSmoother::Double;
    using This = GenericAlgorithm_snap_nodes<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    int spatialDim;
    typename MeshType::MTField *coord_field;
    typename MeshType::MTField *coord_field_current;

    std::vector<typename MeshType::MTNode> nodes;
    double& dmax;

    GenericAlgorithm_snap_nodes(RefMeshSmoother *rms, PerceptMesh *eMesh, double& dm);

    void run()
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index);


  };


  template<>
  GenericAlgorithm_snap_nodes<STKMesh>::
  GenericAlgorithm_snap_nodes(RefMeshSmoother *rms, PerceptMesh *eMesh, double& dm) : m_rms(rms), m_eMesh(eMesh), dmax(dm)
  {
    coord_field = m_eMesh->get_coordinates_field();
    coord_field_current   = coord_field;

    on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
    spatialDim = m_eMesh->get_spatial_dim();

    dmax=0.0;

    // node loop - build list...
    {
      nodes.resize(0);
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
                  nodes.push_back(node);
                }
            }
        }
    }
  }

  template<>
  GenericAlgorithm_snap_nodes<StructuredGrid>::
  GenericAlgorithm_snap_nodes(RefMeshSmoother *rms, PerceptMesh *eMesh, double& dm) : m_rms(rms), m_eMesh(eMesh), dmax(dm)
  {
    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    coord_field = bsg->m_fields["coordinates"].get();
    coord_field_current   = coord_field;

    //on_locally_owned_part = SGridSelector();
    //on_globally_shared_part =  SGridSelector();
    spatialDim = 3;

    dmax=0.0;

    bsg->get_nodes(nodes);
  }

  template<typename MeshType>
  void GenericAlgorithm_snap_nodes<MeshType>::
  operator()(int64_t& index)
  {
    typename MeshType::MTNode node = nodes[index];
    std::pair<bool,int> fixed = m_rms->get_fixed_flag(node);
    if (fixed.first)
      {
        return;
      }

    double coord_current[spatialDim];
    get_field<MeshType>(coord_current, spatialDim, m_eMesh, coord_field_current, node);

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
        m_rms->snap_to(node, &coord_current[0], keep_node_unchanged);
      }

    double dm = 0.0;
    for (int i=0; i < spatialDim; i++)
      {
        dm += (coord_old[i] - coord_project[i])*(coord_old[i] - coord_project[i]);
      }
    dm = std::sqrt(dm);
    dmax = std::max(dmax, dm);
  }

  template<>
  void ReferenceMeshSmootherConjugateGradientImpl< STKMesh >::
  snap_nodes()
  {
    static int anim_step = 0;
    bool save_anim = m_eMesh->getProperty("ReferenceMeshSmootherConjugateGradientImpl.snap_nodes.save_anim") == "true";
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

    double dmax=0.0;

    GenericAlgorithm_snap_nodes<STKMesh> ga(this, m_eMesh, dmax);
    ga.run();

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
        stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & dmax ) );
        if (m_eMesh->get_rank()==0) std::cout << "tmp srk snap_nodes dmax= " << dmax << std::endl;
      }
  }

  template<>
  void ReferenceMeshSmootherConjugateGradientImpl< StructuredGrid >::
  snap_nodes()
  {
    double dmax=0.0;

    GenericAlgorithm_snap_nodes<StructuredGrid> ga(this, m_eMesh, dmax);
    ga.run();
  }

  template<typename MeshType>
  typename ReferenceMeshSmootherConjugateGradientImpl< MeshType >::Double
  ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  total_metric( Double alpha, double multiplicative_edge_scaling, bool& valid, size_t* num_invalid)
  {
	  Double mtot = 0.0;
	  size_t n_invalid=0;

	  GenericAlgorithm_update_coordinates<MeshType> ga1(this,m_eMesh,alpha);
	  ga1.run();

	  valid = true;

	  GenericAlgorithm_total_element_metric<MeshType> ga2(this->m_metric, m_eMesh, valid, num_invalid, mtot, n_invalid);
	  ga2.run();

	  // reset coordinates
	  m_eMesh->copy_field(ga1.coord_field_current, ga1.coord_field_lagged);

	  stk::all_reduce( m_eMesh->parallel() , stk::ReduceSum<1>( & mtot ) );
	  stk::all_reduce( m_eMesh->parallel() , stk::ReduceMin<1>( & valid ) );

	  if (num_invalid)
	  {
		  *num_invalid = n_invalid;
		  stk::all_reduce( m_eMesh->parallel() , stk::ReduceSum<1>( num_invalid ) );
	  }

	  return mtot;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  update_node_positions
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename MeshType>
  struct GenericAlgorithm_update_node_positions
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using Double = typename RefMeshSmoother::Double;
    using This = GenericAlgorithm_update_node_positions<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    int spatialDim;
    typename MeshType::MTField *coord_field;
    typename MeshType::MTField *coord_field_current;
    typename MeshType::MTField *cg_s_field;
    typename MeshType::MTField *cg_edge_length_field;
    std::vector<typename MeshType::MTNode> nodes;
    Double alpha;
    GenericAlgorithm_update_node_positions(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alp);

    void run()
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index);


  };

  template<>
  GenericAlgorithm_update_node_positions<STKMesh>::
  GenericAlgorithm_update_node_positions(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alp) : m_rms(rms), m_eMesh(eMesh), alpha(alp)
  {
    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    spatialDim = eMesh->get_spatial_dim();
    coord_field = eMesh->get_coordinates_field();
    coord_field_current   = coord_field;
    cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    // node loop: update node positions
    nodes.resize(0);
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
                  nodes.push_back(node);
                }
            }
        }
    }

    m_rms->m_dmax = 0.0;
    m_rms->m_dmax_relative = 0.0;

  }

  template<>
  GenericAlgorithm_update_node_positions<StructuredGrid>::
  GenericAlgorithm_update_node_positions(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alp) : m_rms(rms), m_eMesh(eMesh), alpha(alp)
  {
    //on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    //on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    spatialDim = 3;

    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    coord_field                              = bsg->m_fields["coordinates"].get();
    coord_field_current                      = coord_field;
    cg_s_field                               = bsg->m_fields["cg_s"].get();
    cg_edge_length_field                     = bsg->m_fields["cg_edge_length"].get();

    // node loop: update node positions
    bsg->get_nodes(nodes);

    m_rms->m_dmax = 0.0;
    m_rms->m_dmax_relative = 0.0;
  }

  template<typename MeshType>
  void GenericAlgorithm_update_node_positions<MeshType>::
  operator()(int64_t& index)
  {
    PerceptMesh *eMesh = m_eMesh;
    typename MeshType::MTNode node = nodes[index];

    bool fixed = m_rms->get_fixed_flag(node).first;
    bool isGhostNode = MTisGhostNode<MeshType>(eMesh, node);
    if (fixed || isGhostNode)
      {
        return;
      }

    double coord_current[spatialDim];
    double cg_s[spatialDim];
    double cg_edge_length[1];

    get_field<MeshType>(coord_current, spatialDim, m_eMesh, coord_field_current, node);
    get_field<MeshType>(cg_s, spatialDim, m_eMesh, cg_s_field, node);
    get_field<MeshType>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);

    for (int i=0; i < spatialDim; i++)
      {
        Double dt = alpha*cg_s[i];
        m_rms->m_dmax = std::max(std::fabs(dt), m_rms->m_dmax);
        m_rms->m_dmax_relative = std::max(std::fabs(dt)/cg_edge_length[0], m_rms->m_dmax_relative);
        coord_current[i] += dt;
      }
    set_field<MeshType>(coord_current, spatialDim, m_eMesh, coord_field_current, node);
  }

  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  update_node_positions( Double alpha)
  {
    GenericAlgorithm_update_node_positions<MeshType> gaun(this, m_eMesh, alpha);
    gaun.run();

    stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & m_dmax ) );
    stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & m_dmax_relative ) );

    {
      std::vector< const typename MeshType::MTField *> fields;
      fields.push_back(gaun.coord_field);
      MTcommFields<MeshType>(fields, m_eMesh);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  line_search
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename MeshType>
  struct GenericAlgorithm_line_search
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using Double = typename RefMeshSmoother::Double;
    using This = GenericAlgorithm_line_search<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTField *cg_g_field;
    typename MeshType::MTField *cg_s_field;
    typename MeshType::MTField *cg_r_field;

    bool restarted;
    double mfac_mult;
    Double alpha;
    bool extra_print;
    GenericAlgorithm_line_search(RefMeshSmoother *rms, PerceptMesh *eMesh, double mf_mult);
    void run();
  };

  template<>
  GenericAlgorithm_line_search<STKMesh>::
  GenericAlgorithm_line_search(RefMeshSmoother *rms, PerceptMesh *eMesh, double mf_mult) : m_rms(rms), m_eMesh(eMesh), mfac_mult(mf_mult)
  {
    restarted = false;
    extra_print = false;
    if (m_eMesh->getProperty("ReferenceMeshSmootherConjugateGradientImpl::line_search.extra_print") == "true")
      extra_print = true;
    cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
  }

  template<>
  GenericAlgorithm_line_search<StructuredGrid>::
  GenericAlgorithm_line_search(RefMeshSmoother *rms, PerceptMesh *eMesh, double mf_mult) : m_rms(rms), m_eMesh(eMesh), mfac_mult(mf_mult)
  {
    restarted = false;
    extra_print = false;
    if (m_eMesh->getProperty("ReferenceMeshSmootherConjugateGradientImpl::line_search.extra_print") == "true")
      extra_print = true;

    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    cg_g_field                               = bsg->m_fields["cg_g"].get();
    cg_s_field                               = bsg->m_fields["cg_s"].get();
    cg_r_field                               = bsg->m_fields["cg_r"].get();
  }

  template<typename MeshType>
  void GenericAlgorithm_line_search<MeshType>::
  run()
  {
    PerceptMesh *eMesh = m_eMesh;
    m_rms->m_alpha_0 = m_rms->get_alpha_0();
    const Double alpha_fac = 10.0;
    alpha = alpha_fac*m_rms->m_alpha_0;

    bool total_valid = false;
    size_t n_invalid = 0;
    size_t* n_invalid_p = &n_invalid;
    Double metric_0 = m_rms->total_metric( 0.0, 1.0, total_valid, n_invalid_p);
    Double metric = 0.0;
    Double tau = 0.5;
    Double c0 = 1.e-1;
    Double min_alpha_factor = 1.e-12;

    //RMSCG_PRINT_1("metric_0= " << metric_0 << " m_stage= " << m_stage << " m_iter= " << m_iter);

    Double sDotGrad = eMesh->nodal_field_dot(cg_s_field, cg_g_field);
    if (sDotGrad >= 0.0)
      {
        Double sDotGradOld = sDotGrad;
        eMesh->copy_field(cg_s_field, cg_r_field);
        m_rms->m_alpha_0 = m_rms->get_alpha_0();
        alpha = alpha_fac*m_rms->m_alpha_0;
        sDotGrad = eMesh->nodal_field_dot(cg_s_field, cg_g_field);
        RMSCG_PRINT_1("sDotGradOld= " << sDotGradOld << " sDotGrad= " << sDotGrad << " m_stage= " << m_rms->m_stage << " m_iter= " << m_rms->m_iter);
        restarted = true;
      }
    VERIFY_OP_ON(sDotGrad, <, 0.0, "bad sDotGrad");

    Double armijo_offset_factor = c0*sDotGrad;
    bool converged = false;
    total_valid = false;
    int liter = 0, niter = 1000;
    while (!converged && liter < niter)
      {
        metric = m_rms->total_metric(alpha, 1.0, total_valid, n_invalid_p);

        Double mfac = alpha*armijo_offset_factor * mfac_mult;
        converged = (metric < metric_0 + mfac);
        if (m_rms->m_untangled) converged = converged && total_valid;
        if (extra_print) RMSCG_PRINT_1(  "tmp srk alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac)
                                   << " sDotGrad= " << sDotGrad << " mfac= " << mfac << " n_invalid= " << n_invalid << " m_untangled = " << m_rms->m_untangled
                                   << " total_valid= " << total_valid << " converged= " << converged);
        if (!converged)
          alpha *= tau;
        if (alpha < min_alpha_factor*m_rms->m_alpha_0)
          break;
        ++liter;
      }

    if (!converged)
      {
        restarted = true;
        eMesh->copy_field(cg_s_field, cg_r_field);
        m_rms->m_alpha_0 = m_rms->get_alpha_0();
        alpha = alpha_fac*m_rms->m_alpha_0;
        sDotGrad = eMesh->nodal_field_dot(cg_s_field, cg_g_field);
        RMSCG_PRINT_1("not converged, trying restart, sDotGrad new= " << sDotGrad << " m_stage= " << m_rms->m_stage << " m_iter= " << m_rms->m_iter);
        VERIFY_OP_ON(sDotGrad, <, 0.0, "bad sDotGrad 2nd time");
      }

    liter = 0;
    while (!converged && liter < niter)
      {
        metric = m_rms->total_metric(alpha, 1.0, total_valid);

        Double mfac = alpha*armijo_offset_factor;
        converged = (metric < metric_0 + mfac);
        if (m_rms->m_untangled)
          {
            converged = converged && total_valid;
            if (total_valid && m_rms->m_dmax_relative < m_rms->gradNorm)
              converged = true;
          }

        RMSCG_PRINT_1(  "alpha 2nd time= " << alpha << " alpha_0= " << m_rms->m_alpha_0 << " sDotGrad= " << sDotGrad << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac)
                  << " m_untangled = " << m_rms->m_untangled << " m_stage= " << m_rms->m_stage
                  << " total_valid= " << total_valid );
        if (!converged)
          alpha *= tau;
        if (alpha < min_alpha_factor*m_rms->m_alpha_0)
          break;
        ++liter;
      }

    if (!converged)
      {
        RMSCG_PRINT_1( "WARNING: can't reduce metric 2nd time = " << metric << " metric_0 + armijo_offset " << metric_0+alpha*armijo_offset_factor << " sDotGrad = " << sDotGrad << " alpha_0= " << m_rms->m_alpha_0 << " alpha= " << alpha);
        throw std::runtime_error("can't reduce metric");
      }
    else
      {
        Double a1 = alpha/2.;
        Double a2 = alpha;
        Double f0 = metric_0, f1 = m_rms->total_metric( a1, 1.0, total_valid), f2 = m_rms->total_metric( a2, 1.0, total_valid);
        Double den = 2.*(a2*(-f0 + f1) + a1*(f0 - f2));
        Double num = a2*a2*(f1-f0)+a1*a1*(f0-f2);
        if (std::fabs(den) > 1.e-10)
          {
            Double alpha_quadratic = num/den;
            if (alpha_quadratic > 1.e-10 && alpha_quadratic < 2*alpha)
              {
                Double fm=m_rms->total_metric( alpha_quadratic, 1.0, total_valid);
                if (fm < f2 && (m_rms->m_stage==0 || total_valid))
                  {
                    alpha = alpha_quadratic;
                    if (extra_print) RMSCG_PRINT_1( "\ntmp srk alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 << " f0= " << f0 << " f2= " << f2 << " fq= " << fm << "\n");
                  }
                if (fm < f2 && (m_rms->m_stage!=0 && !total_valid))
                  {
                    RMSCG_PRINT_1( "WARNING: !total_valid alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 );
                  }
              }
          }
      }

  }

  template<typename MeshType>
  typename ReferenceMeshSmootherConjugateGradientImpl< MeshType >::Double
  ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  line_search(bool& restarted, double mfac_mult)
  {
    GenericAlgorithm_line_search<MeshType> gal(this, m_eMesh, mfac_mult);
    gal.run();
    restarted = gal.restarted;
    return gal.alpha;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_surface_normals
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(STK_PERCEPT_HAS_GEOMETRY)
  template<typename MeshType>
  struct GenericAlgorithm_get_surface_normals
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using Double = typename RefMeshSmoother::Double;
    using This = GenericAlgorithm_get_surface_normals<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTField *cg_normal_field;
    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    std::vector<double> norm;
    std::vector<typename MeshType::MTNode> nodes;

    GenericAlgorithm_get_surface_normals(RefMeshSmoother *rms, PerceptMesh *eMesh);

    void run() const
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      typename MeshType::MTNode node = nodes[index];
      double cg_normal[m_eMesh->get_spatial_dim()];
      get_field<MeshType>(cg_normal, m_eMesh->get_spatial_dim(), m_eMesh, cg_normal_field, node);

      std::pair<bool,int> fixed = m_rms->get_fixed_flag(node);
      if (!fixed.first && (fixed.second == MS_SURFACE || fixed.second == MS_ON_BOUNDARY))
        {
          m_rms->m_meshGeometry->normal_at(m_eMesh, node, norm);
          for (int ii=0; ii < m_eMesh->get_spatial_dim(); ++ii)
            {
              cg_normal[ii] = norm[ii];
            }
          set_field<MeshType>(cg_normal, m_eMesh->get_spatial_dim(), m_eMesh, cg_normal_field, node);
        }
    }

  };

  template<>
  GenericAlgorithm_get_surface_normals<STKMesh>::
  GenericAlgorithm_get_surface_normals(RefMeshSmoother *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    norm.resize(3,0.0);
    cg_normal_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_normal");
    VERIFY_OP_ON(cg_normal_field, !=, 0, "must have cg_normal_field");

    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    nodes.resize(0);
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
                nodes.push_back(node);
              }
          }
      }
  }

  template<>
  GenericAlgorithm_get_surface_normals<StructuredGrid>::
  GenericAlgorithm_get_surface_normals(RefMeshSmoother *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    norm.resize(3,0.0);
    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    cg_normal_field                          = bsg->m_fields["cg_normal"].get();
    VERIFY_OP_ON(cg_normal_field, !=, 0, "must have cg_normal_field");

    //on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    //on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    bsg->get_nodes(nodes);
  }
#endif

  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_surface_normals(PerceptMesh * eMesh)
  {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    GenericAlgorithm_get_surface_normals< MeshType > ga_gsn(this, eMesh);
    ga_gsn.run();
#else
    VERIFY_MSG("you must configure Percept with geometry (STK_PERCEPT_HAS_GEOMETRY=1) to use get_surface_normals");
#endif
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_edge_lengths
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<typename MeshType>
  struct GenericAlgorithm_get_edge_lengths
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using Double = typename RefMeshSmoother::Double;
    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;
    typename MeshType::MTField *cg_edge_length_field;
    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    std::vector<typename MeshType::MTNode> nodes;

    using This = GenericAlgorithm_get_edge_lengths<MeshType>;

    GenericAlgorithm_get_edge_lengths(RefMeshSmoother *rms, PerceptMesh *eMesh);


    void run() const
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      typename MeshType::MTNode node = nodes[index];
      double cg_edge_length[1];
      get_field<MeshType>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);

      //if (on_locally_owned_part(node) || on_globally_shared_part(node))
      {
        Double edge_length_ave = m_rms->nodal_edge_length_ave(node);
        cg_edge_length[0] = edge_length_ave;
        set_field<MeshType>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);
      }

    }
  };

  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_edge_lengths(PerceptMesh * eMesh)
  {
    GenericAlgorithm_get_edge_lengths<MeshType> gae(this, eMesh);
    gae.run();
  }

  template<>
  GenericAlgorithm_get_edge_lengths<STKMesh>::
  GenericAlgorithm_get_edge_lengths(RefMeshSmoother *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    nodes.resize(0);

    cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");
    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

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
                nodes.push_back(node);
              }
          }
      }
  }

  template<>
  GenericAlgorithm_get_edge_lengths<StructuredGrid>::
  GenericAlgorithm_get_edge_lengths(RefMeshSmoother *rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    cg_edge_length_field                     = bsg->m_fields["cg_edge_length"].get();

    //on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    //on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    bsg->get_nodes(nodes);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_alpha_0
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  template<typename MeshType>
  struct GenericAlgorithm_get_alpha_0
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTField *cg_s_field;
    typename MeshType::MTField *cg_edge_length_field;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    int spatialDim;

    using Double = typename RefMeshSmoother::Double;

    Double alpha;
    bool alpha_set;

    std::vector<typename MeshType::MTNode> nodes;
    using This = GenericAlgorithm_get_alpha_0<MeshType>;

    GenericAlgorithm_get_alpha_0(RefMeshSmoother *rms, PerceptMesh *eMesh);

    void run()
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }

      stk::all_reduce( m_eMesh->parallel() , stk::ReduceMax<1>( & alpha_set ) );
      if (!alpha_set)
        alpha = 1.0;

      stk::all_reduce( m_eMesh->parallel() , stk::ReduceMin<1>( & alpha ) );
      VERIFY_OP_ON(alpha, > , 0.0, "bad alpha");

    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      typename MeshType::MTNode node = nodes[index];

      //VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");

      double cg_edge_length[1];
      get_field<MeshType>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);
      Double edge_length_ave = cg_edge_length[0];

      bool isGhostNode = MTisGhostNode<MeshType>(m_eMesh, node);
      VERIFY_OP_ON(isGhostNode, ==, false, "hmmmm");
      bool fixed = m_rms->get_fixed_flag(node).first;
      RMSCG_PRINT("tmp srk1 index= " << index << " edge_length_ave= " << edge_length_ave << " fixed= " << fixed);
      if (fixed || isGhostNode)
        return;

      double cg_s[spatialDim];
      get_field<MeshType>(cg_s, spatialDim, m_eMesh, cg_s_field, node);
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

  };

  /// gets a global scale factor so that local gradient*scale is approximately the size of the local mesh edges
  /// also uses the reference mesh to compute a local scaled gradient norm for convergence checks
  template<typename MeshType>
  typename ReferenceMeshSmootherConjugateGradientImpl< MeshType >::Double
  ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_alpha_0()
  {
    GenericAlgorithm_get_alpha_0<MeshType> ga0(this, m_eMesh);
    ga0.run();
    return ga0.alpha;
  }

  template<>
  GenericAlgorithm_get_alpha_0<STKMesh>::
  GenericAlgorithm_get_alpha_0(RefMeshSmoother * rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    spatialDim = eMesh->get_spatial_dim();

    alpha = std::numeric_limits<double>::max();
    alpha_set = false;

    nodes.resize(0);

    {
      const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          // FIXME
          if (MeshSmootherImpl<STKMesh>::select_bucket(**k, m_eMesh) && (on_locally_owned_part(**k) || on_globally_shared_part(**k)))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");
                  nodes.push_back(node);
                }
            }
        }
    }
    RMSCG_PRINT("tmp srk1 nodes.size= " << nodes.size());
  }

  template<>
  GenericAlgorithm_get_alpha_0<StructuredGrid>::
  GenericAlgorithm_get_alpha_0(RefMeshSmoother * rms, PerceptMesh *eMesh) : m_rms(rms), m_eMesh(eMesh)
  {
    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    cg_s_field                               = bsg->m_fields["cg_s"].get();
    cg_edge_length_field                     = bsg->m_fields["cg_edge_length"].get();

    //on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    //on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    spatialDim = 3;

    alpha = std::numeric_limits<double>::max();
    alpha_set = false;
    bsg->get_nodes(nodes);
    RMSCG_PRINT("tmp srk1 nodes.size= " << nodes.size());

  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  check_convergence
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<typename MeshType>
  bool ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  check_convergence()
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  run_one_iter
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  template<typename MeshType>
  struct GenericAlgorithm_run_one_iter
  {
    PerceptMesh *m_eMesh;

    typename MeshType::MTField *cg_g_field;
    typename MeshType::MTField *cg_r_field;
    typename MeshType::MTField *cg_d_field;
    typename MeshType::MTField *cg_s_field;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;

    GenericAlgorithm_run_one_iter(PerceptMesh *eMesh);

  };

  template<>
  GenericAlgorithm_run_one_iter<STKMesh>::
  GenericAlgorithm_run_one_iter(PerceptMesh *eMesh) : m_eMesh(eMesh)
  {
    cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
    cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");
    cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");

    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

  };

  template<>
  GenericAlgorithm_run_one_iter<StructuredGrid>::
  GenericAlgorithm_run_one_iter(PerceptMesh *eMesh) : m_eMesh(eMesh)
  {
    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    cg_g_field                               = bsg->m_fields["cg_g"].get();
    cg_r_field                               = bsg->m_fields["cg_r"].get();
    cg_d_field                               = bsg->m_fields["cg_d"].get();
    cg_s_field                               = bsg->m_fields["cg_s"].get();

    //on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    //on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

  };

  template<typename MeshType>
  double ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  run_one_iteration()
  {
	stk::diag::Timer cumulative_timer(m_eMesh->getProperty("in_filename"), rootTimerStructured());
	stk::diag::Timer runOneIter("RMSCGI::run_one_iteration()",cumulative_timer);
	stk::diag::TimeBlock my_time_block(runOneIter);

    PerceptMesh *m_eMesh = Base::m_eMesh;

    GenericAlgorithm_run_one_iter<MeshType> gab(m_eMesh);

    typename MeshType::MTField *cg_g_field    = gab.cg_g_field;
    typename MeshType::MTField *cg_r_field    = gab.cg_r_field;
    typename MeshType::MTField *cg_d_field    = gab.cg_d_field;
    typename MeshType::MTField *cg_s_field    = gab.cg_s_field;

    // typename MeshType::MTSelector on_locally_owned_part =  gab.on_locally_owned_part;
    // typename MeshType::MTSelector on_globally_shared_part =  gab.on_globally_shared_part;

    bool total_valid=false;

    if (Base::m_iter == 0)
      {
        get_edge_lengths(m_eMesh);
        m_eMesh->nodal_field_set_value(cg_g_field, 0.0);
        m_eMesh->nodal_field_set_value(cg_r_field, 0.0);
        m_eMesh->nodal_field_set_value(cg_d_field, 0.0);
        m_eMesh->nodal_field_set_value(cg_s_field, 0.0);
      }

    /// f'(x)
    get_gradient();

    /// r = -g
    m_eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);

    Base::m_dold = m_eMesh->nodal_field_dot(cg_d_field, cg_d_field);
    Base::m_dmid = m_eMesh->nodal_field_dot(cg_r_field, cg_d_field);
    Base::m_dnew = m_eMesh->nodal_field_dot(cg_r_field, cg_r_field);

    if (Base::m_iter == 0)
      {
        Base::m_d0 = Base::m_dnew;
      }

    Double metric_check = total_metric( 0.0, 1.0, total_valid); //it seems like every function like this is reinstantiating its corresponding structs with each iter
    Base::m_total_metric = metric_check;

    if (check_convergence() || metric_check == 0.0)
      {
        RMSCG_PRINT_1( "INFO: already converged m_dnew= " << Base::m_dnew << " gradNorm= " << Base::gradNorm << " metric_check= " << metric_check );
        //update_node_positions
        return total_metric(0.0,1.0, total_valid);
      }

    Double cg_beta = 0.0;
    if (Base::m_dold == 0.0)
      cg_beta = 0.0;
    else if (Base::m_iter > 0)
      cg_beta = (Base::m_dnew - Base::m_dmid) / Base::m_dold;

    RMSCG_PRINT("tmp srk beta = " << cg_beta);

    size_t N = Base::m_num_nodes;
    if (Base::m_iter % N == 0 || cg_beta <= 0.0)
      {
        /// s = r
        m_eMesh->copy_field(cg_s_field, cg_r_field);
      }
    else
      {
        /// s = r + beta * s
        m_eMesh->nodal_field_axpby(1.0, cg_r_field, cg_beta, cg_s_field);
      }

    m_eMesh->copy_field(cg_d_field, cg_r_field);

    bool restarted = false;
    Double alpha = line_search(restarted);
    Double snorm = m_eMesh->nodal_field_dot(cg_s_field, cg_s_field);
    Base::m_grad_norm_scaled = Base::m_alpha_0*std::sqrt(snorm)/Double(Base::m_num_nodes);
    RMSCG_PRINT("tmp srk1 alpha= " << alpha << " snorm= " << snorm << " m_grad_norm_scaled= " << Base::m_grad_norm_scaled << " m_num_nodes= " << Base::m_num_nodes << " m_alpha_0= " << Base::m_alpha_0);

    /// x = x + alpha*d
    Base::m_alpha = alpha;
    update_node_positions(alpha);

    // check if re-snapped geometry is acceptable
    if (m_eMesh->get_smooth_surfaces())
      {
        snap_nodes();
        if (Base::m_stage != 0)
          {
            bool total_valid_0=true;
            total_metric( 0.0, 1.0, total_valid_0);
            VERIFY_OP_ON(total_valid_0, ==, true, "bad mesh after snap_node_positions...");
          }
      }

    Double tm = total_metric(0.0,1.0, total_valid);

    return tm;
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////
  ////////  get_gradient
  ////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  template<class MeshType>
  struct GenericAlgorithmBase_get_gradient {

    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;
    typename MeshType::MTField *coord_field;
    typename MeshType::MTField *coord_field_current;
    typename MeshType::MTField *cg_g_field;
    typename MeshType::MTField *cg_r_field;
    //typename MeshType::MTField *cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");
    typename MeshType::MTField *cg_edge_length_field;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;

    int spatialDim;

    GenericAlgorithmBase_get_gradient(RefMeshSmoother *rms, PerceptMesh *eMesh);

  };

  template<>
  GenericAlgorithmBase_get_gradient<STKMesh>::
  GenericAlgorithmBase_get_gradient(RefMeshSmoother *rms, PerceptMesh *eMesh) :  m_rms(rms), m_eMesh(eMesh), spatialDim(eMesh->get_spatial_dim())
  {
    rms->m_scale = 1.e-10;

    coord_field = eMesh->get_coordinates_field();
    coord_field_current   = coord_field;
    cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
    cg_r_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_r");
    cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

  }

  template<>
  GenericAlgorithmBase_get_gradient<StructuredGrid>::
  GenericAlgorithmBase_get_gradient(RefMeshSmoother *rms, PerceptMesh *eMesh) :  m_rms(rms), m_eMesh(eMesh), spatialDim(eMesh->get_spatial_dim())
  {
    rms->m_scale = 1.e-10;

    //on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    //on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
    coord_field                              = bsg->m_fields["coordinates"].get();
    coord_field_current                      = coord_field;
    cg_g_field                               = bsg->m_fields["cg_g"].get();
    cg_r_field                               = bsg->m_fields["cg_r"].get();
    cg_edge_length_field                     = bsg->m_fields["cg_edge_length"].get();

  }

  template<class MeshType>
  struct GenericAlgorithm_get_gradient_1 : public GenericAlgorithmBase_get_gradient<MeshType> {

    using Base = GenericAlgorithmBase_get_gradient<MeshType> ;
    using This = GenericAlgorithm_get_gradient_1<MeshType>;

    using Base::m_eMesh;
    using Base::coord_field;
    using Base::coord_field_current;
    using Base::cg_g_field;
    using Base::cg_r_field;
    //typename MeshType::MTField *cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");
    using Base::cg_edge_length_field;

    using Base::on_locally_owned_part;
    using Base::on_globally_shared_part;
    using Base::spatialDim;

    using RefMeshSmoother = ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    GenericAlgorithm_get_gradient_1(RefMeshSmoother *rms, PerceptMesh *eMesh) : Base(rms, eMesh)
    {
    }

    std::vector<typename MeshType::MTElement> elements;
    std::vector<const typename MeshType::MTCellTopology *> topos;

    void MTSetTopo(int64_t index) const;

    void run() const
    {
      for (int64_t index = 0; index < (int64_t)elements.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {

      const double macheps = std::numeric_limits<double>::epsilon();
      const double sqrt_eps = std::sqrt(macheps);

      typename MeshType::MTElement element = elements[index];
      MTSetTopo(index);

      unsigned num_node = get_num_nodes<MeshType>(m_eMesh, element);
      std::vector<typename MeshType::MTNode> nodes_buffer;
      const typename MeshType::MTNode *elem_nodes = get_nodes<MeshType>(m_eMesh, element, &nodes_buffer);

      using Double = typename RefMeshSmoother::Double;
      Double edge_length_ave = 0;

      const bool use_analytic_grad = true;
      double analytic_grad[8][4];
      const bool test_analytic_grad = false;
      VERIFY_OP_ON(Base::m_rms->m_metric->has_gradient(), ==, true, "bad has_gradient");
      if ((test_analytic_grad || use_analytic_grad) && Base::m_rms->m_stage >= 0 && Base::m_rms->m_metric->has_gradient())
        {
          bool gmvalid = true;
          Double gm = Base::m_rms->m_metric->grad_metric(element, gmvalid, analytic_grad);
          (void)gm;
          if ((gmvalid || Base::m_rms->m_stage == 0) && !test_analytic_grad)
            {
              for (unsigned inode=0; inode < num_node; inode++)
                {
                  typename MeshType::MTNode node = elem_nodes[ inode ];

                  bool isGhostNode = MTisGhostNode<MeshType>(m_eMesh, node);
                  bool node_locally_owned = MTnode_locally_owned<MeshType>(m_eMesh, node);
                  bool fixed = Base::m_rms->get_fixed_flag(node).first;
                  // if (fixed)
                  //   std::cout << "inode= " << inode << " is fixed" << std::endl;

                  if (fixed || isGhostNode)
                    continue;

                  VERIFY_OP_ON(Base::spatialDim, ==, spatialDim, "bad spatialDim");
                  VERIFY_OP_ON(Base::spatialDim, >=, 2, "bad spatialDim");
                  double cg_g[Base::spatialDim];
                  get_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);
                  for (int jdim=0; jdim < spatialDim; jdim++)
                    {
                      if (node_locally_owned)
                        cg_g[jdim] += analytic_grad[inode][jdim];
                      else
                        cg_g[jdim] = 0.0;
                    }
                  // if (1 || inode == 0)
                  //   std::cout << "inode= " << inode << " cg_g= " << Math::print_3d(cg_g) << std::endl;
                  set_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);
                }
            }
          if (!test_analytic_grad)
            return;
        }

      if (1)
        {
          // finite-different grad
          for (unsigned inode=0; inode < num_node; inode++)
            {
              typename MeshType::MTNode node = elem_nodes[ inode ];

              bool isGhostNode = MTisGhostNode<MeshType>(m_eMesh, node);
              bool node_locally_owned = MTnode_locally_owned<MeshType>(m_eMesh, node);
              bool fixed = Base::m_rms->get_fixed_flag(node).first;
              if (fixed || isGhostNode)
                continue;

              double cg_edge_length[1];
              get_field<MeshType>(cg_edge_length, 1, Base::m_rms->m_eMesh, cg_edge_length_field, node);

              edge_length_ave = cg_edge_length[0];

              Base::m_rms->m_metric->set_node(node);
              double coord_current[spatialDim];
              get_field<MeshType>(coord_current, spatialDim, Base::m_rms->m_eMesh, coord_field_current, node);
              double cg_g[spatialDim];
              get_field<MeshType>(cg_g, spatialDim, Base::m_rms->m_eMesh, cg_g_field, node);

              //Double eps1 = cbrt_eps*edge_length_ave;
              Double eps1 = sqrt_eps*edge_length_ave;

              double gsav[3]={0,0,0};
              for (int idim=0; idim < spatialDim; idim++)
                {
                  Double cc = coord_current[idim];
                  coord_current[idim] += eps1;
                  set_field<MeshType>(coord_current, spatialDim,  idim, Base::m_rms->m_eMesh, coord_field_current, node);
                  bool pvalid=false, mvalid=false;
                  Double mp = Base::m_rms->m_metric->metric(element, pvalid);
                  const bool second_order = true;
                  if (second_order)
                    coord_current[idim] -= 2.0*eps1;
                  else
                    coord_current[idim] = cc;
                  set_field<MeshType>(coord_current, spatialDim,  idim, Base::m_rms->m_eMesh, coord_field_current, node);
                  Double mm = Base::m_rms->m_metric->metric(element, mvalid);
                  coord_current[idim] = cc;
                  set_field<MeshType>(coord_current, spatialDim,  idim, Base::m_rms->m_eMesh, coord_field_current, node);
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
                  set_field<MeshType>(cg_g, spatialDim, idim, Base::m_rms->m_eMesh, cg_g_field, node);
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
                    }

                }
            }
        }
    }
  };


  template<>
  void GenericAlgorithm_get_gradient_1<STKMesh>::
  MTSetTopo(int64_t index) const
  {
    Base::m_rms->m_metric->m_topology_data = topos[index];
  }

  template<>
  void GenericAlgorithm_get_gradient_1<StructuredGrid>::
  MTSetTopo(int64_t index) const
  {
    //Base::m_rms->m_metric->m_topology_data = topos[index];
  }

  template<>
  GenericAlgorithm_get_gradient_1<STKMesh>::
  GenericAlgorithm_get_gradient_1(ReferenceMeshSmootherConjugateGradientImpl<STKMesh> *rms, PerceptMesh *eMesh)
    : GenericAlgorithmBase_get_gradient<STKMesh>(rms, eMesh)
  {
    eMesh->nodal_field_set_value(cg_g_field, 0.0);

    // get elements
    if (1)
      {
        const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (MeshSmootherImpl<STKMesh>::select_bucket(**k, m_eMesh))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();
                Base::m_rms->m_metric->m_topology_data = m_eMesh->get_cell_topology(bucket);

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];

                    if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
                      continue;
                    elements.push_back(element);
                    topos.push_back(m_eMesh->get_cell_topology(bucket));
                  }
              }
          }
      }
  }

  template<>
  GenericAlgorithm_get_gradient_1<StructuredGrid>::
  GenericAlgorithm_get_gradient_1(ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> *rms, PerceptMesh *eMesh) : GenericAlgorithmBase_get_gradient<StructuredGrid>(rms, eMesh)
  {
    eMesh->nodal_field_set_value(cg_g_field, 0.0);

    // get elements
    //FIXME if (MeshSmootherImpl<STKMesh>::select_bucket(**k, m_eMesh))
    std::shared_ptr<BlockStructuredGrid> bsg = eMesh->get_block_structured_grid();
    bsg->get_elements(elements);
    topos.resize(elements.size(), static_cast<const typename StructuredGrid::MTCellTopology *>(0));
  }

  template<typename MeshType>
  struct GenericAlgorithm_get_gradient_2 : public GenericAlgorithmBase_get_gradient<MeshType>
  {
    std::vector<typename MeshType::MTNode> nodes;

    using Base = GenericAlgorithmBase_get_gradient<MeshType>;
    using This = GenericAlgorithm_get_gradient_2<MeshType>;

    using Base::m_eMesh;
    using Base::spatialDim;
    using Base::cg_g_field;
    using Base::cg_r_field;

    GenericAlgorithm_get_gradient_2(ReferenceMeshSmootherConjugateGradientImpl<MeshType> *rms, PerceptMesh *eMesh);

    void run() const
    {
      for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<This *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      // project deltas to surface
      if (m_eMesh->get_smooth_surfaces())
        {
          typename MeshType::MTNode node = nodes[index];

          std::pair<bool,int> fixed = Base::m_rms->get_fixed_flag(node);
          if (!fixed.first)
            {
              std::vector<double> cg_g(spatialDim);
              if (fixed.second == MS_SURFACE)
                {
                  double cg_g[spatialDim];
                  get_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);
                  Base::m_rms->project_delta_to_tangent_plane(node, cg_g);
                  set_field<MeshType>(cg_g, spatialDim, m_eMesh, cg_g_field, node);
                }
            }
        }
    }
  };

  template<>
  GenericAlgorithm_get_gradient_2<STKMesh>::
  GenericAlgorithm_get_gradient_2(ReferenceMeshSmootherConjugateGradientImpl<STKMesh> *rms, PerceptMesh *eMesh) :  GenericAlgorithmBase_get_gradient<STKMesh>(rms, eMesh)
  {
    nodes.clear();

    const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                nodes.push_back(bucket[i_node]);
              }
          }
      }
  }

  template<>
  GenericAlgorithm_get_gradient_2<StructuredGrid>::
  GenericAlgorithm_get_gradient_2(ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> *rms, PerceptMesh *eMesh) :  GenericAlgorithmBase_get_gradient<StructuredGrid>(rms, eMesh)
  {
    eMesh->get_block_structured_grid()->get_nodes(nodes);
  }

  /// fills cg_g_field with f'(x)
  template<typename MeshType>
  void ReferenceMeshSmootherConjugateGradientImpl< MeshType >::
  get_gradient()
  {
    GenericAlgorithm_get_gradient_1<MeshType> ga_1(this, m_eMesh);
    ga_1.run();

    //Double gnorm = m_eMesh->nodal_field_dot(ga_1.cg_g_field, ga_1.cg_g_field);
    //std::cout << "gnorm= " << gnorm << std::endl;

    std::vector<typename MeshType::MTField *> fields_0(1, ga_1.cg_g_field);
    MTsum_fields<MeshType>(fields_0, m_eMesh);

    m_eMesh->copy_field(ga_1.cg_r_field, ga_1.cg_g_field);

    //Double rnorm = m_eMesh->nodal_field_dot(ga_1.cg_r_field, ga_1.cg_r_field);
    //std::cout << "rnorm= " << rnorm << std::endl;

    GenericAlgorithm_get_gradient_2<MeshType> ga_2(this, m_eMesh);
    ga_2.run();

    {
      std::vector<  const  typename MeshType::MTField *> fields;
      fields.push_back(ga_1.cg_g_field);
      fields.push_back(ga_1.cg_r_field);
      MTcommFields<MeshType>(fields, m_eMesh);
    }

  }

#if defined (WITH_KOKKOS) //&& defined(KOKKOS_HAVE_OPENMP)

  template<>
  templeKokkos<int>::templeKokkos(unsigned into)
  {
  	size = into;
  	Kokkos::View<double*, DataLayout , MemSpace > temp("thaName",size);
  	templeDubs = temp;
  }

  template<>
  void templeKokkos<int>::operator()(int64_t index) const
  		{
  			templeDubs(index) *= 2.8;
  			templeDubs(index) *= 2.7;
  			templeDubs(index) *= 2.8;
  			templeDubs(index) += 2.8;
  		}

  template<>
  void templeKokkos<int>::run()
  {
	stk::diag::Timer cumulativePF("templeKokkosPF",rootTimerStructured());
	stk::diag::TimeBlock pf(cumulativePF);
  	Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,(int64_t)size),*this);
  }

#endif

}//percept

#endif
