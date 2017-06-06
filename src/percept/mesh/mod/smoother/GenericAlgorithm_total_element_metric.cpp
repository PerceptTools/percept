#include <percept/mesh/mod/smoother/GenericAlgorithm_total_element_metric.hpp>

#include <percept/PerceptUtils.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/structured/BlockStructuredGrid.hpp>

#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>

namespace percept {

template<>
GenericAlgorithm_total_element_metric<STKMesh>::
GenericAlgorithm_total_element_metric(SmootherMetricImpl<STKMesh> *metric, PerceptMesh *eMesh, bool& valid_in, size_t *num_invalid_in, Double& mtot_in, size_t& n_invalid_in)
  : m_metric(metric), m_eMesh(eMesh), valid(valid_in), num_invalid(num_invalid_in), mtot(mtot_in), n_invalid(n_invalid_in)
{
  coord_field = m_eMesh->get_coordinates_field();
  coord_field_current   = coord_field;
  coord_field_lagged  = m_eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");

  cg_s_field    = m_eMesh->get_field(stk::topology::NODE_RANK, "cg_s");

  on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
  on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
  spatialDim = m_eMesh->get_spatial_dim();

  {
    elements.resize(0);
    topos.resize(0);

    // element loop
    const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->element_rank() );

    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (MeshSmootherImpl<STKMesh>::select_bucket(**k, m_eMesh) && on_locally_owned_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_elements_in_bucket = bucket.size();

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
#if defined(WITH_KOKKOS)
  Kokkos::Experimental::resize(element_invalid_flags, elements.size());
#else
  element_invalid_flags.resize(elements.size());
#endif
}

template<>
GenericAlgorithm_total_element_metric<StructuredGrid>::
GenericAlgorithm_total_element_metric(SmootherMetricImpl<StructuredGrid> *metric, PerceptMesh *eMesh, bool& valid_in, size_t *num_invalid_in, Double& mtot_in, size_t& n_invalid_in)
  : m_metric(metric), m_eMesh(eMesh), valid(valid_in), num_invalid(num_invalid_in), mtot(mtot_in), n_invalid(n_invalid_in)
{
  std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
  coord_field                              = bsg->m_fields["coordinates"].get();
  coord_field_current                      = coord_field;
  coord_field_lagged                       = bsg->m_fields["coordinates_lagged"].get();
  cg_s_field                               = bsg->m_fields["cg_s"].get();

  //on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
  //on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
  spatialDim = 3;

  bsg->get_elements(elements);
  topos.resize(elements.size(), static_cast<const typename StructuredGrid::MTCellTopology *>(0));

#if defined(WITH_KOKKOS)
  Kokkos::Experimental::resize(element_invalid_flags, elements.size());
#else
  element_invalid_flags.resize(elements.size());
#endif
}

template<typename MeshType>
void GenericAlgorithm_total_element_metric<MeshType>::
run()
{
#if defined WITH_KOKKOS && !(KOKKOS_HAVE_CUDA)
  stk::diag::Timer root_timer("GATM", rootTimerStructured());
  stk::diag::TimeBlock root_block(root_timer);

  {
    stk::diag::Timer reduce1("GATM pr 1", root_timer);
    stk::diag::TimeBlock reduce1_block(reduce1);

    // main loop computes local and global metrics
    Kokkos::parallel_reduce(elements.size(), *this, mtot);
  }

  {
  	stk::diag::Timer reduce1("GATM pr 2", root_timer);
  	stk::diag::TimeBlock reduce1_block(reduce1);

    // second loop computes n_invalid and valid data
    Kokkos::parallel_reduce(elements.size(), KOKKOS_LAMBDA (unsigned index, size_t &local_n_invalid) {
        local_n_invalid += element_invalid_flags(index);
      }, n_invalid);
  }
#else
  for (unsigned index=0; index<elements.size(); index++) {
    operator()(index, mtot);
  }

  for (unsigned index=0; index<elements.size(); index++) {
    n_invalid += element_invalid_flags[index];
  }
#endif
  valid = (n_invalid==0);
}


// ETI
template
void GenericAlgorithm_total_element_metric<STKMesh>::
run();

template
void GenericAlgorithm_total_element_metric<StructuredGrid>::
run();

///////////////////////////////////////////////////////////////////////

template<typename MetricType>
SGridGenericAlgorithm_total_element_metric<MetricType>::
SGridGenericAlgorithm_total_element_metric(PerceptMesh *eMesh, long double& mtot_in, size_t& n_invalid_in)
  : m_metric(eMesh), mtot(mtot_in), n_invalid(n_invalid_in)
{
  std::shared_ptr<BlockStructuredGrid> bsg = eMesh->get_block_structured_grid();
  bsg->get_elements(elements);
  Kokkos::Experimental::resize(element_invalid_flags, elements.size());
}

// ETI
template
SGridGenericAlgorithm_total_element_metric<SmootherMetricUntangleImpl<StructuredGrid> >::
SGridGenericAlgorithm_total_element_metric(PerceptMesh *eMesh, long double& mtot_in, size_t& n_invalid_in);

template<typename MetricType>
void SGridGenericAlgorithm_total_element_metric<MetricType>::
run()
{
  stk::diag::Timer root_timer("GATM", rootTimerStructured());
  stk::diag::TimeBlock root_block(root_timer);

  {
    stk::diag::Timer reduce1(std::string("GATM pr 1 ") + std::to_string((int)std::cbrt((double)elements.size())), root_timer);
    stk::diag::TimeBlock reduce1_block(reduce1);

    // main loop computes local and global metrics
    Kokkos::parallel_reduce(elements.size(), *this, mtot);
  }

  {
    // second loop computes n_invalid and valid data
    Kokkos::parallel_reduce(elements.size(), KOKKOS_LAMBDA (unsigned index, size_t &local_n_invalid) {
        local_n_invalid += element_invalid_flags(index);
      }, n_invalid);
  }
}

// ETI
template
void SGridGenericAlgorithm_total_element_metric<SmootherMetricUntangleImpl<StructuredGrid> >::run();

template<typename MetricType>
KOKKOS_INLINE_FUNCTION
void SGridGenericAlgorithm_total_element_metric<MetricType>::
operator()(const unsigned& index, long double& mtot_loc) const
{
  typename StructuredGrid::MTElement element = elements[index];

  bool local_valid = false;
  long double mm = m_metric.metric(element, local_valid);
  element_invalid_flags(index) = !local_valid;

  mtot_loc += mm;
}

// ETI
template
void SGridGenericAlgorithm_total_element_metric<SmootherMetricUntangleImpl<StructuredGrid> >::
operator()(const unsigned& index, long double& mtot_loc) const;

///////////////////////////////////////////////////////////////////////

template<>
#if defined(WITH_KOKKOS)
KOKKOS_INLINE_FUNCTION
#else
inline
#endif
void GenericAlgorithm_total_element_metric<STKMesh>::
operator()(const unsigned& index, Double& mtot_loc)
{
  typename STKMesh::MTElement element = elements[index];
  // FIXME
  m_metric->m_topology_data = topos[index];

  bool local_valid = false;
  Double mm = m_metric->metric(element, local_valid);
  element_invalid_flags(index) = !local_valid;

  mtot_loc += mm;
}

template<>
#if defined(WITH_KOKKOS)
KOKKOS_INLINE_FUNCTION
#else
inline
#endif
void GenericAlgorithm_total_element_metric<StructuredGrid>::

operator()(const unsigned& index, Double& mtot_loc)
{				//madbrew: the KKIF doesn't mind having this reference being passed?

  typename StructuredGrid::MTElement element = elements[index];


  bool local_valid = false;
  Double mm = m_metric->metric(element, local_valid);
  element_invalid_flags(index) = !local_valid;

  mtot_loc += mm;

}

} // namespace percept
