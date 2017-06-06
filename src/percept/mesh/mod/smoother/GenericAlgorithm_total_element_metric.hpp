// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GenericAlgorithm_total_element_metric_hpp
#define GenericAlgorithm_total_element_metric_hpp

#include <percept/MeshType.hpp>

namespace percept {

template<typename MeshType>
class SmootherMetricImpl;
class PerceptMesh;

///////////////////////////////////////////////////////////////////////

    template<typename MetricType>
    struct SGridGenericAlgorithm_total_element_metric
    {
      mutable MetricType m_metric;

      std::vector<typename StructuredGrid::MTElement> elements;
      Kokkos::View<int*, DataLayout , MemSpace > element_invalid_flags;

      long double& mtot;
      size_t& n_invalid;

      SGridGenericAlgorithm_total_element_metric(PerceptMesh *eMesh, long double& mtot_in, size_t& n_invalid_in);

      void run();

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned& index, long double& mtot_loc) const;
    };

    ///////////////////////////////////////////////////////////////////////

    template<typename MeshType>
    struct GenericAlgorithm_total_element_metric
    {
      typedef long double Double;
      using This = GenericAlgorithm_total_element_metric<MeshType>;

      SmootherMetricImpl<MeshType> *m_metric;

      PerceptMesh *m_eMesh;

      typename MeshType::MTSelector on_locally_owned_part;
      typename MeshType::MTSelector on_globally_shared_part;
      int spatialDim;
      typename MeshType::MTField *coord_field;
      typename MeshType::MTField *coord_field_current;
      typename MeshType::MTField *coord_field_lagged;

      typename MeshType::MTField *cg_s_field;

      std::vector<typename MeshType::MTElement> elements;
      Kokkos::View<int*, DataLayout , MemSpace > element_invalid_flags;
      std::vector<const typename MeshType::MTCellTopology *> topos;
      bool& valid;
      size_t *num_invalid;
      Double& mtot;
      size_t& n_invalid;

      GenericAlgorithm_total_element_metric(SmootherMetricImpl<MeshType> *metric, PerceptMesh *eMesh, bool& valid_in, size_t *num_invalid_in, Double& mtot_in, size_t& n_invalid_in);

      void run();

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned& index, Double& mtot_loc) const
      {
        const_cast<This *>(this)->operator()(index, mtot_loc);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned& index, Double& mtot_loc);
    };

} // namespace percept

#endif
