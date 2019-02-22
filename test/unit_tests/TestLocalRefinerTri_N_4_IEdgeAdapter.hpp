// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TestLocalRefinerTri_N_4_IEdgeAdapter_hpp
#define adapt_TestLocalRefinerTri_N_4_IEdgeAdapter_hpp

#include <adapt/IEdgeAdapter.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation as a use case for IEdgeAdapter
     */
    class TestLocalRefinerTri_N_4_IEdgeAdapter : public IEdgeAdapter
    {
    public:
      TestLocalRefinerTri_N_4_IEdgeAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0)
      : IEdgeAdapter(eMesh, bp, proc_rank_field), location(0.0), threshold(0.0), max_length(0.0) {}
      void set_location(double loc) { location = loc; }
      void set_threshold(double thresh) { threshold = thresh; }
      void set_max_length(double max_len) { max_length = max_len; }
    protected:

      /// Client supplies these methods - given an element, which edge, and the nodes on the edge, return instruction on what to do to the edge,
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int markEdge(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
                           double *coord0, double *coord1, std::vector<int>* existing_edge_marks) ;
      double location;
      double threshold;
      double max_length;
    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_4_IEdgeAdapter in UnitTestLocalRefiner.cpp)

    inline
    int TestLocalRefinerTri_N_4_IEdgeAdapter::
    markEdge(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
             double *coord0, double *coord1, std::vector<int>* existing_edge_marks)
    {
      // Edge length
      const stk::mesh::MetaData & stk_meta = *m_eMesh.get_fem_meta_data();
      const UInt spatial_dimension = stk_meta.spatial_dimension();
      double length = 0.0;
      for (UInt dim=0; dim<spatial_dimension; ++dim)
      {
        length += (coord0[dim]-coord1[dim])*(coord0[dim]-coord1[dim]);
      }
      length = std::sqrt(length);

      const double ls0 = coord0[0] - location;
      const double ls1 = coord1[0] - location;

      if ((ls0 > threshold && ls1 > threshold) || (ls0 < -threshold && ls1 < -threshold))
      {
        return percept::DO_UNREFINE;
      }
      else
      {
        if (length > max_length)
        {
          return percept::DO_REFINE;
        }
      }
      return percept::DO_NOTHING;
    }

  }

#endif
