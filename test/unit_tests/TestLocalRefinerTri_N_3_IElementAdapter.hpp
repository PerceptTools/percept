// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TestLocalRefinerTri_N_3_IElementAdapter_hpp
#define adapt_TestLocalRefinerTri_N_3_IElementAdapter_hpp

#include <adapt/IElementAdapter.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation as a use case for IElementAdapter
     */
    class TestLocalRefinerTri_N_3_IElementAdapter : public IElementAdapter
    {
    public:
      TestLocalRefinerTri_N_3_IElementAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      //virtual ElementUnrefineCollection  buildUnrefList();

    protected:

      /// Client supplies these methods - given an element, which edge, and the nodes on the edge, return instruction on what to do to the edge,
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int markElement(const stk::mesh::Entity element);

    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_3_IElementAdapter in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTri_N_3_IElementAdapter::TestLocalRefinerTri_N_3_IElementAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      IElementAdapter(eMesh, bp, proc_rank_field)
    {
    }


    int TestLocalRefinerTri_N_3_IElementAdapter::
    markElement(const stk::mesh::Entity element)
    {
      int mark=0;

      // refine test
      {
        // vertical line position
        const double vx = 0.21;

        // horizontal line position
        const double vy = 1.21;

        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        for (unsigned inode=0; inode < elem_nodes.size(); inode++)
          {
            stk::mesh::Entity node = elem_nodes[inode].entity();
            stk::mesh::Entity node1 = elem_nodes[(inode + 1) % 3].entity();
            double *coord0 = stk::mesh::field_data( *m_eMesh.get_coordinates_field(), node );
            double *coord1 = stk::mesh::field_data( *m_eMesh.get_coordinates_field(), node1 );

            // choose to refine or not
            if (
                ( std::fabs(coord0[0]-coord1[0]) > 1.e-3 &&
                  ( (coord0[0] < vx && vx < coord1[0]) || (coord1[0] < vx && vx < coord0[0]) )
                  )
                ||
                ( std::fabs(coord0[1]-coord1[1]) > 1.e-3 &&
                  ( (coord0[1] < vy && vy < coord1[1]) || (coord1[1] < vy && vy < coord0[1]) )
                  )
                )
              {
                mark |= DO_REFINE;
                break;
              }
          }
      }

      // unrefine test
      {
        const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

        bool found = true;
        for (unsigned inode=0; inode < elem_nodes.size(); inode++)
          {
            stk::mesh::Entity node = elem_nodes[inode].entity();
            double *coord = stk::mesh::field_data( *m_eMesh.get_coordinates_field(), node );
            //if (coord[0] > 2.1 || coord[1] > 2.1)
            if (coord[0] > 1.0001 || coord[1] > 1.0001)
              {
                found = false;
                break;
              }
          }

        if (found)
          {
            mark |= DO_UNREFINE;
          }
      }
      return mark;
    }

  }

#endif
