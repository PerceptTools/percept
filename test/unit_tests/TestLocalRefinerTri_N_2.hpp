// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TestLocalRefinerTri_N_2_hpp
#define adapt_TestLocalRefinerTri_N_2_hpp

#include <adapt/Refiner.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that marks some edges randomly to test RefinerPattern_Tri3_Tri3_N_2
     */
    class TestLocalRefinerTri_N_2 : public Refiner
    {
    public:
      TestLocalRefinerTri_N_2(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      void  buildTestUnrefineList(ElementUnrefineCollection&);

    protected:

      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, 
            vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data);


    };



  }

#endif
