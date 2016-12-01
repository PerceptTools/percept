// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TestLocalRefiner_hpp
#define adapt_TestLocalRefiner_hpp

#include <adapt/Refiner.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that does uniform refinement but uses non-uniform methods
     */
    class TestLocalRefiner : public Refiner
    {
    public:
      TestLocalRefiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

    protected:
      virtual unsigned
      doForAllElements(unsigned irank, std::string function_info,
                       stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
                       unsigned elementType,
                       vector<NeededEntityType>& needed_entity_ranks,
                       bool only_count=false, bool doAllElements=true);


      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
            vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data);


    };



  }

#endif
