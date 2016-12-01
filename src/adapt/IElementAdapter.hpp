// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_IElementAdapter_hpp
#define adapt_IElementAdapter_hpp

#include <adapt/IAdapter.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * An IElementAdapter is an abstract base class for derived classes that are required to overload the mark method,
     *   which supplies the derived class with the element to be marked for refine, unrefine, or both (@see IAdapter::AdaptInstruction)
     */
    class IElementAdapter : public IAdapter
    {
    public:
      IElementAdapter(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0)
      : IAdapter(eMesh, bp, proc_rank_field) {}

      virtual void  buildUnrefineList(ElementUnrefineCollection&) ;

    protected:

      /// Client supplies this method - given an element return instruction on what to do to the element:
      ///    DO_NOTHING (nothing), DO_REFINE (refine), DO_UNREFINE
      virtual int markElement(const stk::mesh::Entity element) = 0;

      virtual void
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
                                              vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data);

    };
  }

#endif
