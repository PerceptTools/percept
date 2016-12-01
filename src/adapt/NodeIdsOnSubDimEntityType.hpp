// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_NodeIdsOnSubDimEntityType_hpp
#define adapt_NodeIdsOnSubDimEntityType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include <stk_mesh/base/Entity.hpp>
#include <percept/stk_mesh.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <percept/NoMallocArray.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>

#include <percept/PerceptBoostArray.hpp>
#include <adapt/SubDimCell.hpp>
#include <adapt/NIDQuantum.hpp>

  namespace percept {


    // type defining what is stored on the edge
    typedef stk::mesh::Entity NodeIdsOnSubDimEntityTypeQuantum;
    //typedef NIDQuantum NodeIdsOnSubDimEntityTypeQuantum;

    /// data on a sub-dim entity (global node ids on the entity, the owning element's id)
    // FIXME - don't inherit from vector
    struct NodeIdsOnSubDimEntityType : public std::vector<NodeIdsOnSubDimEntityTypeQuantum>
    {
      typedef std::vector<NodeIdsOnSubDimEntityTypeQuantum> base_type;
      typedef std::vector<stk::mesh::EntityId> entity_id_vector_type;
      entity_id_vector_type m_entity_id_vector;
      unsigned m_mark;

      NodeIdsOnSubDimEntityType(unsigned sz=1, NodeIdsOnSubDimEntityTypeQuantum allValues = NodeIdsOnSubDimEntityTypeQuantum(),
                                unsigned mark=0u) : base_type(sz,allValues),
                                                    m_entity_id_vector(sz,0u),
                                                    m_mark(mark)  {}
      void resize(size_t sz)
      {
        m_entity_id_vector.resize(sz);
        base_type::resize(sz);
      }

      void pack(percept::PerceptMesh& eMesh, stk::CommBuffer& buff)
      {
        buff.pack< unsigned > ( this->size() );
        m_entity_id_vector.resize( this->size() );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            if ( ! eMesh.is_valid((*this)[ii]) )
              throw std::logic_error("logic err in NodeIdsOnSubDimEntityType::pack");
            stk::mesh::EntityId id = eMesh.identifier((*this)[ii]);
            VERIFY_OP_ON(id, != , 0, "logic err 2 in NodeIdsOnSubDimEntityType::pack");
            m_entity_id_vector[ii] = id;
            buff.pack<stk::mesh::EntityId>( id );
          }
      }
      void unpack(percept::PerceptMesh& eMesh, stk::CommBuffer& buff)
      {
        unsigned sz;
        buff.unpack< unsigned > ( sz );
        this->resize( sz );
        m_entity_id_vector.resize( sz );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            stk::mesh::EntityId id=0;
            buff.unpack<stk::mesh::EntityId>( id );
            m_entity_id_vector[ii] = id;
            VERIFY_OP_ON(id, != , 0, "logic err 3 in NodeIdsOnSubDimEntityType::unpack");
          }
      }
    };

  }

#endif
