// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_NIDQuantum_hpp
#define adapt_NIDQuantum_hpp

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

  namespace percept {

    class NIDQuantum {
      static percept::PerceptMesh *s_pMesh;
      stk::mesh::EntityId m_id;
    public:
      NIDQuantum(stk::mesh::EntityId id = 0ULL) : m_id(id) {
        s_pMesh = percept::PerceptMesh::get_static_instance();
        if (!s_pMesh) throw std::runtime_error("no get_static_instance...");
      }
      NIDQuantum(const NIDQuantum& nid) : m_id(nid.m_id) {
      }
      NIDQuantum(const stk::mesh::Entity& entity) : m_id(s_pMesh->identifier(entity)) {
      }

      stk::mesh::EntityId get_id() { return m_id; }
      const stk::mesh::EntityId get_id() const { return m_id; }
      void set_id(stk::mesh::EntityId id) { m_id = id; }
      stk::mesh::Entity get_entity() { return s_pMesh->get_bulk_data()->get_entity(s_pMesh->node_rank(), m_id); }
      const stk::mesh::Entity get_entity() const { return s_pMesh->get_bulk_data()->get_entity(s_pMesh->node_rank(), m_id); }
      //const stk::mesh::Entity get_entity() const { return s_pMesh->get_bulk_data()->get_entity(s_pMesh->node_rank(), m_id); }
      operator stk::mesh::Entity () { return get_entity(); }
      operator const stk::mesh::Entity () const { return get_entity(); }
      NIDQuantum& operator=(const stk::mesh::Entity& entity) {
        m_id = s_pMesh->identifier(entity);
        return *this;
      }
      NIDQuantum& operator=(const stk::mesh::EntityId& id) {
        m_id = id;
        return *this;
      }
      bool operator==(const NIDQuantum& b) const { return m_id == b.get_id(); }
      bool operator!=(const NIDQuantum& b) const { return m_id != b.get_id(); }
      bool operator<(const NIDQuantum& b) const { return m_id < b.get_id(); }

      size_t local_offset() const { return get_entity().local_offset(); }

    };

    std::ostream & operator << ( std::ostream & os, const NIDQuantum & nidq);
  }
#endif
