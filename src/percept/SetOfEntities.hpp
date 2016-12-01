// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SetOfEntities_hpp
#define SetOfEntities_hpp

#define PERCEPT_USE_STD_SET 1
#define PERCEPT_USE_BOOST_SET 0

#if PERCEPT_USE_BOOST_SET
#include <boost/unordered_set.hpp>
#endif


  namespace percept {

#if PERCEPT_USE_STD_SET
    typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> SetOfEntities;
#endif
#if PERCEPT_USE_BOOST_SET
    struct myEsetHash : public std::unary_function< stk::mesh::Entity, std::size_t>
    {
      stk::mesh::BulkData * m_bulk;
      myEsetHash(stk::mesh::BulkData& bulk) : m_bulk(&bulk) {}
      std::size_t operator()(const stk::mesh::Entity& e) const { return m_bulk->identifier(e); }
    };
    struct myEsetEqual : public std::binary_function<stk::mesh::Entity, stk::mesh::Entity, bool>
    {
      stk::mesh::BulkData * m_bulk;
      myEsetEqual(stk::mesh::BulkData& bulk):m_bulk(&bulk) {}
      bool operator()(const stk::mesh::Entity& e1, const stk::mesh::Entity& e2) const { return e1 == e2; }
    };

#if 0
    typedef boost::unordered_set<stk::mesh::Entity, myEsetHash, myEsetEqual> SetOfEntitiesBase;
    struct SetOfEntities : public SetOfEntitiesBase
    {
      //stk::mesh::BulkData& m_bulk;
      SetOfEntities(stk::mesh::BulkData& bulk) : SetOfEntitiesBase(0, myEsetHash(bulk), myEsetEqual(bulk))
      {
      }
    };
#else
    typedef boost::unordered_set<stk::mesh::Entity> SetOfEntitiesBase;
    struct SetOfEntities : public SetOfEntitiesBase
    {
      //stk::mesh::BulkData& m_bulk;
      SetOfEntities(stk::mesh::BulkData& bulk) : SetOfEntitiesBase()
      {
      }
    };
#endif

#endif
    //typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> ElementUnrefineCollection;
    //typedef elements_to_be_destroyed_type ElementUnrefineCollection;
  }

#endif
