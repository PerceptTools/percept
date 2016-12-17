// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefinerPattern_Wedge18_Wedge18_8_sierra_hpp
#define adapt_UniformRefinerPattern_Wedge18_Wedge18_8_sierra_hpp

#include <adapt/sierra_element/RefinementTopology.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>

#define FACE_BREAKER_W18_W18_8 0
#if FACE_BREAKER_W18_W18_8
#include "UniformRefinerPattern_Quad8_Quad8_4_sierra.hpp"
#include "UniformRefinerPattern_Tri6_Tri6_4_sierra.hpp"
#endif

#include <percept/PerceptBoostArray.hpp>

  namespace percept {

    template <>
    class UniformRefinerPattern<shards::Wedge<18>, shards::Wedge<18>, 8, SierraPort > : public URP<shards::Wedge<18>,shards::Wedge<18>  >
    {

#if FACE_BREAKER_W18_W18_8
      UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > * m_face_breaker;
      UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > * m_face_breaker_tri;
#endif

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Wedge<18>, shards::Wedge<18>  >(eMesh)
      {
        m_primaryEntityRank = stk::topology::ELEMENT_RANK;

        setNeededParts(eMesh, block_names, true);
        Elem::StdMeshObjTopologies::bootstrap();

#if FACE_BREAKER_W18_W18_8

        m_face_breaker =  new UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > (eMesh, block_names) ;
        m_face_breaker_tri = new UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > (eMesh, block_names);
#endif

      }
      ~UniformRefinerPattern()
      {
#if FACE_BREAKER_W18_W18_8
        if (m_face_breaker) delete m_face_breaker;
        if (m_face_breaker_tri) delete m_face_breaker_tri;
#endif
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;

#if FACE_BREAKER_W18_W18_8
        bp.resize(3);
        bp[0] = this;
        bp[1] = m_face_breaker;
        bp[2] = m_face_breaker_tri;
#else
        bp.resize(1);
        bp[0] = this;
#endif
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        // FIXME need to take into account the mixed topology nature of the faces
        needed_entities.resize(3);

        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 9u);
        needed_entities[2] = NeededEntityType(stk::topology::ELEMENT_RANK, 18u);

      }

      virtual unsigned getNumNewElemPerElem() { return 8; }


      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        vector<stk::mesh::Entity>::iterator& ft_element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
#if 0
        static bool s_not_printed = true;

        if (s_not_printed)
          {
            s_not_printed = false;
            printRefinementTopoX_Table();
          }
#endif
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool, ft_element_pool, 
                                        proc_rank_field);
      }

    };

  }

#endif