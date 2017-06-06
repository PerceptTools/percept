// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TEA_SelectIfRefined_hpp
#define adapt_TEA_SelectIfRefined_hpp

#include <adapt/Refiner.hpp>


namespace percept {

  class TEA_SelectIfRefined : public RefinerSelector
  {

  public:
    PerceptMesh& m_eMesh;
    //RefineFieldType *refine_field = eMesh.get_fem_meta_data()-> template get_field<RefineFieldType>(stk::topology::ELEMENT_RANK, "refine_field");
    SetOfEntities m_entities[4]; // elements or sides or edges
    TEA_SelectIfRefined(PerceptMesh& eMesh) : m_eMesh(eMesh)
    {
      //inititialize();
    }

    void initialize()
    {
      PerceptMesh& eMesh = m_eMesh;
      SetOfEntities node_set;  // participating nodes - elems or sides touching these stay in the list

      for (stk::mesh::EntityRank rank = eMesh.edge_rank(); rank <= eMesh.element_rank(); ++rank)
        {
          m_entities[rank].clear();
        }

      std::vector<stk::mesh::Entity> non_participating_elems;

      // build list of participating nodes - include ghost elements, also possible list of non-participating elements
      {
        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity entity = bucket[ientity];
                int *refine_field_elem = stk::mesh::field_data( *m_eMesh.m_refine_field , entity );
                if (refine_field_elem[0] > 0)
                  {
                    m_entities[eMesh.element_rank()].insert(entity);
                    unsigned nnode= eMesh.get_bulk_data()->num_nodes(entity);
                    stk::mesh::Entity const *elem_nodes = eMesh.get_bulk_data()->begin_nodes(entity);
                    for (unsigned ii=0; ii < nnode; ii++)
                      {
                        node_set.insert(elem_nodes[ii]);
                      }
                  }
                else
                  {
                    non_participating_elems.push_back(entity);
                  }
              }
          }
      }

      // non participating elements
      {
        for (stk::mesh::Entity entity : non_participating_elems)
          {
            unsigned nnode= eMesh.get_bulk_data()->num_nodes(entity);
            stk::mesh::Entity const *elem_nodes = eMesh.get_bulk_data()->begin_nodes(entity);
            for (unsigned ii=0; ii < nnode; ii++)
              {
                if (node_set.find(elem_nodes[ii]) != node_set.end())
                  {
                    m_entities[eMesh.element_rank()].insert(entity);
                    break;
                  }
              }
          }
      }

      // edges and sides...
      {
        for (stk::mesh::EntityRank rank = eMesh.edge_rank(); rank <= eMesh.side_rank(); ++rank)
          {
            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank );
            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_entity_in_bucket = bucket.size();
                for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                  {
                    stk::mesh::Entity entity = bucket[ientity];
                    unsigned nnode= eMesh.get_bulk_data()->num_nodes(entity);
                    stk::mesh::Entity const *elem_nodes = eMesh.get_bulk_data()->begin_nodes(entity);
                    for (unsigned ii=0; ii < nnode; ii++)
                      {
                        if (node_set.find(elem_nodes[ii]) != node_set.end())
                          {
                            m_entities[rank].insert(entity);
                            break;
                          }
                      }
                  }
              }
          }
      }
    }

    virtual bool use_batch_filter() { return true; }
    virtual void batch_filter(stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& elements)
    {
      elements.assign(m_entities[rank].begin(), m_entities[rank].end());
    }
  };
}

#endif
