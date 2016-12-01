// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <adapt/IElementAdapter.hpp>
#include <percept/FieldTypes.hpp>

  namespace percept {

    void IElementAdapter::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element,
                                            vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      //CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

      int markInfo = markElement(element);

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_rank == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_rank == stk::topology::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          bool needNodes = ( markInfo & DO_REFINE);
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  bool doMark = true;
                  if (needed_entity_ranks[ineed_ent].third.size())
                    {
                      VERIFY_OP_ON(needed_entity_ranks[ineed_ent].third.size(), ==, numSubDimNeededEntities, "bad size");
                      if (!needed_entity_ranks[ineed_ent].third[iSubDimOrd])
                        doMark = false;
                    }
                  if (doMark)
                    (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, needNodes,bucket_topo_data);
                } // iSubDimOrd
            }
        } // ineed_ent
    }


    void IElementAdapter::buildUnrefineList(ElementUnrefineCollection& elements_to_unref)
    {
      //ElementUnrefineCollection elements_to_unref(*m_eMesh.get_bulk_data());
      elements_to_unref.clear();

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];

                // FIXME
                // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error is isParentElement)
                const bool check_for_family_tree = false;
                bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);

                if (isParent)
                  continue;

                const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

                if (elem_nodes.size())
                  {
                    int markInfo = markElement(element);
                    if (markInfo & DO_UNREFINE)
                      {
                        elements_to_unref.insert(element);
                      }
                  }
              }
          }
        }

      //return elements_to_unref;
    }

  }


