// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <unit_tests/TestLocalRefiner.hpp>

//#define STK_PERCEPT_HAS_GEOMETRY
#undef STK_PERCEPT_HAS_GEOMETRY
#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <adapt/geometry/MeshGeometry.hpp>
#include <adapt/geometry/GeometryFactory.hpp>
#endif

// FIXME
// #include <stk_mesh/baseImpl/EntityImpl.hpp>
// #include <stk_mesh/base/Entity.hpp>
// FIXME

  namespace percept {


    TestLocalRefiner::TestLocalRefiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      Refiner(eMesh, bp, proc_rank_field)
    {
    }

    // test uniform refinement using bulk data and buckets directly (not using the element color vectors)
    unsigned TestLocalRefiner::
    doForAllElements(unsigned irank, std::string function_info,
                     stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
                     unsigned elementType,
                     vector<NeededEntityType>& needed_entity_ranks,
                     bool doAllElements)
    {
      EXCEPTWATCH;
      unsigned num_elem = 0;

      percept::PerceptMesh& eMesh = m_eMesh;

      stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();

      const stk::mesh::BucketVector & buckets = bulkData.buckets( rank );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (in_surface_selector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            // in case the cell topology is needed
            const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(bucket);
            shards::CellTopology cell_topo(cell_topo_data);

            if (cell_topo.getKey() != elementType)
              continue;

            const unsigned num_elements_in_bucket = bucket.size();

            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                const stk::mesh::Entity element = bucket[iElement];

                //const stk::mesh::PairIterRelation& elem_nodes = element.relations( stk::topology::NODE_RANK );

                //const stk::mesh::Entity element = * element_p;

                bool elementIsGhost = m_eMesh.isGhostElement(element);
                if (!elementIsGhost)
                  ++num_elem;

                if (doAllElements || elementIsGhost)
                  {
                    refineMethodApply(function, element, needed_entity_ranks, cell_topo_data);
                  }
              }
          }
        }

      return num_elem;
    }

    void TestLocalRefiner::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

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

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              //SubDimCell_SDCEntityType subDimEntity;
              //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
              //bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
              //if(1||!is_empty)

              (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true,bucket_topo_data);

            } // iSubDimOrd
        } // ineed_ent
    }


  } // namespace percept
