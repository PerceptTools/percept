// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
#include <adapt/geometry/GeometryKernelOpenNURBS.hpp>
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
                     bool only_count, bool doAllElements)
    //bool only_count=false, bool doAllElements=true)
    {
      EXCEPTWATCH;
      unsigned num_elem = 0;

      int progress_meter_num_total = 0;
      if (m_doProgress)
        {
          m_doProgress = false;
          progress_meter_num_total = doForAllElements(irank, function_info, rank, function, elementType, needed_entity_ranks,  true, doAllElements);
          m_doProgress = true;
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, "NodeRegistry passes");
          notifyObservers(&pd);
        }
      int progress_meter_when_to_post = progress_meter_num_total / m_progress_meter_frequency;
      if (0 == progress_meter_when_to_post)
        progress_meter_when_to_post = 1;
      double d_progress_meter_num_total = progress_meter_num_total;

      percept::PerceptMesh& eMesh = m_eMesh;

#if 0
      stk::mesh::MetaData& metaData = *eMesh.get_fem_meta_data();
      const std::vector< stk::mesh::Part * > & parts = metaData.get_parts();
      unsigned nparts = parts.size();
      if (1) std::cout << "Number of parts = " << nparts << std::endl;
      CoordinatesFieldType* coordField = eMesh.get_coordinates_field();
#endif
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

                if (!only_count && (doAllElements || elementIsGhost))
                  {
                    refineMethodApply(function, element, needed_entity_ranks, cell_topo_data);
                  }

                if (m_doProgress && (num_elem % progress_meter_when_to_post == 0) )
                  {
                    double progress_meter_percent = 100.0*((double)num_elem)/d_progress_meter_num_total;
                    ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, "NodeRegistry passes");
                    notifyObservers(&pd);
                  }
              }

          }
        }

      if (m_doProgress)
        {
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, "NodeRegistry passes");
          notifyObservers(&pd);
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
