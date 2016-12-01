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

#include <unit_tests/TestLocalRefinerTri_N.hpp>

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


    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTri_N::TestLocalRefinerTri_N(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      Refiner(eMesh, bp, proc_rank_field)
    {
    }

    void TestLocalRefinerTri_N::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      //CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

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

          // see how many edges are already marked
          int num_marked=0;
          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
                {
                  bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
                  if (!is_empty) ++num_marked;
                }
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              //SubDimCell_SDCEntityType subDimEntity;
              //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
              //bool is_empty = m_nodeRegistry->is_empty( element, needed_entity_rank, iSubDimOrd);
              //if(1||!is_empty)

              if (needed_entity_rank == m_eMesh.edge_rank())
                {
                  stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  //double * const coord0 =stk::mesh::field_data( *coordField , node0 );
                  //double * const coord1 =stk::mesh::field_data( *coordField , node1 );

                  // choose to refine or not
                  //if (0 || (m_eMesh.identifier(node0) < m_eMesh.identifier(node1) || (m_eMesh.identifier(node0) % 2 == 0)))
                  if (m_eMesh.identifier(node0) + m_eMesh.identifier(node1) < 40)
                    {
                      int delta = (int)m_eMesh.identifier(node0) - (int)m_eMesh.identifier(node1);
                      if (delta < 0) delta = -delta;
                      if ( delta > 4 && (m_eMesh.identifier(node0) + m_eMesh.identifier(node1) < 30))
                        {
                          (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true,bucket_topo_data);
                        }
                      if ( (m_eMesh.identifier(node0) + m_eMesh.identifier(node1) > 20))
                        {
                          (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true,bucket_topo_data);
                        }
                    }
                }

            } // iSubDimOrd
        } // ineed_ent
    }

    void TestLocalRefinerTri_N::buildTestUnrefineList(ElementUnrefineCollection& elements_to_unref)
    {
      //ElementUnrefineCollection elements_to_unref(*m_eMesh.get_bulk_data());
      elements_to_unref.clear();

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_entity_in_bucket = bucket.size();
            for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
              {
                stk::mesh::Entity element = bucket[ientity];
                bool elementIsGhost = m_eMesh.isGhostElement(element);
                if (elementIsGhost)
                  continue;

                const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

                if (elem_nodes.size() && m_eMesh.isChildElement(element))
                  {
                    bool found = true;
                    for (unsigned inode=0; inode < elem_nodes.size(); inode++)
                      {
                        stk::mesh::Entity node = elem_nodes[inode].entity();
                        double *coord = stk::mesh::field_data( *m_eMesh.get_coordinates_field(), node );
                        if (coord[0] > 2.1 || coord[1] > 2.1)
                          {
                            found = false;
                            break;
                          }
                      }
                    if (found)
                      {
                        elements_to_unref.insert(element);
                        //std::cout << "tmp element id= " << m_eMesh.identifier(element) << " ";
                        //m_eMesh.print_entity(std::cout, element);
                      }
                  }
              }
          }
        }

      //return elements_to_unref;
    }


  } // namespace percept

