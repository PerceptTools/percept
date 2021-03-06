// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_TestLocalRefinerTet_N_3_hpp
#define adapt_TestLocalRefinerTet_N_3_hpp

#include <adapt/Refiner.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    /**
     * A test implementation that marks some edges randomly to test RefinerPattern_Tri3_Tri3_N_3
     */
    class TestLocalRefinerTet_N_3 : public Refiner
    {
    public:
      TestLocalRefinerTet_N_3(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

      // ElementUnrefineCollection  buildTestUnrefineList();

    protected:


      virtual void 
      refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, 
            vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data);


    };

    // This is a very specialized test that is used in unit testing only (see unit_localRefiner/break_tri_to_tri_N_3 in UnitTestLocalRefiner.cpp)

    TestLocalRefinerTet_N_3::TestLocalRefinerTet_N_3(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      Refiner(eMesh, bp, proc_rank_field)
    {
    }


    void TestLocalRefinerTet_N_3::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
    {
      //static int n_seq = 400;

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
#if 0
                  stk::mesh::Entity node0 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
                  stk::mesh::Entity node1 = *elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
                  double * const coord0 =stk::mesh::field_data( *coordField , node0 );
                  double * const coord1 =stk::mesh::field_data( *coordField , node1 );
                  
                  // vertical line position
                  const double vx = 0.21;

                  // horizontal line position
                  const double vy = 1.21;

                  // choose to refine or not 
                  if (
                      ( std::fabs(coord0[0]-coord1[0]) > 1.e-3 &&
                        ( (coord0[0] < vx && vx < coord1[0]) || (coord1[0] < vx && vx < coord0[0]) )
                        )
                      ||
                      ( std::fabs(coord0[1]-coord1[1]) > 1.e-3 &&
                        ( (coord0[1] < vy && vy < coord1[1]) || (coord1[1] < vy && vy < coord0[1]) )
                        )
                      )
                    {
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd,  true,bucket_topo_data);
                    }

#endif
                  // mark first m_edge_mark_bitcode edges 
                  unsigned edge_mark_bitcode = (m_eMesh.identifier(element) - 1);

                  if ( ((1 << iSubDimOrd) & edge_mark_bitcode ) )
                    {
                      //if (edge_mark_bitcode == 4)
                      if (0)
                        {
                          std::cout << "tmp TestLocalRefinerTet_N_3 m_eMesh.identifier(element) = " << m_eMesh.identifier(element) 
                                    << " edge_mark_bitcode = " << edge_mark_bitcode << "  iSubDimOrd= " << iSubDimOrd << std::endl;
                        }
                      (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd,  true,bucket_topo_data);
                    }

                }

            } // iSubDimOrd
        } // ineed_ent
    }


  }

#endif
