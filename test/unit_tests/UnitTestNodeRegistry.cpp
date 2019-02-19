// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_io/IossBridge.hpp>

#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/fixtures/Fixture.hpp>

#include <adapt/NodeRegistry.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/SerializeNodeRegistry.hpp>
#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <stk_util/parallel/Parallel.hpp>
#include <math.h>
#include <stk_mesh/base/MeshUtils.hpp>

namespace percept {
namespace unit_tests {

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================

TEST(nodeRegistry, createAddNodes_serial_and_1st_parallel)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  percept::PerceptMesh eMesh(3u);
  eMesh.new_mesh(percept::GMeshSpec("3x3x12|bbox:0,0,0,1,1,1"));  // create a 3x3x12 hex mesh in the unit cube
  if (1)
    {
      stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part_with_topology("refine_new_nodes_part", stk::topology::NODE);
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      meta.declare_attribute_no_delete(part, &percept::auto_part);
    }
  //int scalarDimension = 0; // a scalar
  //int vectorDimension = 3;
  eMesh.commit();
  eMesh.print_info();

  unsigned p_size = eMesh.get_bulk_data()->parallel_size();
  unsigned p_rank = eMesh.get_bulk_data()->parallel_rank();

  std::cout << "TEST::nodeRegistry::createAddNodes_serial: p_size = "<< p_size << " rank= " << p_rank << std::endl;
  if (p_size >= 1)  // FIXME
    return;

  if (p_size >= 3)
    return;

  //stk::CommSparse comm_all(eMesh.get_bulk_data()->parallel());
  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size == 2) return; // FIXME
  if (p_size == 2)
  {
    unsigned elem_num=(12/p_size)*3*3;

    const stk::mesh::Entity element_1_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num);
    dw() << "P["<<p_rank<<"] element_1_p = " << eMesh.identifier(element_1_p) << DWENDL;
    stk::mesh::EntityRank stk_mesh_Edge = stk::topology::EDGE_RANK;

    NeededEntityType entity_rank(stk_mesh_Edge, 1u);
    unsigned iSubDimOrd = 0u;

    const stk::mesh::Entity element_1 = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num);

    //nodeRegistry.getSubDimEntity(subDimEntity, element_1, entity_rank, iSubDimOrd);
    nodeRegistry.registerNeedNewNode(element_1, entity_rank, iSubDimOrd, true, 0);

    // must be a ghost on one processor
    int n_check = 0;
    if (eMesh.isGhostElement(element_1))
    {
      std::cout  << "P["<<p_rank<<"] element no= " << elem_num
                 << " TEST::nodeRegistry::createAddNodes_serial: is ghost on rank= " << p_rank << std::endl;
      EXPECT_EQ(nodeRegistry.local_size(), 1u);
      ++n_check;
    }
    else
    {
      std::cout  << "P["<<p_rank<<"] element no= " << elem_num
                 << " TEST::nodeRegistry::createAddNodes_serial: is not ghost on rank= " << p_rank << std::endl;
      EXPECT_EQ(nodeRegistry.local_size(), 1u);
      ++n_check;
    }
    EXPECT_EQ(n_check, 1);
    //Util::pause(true, "tmp");
    //exit(1);
  }
  else if (p_size == 1)
  {
    unsigned elem_num=1;

    const stk::mesh::Entity element_1_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num);
    std::cout << "P["<<p_rank<<"] element_1_p = " << element_1_p << std::endl;
    stk::mesh::EntityRank stk_mesh_Edge = stk::topology::EDGE_RANK;

    NeededEntityType entity_rank(stk_mesh_Edge, 1u);
    unsigned iSubDimOrd = 0u;

    const stk::mesh::Entity element_1 = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num);

    nodeRegistry.registerNeedNewNode(element_1, entity_rank, iSubDimOrd, true, 0);

    // must be a ghost on one processor
    if (eMesh.isGhostElement(element_1))
    {
      std::cout << "TEST::nodeRegistry::createAddNodes_serial: is ghost on rank= " << p_rank << std::endl;
      EXPECT_EQ(nodeRegistry.local_size(), 0u);
    }
    else
    {
      std::cout << "TEST::nodeRegistry::createAddNodes_serial: is not ghost on rank= " << p_rank << std::endl;
      EXPECT_EQ(nodeRegistry.local_size(), 1u);
    }
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

TEST(nodeRegistry, test_parallel_0)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  /* output snippets from this test - note that edge 4 of elem 54 is same as edge 0 of elem 63 (used in next test)

     num_elements_in_bucket = 9 element ids =
     46 47 48 49 50 51 52 53 54
     num_elements_in_bucket = 54 element ids =
     55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74
     75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94
     95 96 97 98 99 100 101 102 103 104 105 106 107 108

     P[1] element_local = Elem: 63 nodes: 107 108 112 111 123 124 128 127
     P[1] element_ghost = Elem: 54 nodes: 91 92 96 95 107 108 112 111
     P[1] element_local edge 0 = 107 108
     P[1] element_local edge 1 = 108 112
     P[1] element_local edge 2 = 111 112
     P[1] element_local edge 3 = 107 111
     P[1] element_local edge 4 = 123 124
     P[1] element_local edge 5 = 124 128
     P[1] element_local edge 6 = 127 128
     P[1] element_local edge 7 = 123 127
     P[1] element_local edge 8 = 107 123
     P[1] element_local edge 9 = 108 124
     P[1] element_local edge 10 = 112 128
     P[1] element_local edge 11 = 111 127

     P[0] element_local = Elem: 54 nodes: 91 92 96 95 107 108 112 111
     P[0] element_ghost = Elem: 63 nodes: 107 108 112 111 123 124 128 127
     P[0] element_local edge 0 = 91 92
     P[0] element_local edge 1 = 92 96
     P[0] element_local edge 2 = 95 96
     P[0] element_local edge 3 = 91 95
     P[0] element_local edge 4 = 107 108
     P[0] element_local edge 5 = 108 112
     P[0] element_local edge 6 = 111 112
     P[0] element_local edge 7 = 107 111
     P[0] element_local edge 8 = 91 107
     P[0] element_local edge 9 = 92 108
     P[0] element_local edge 10 = 96 112
     P[0] element_local edge 11 = 95 111

  */

  // start_demo_nodeRegistry_test_parallel_0
  percept::PerceptMesh eMesh(3u);
  eMesh.new_mesh(percept::GMeshSpec("3x3x12|bbox:0,0,0,1,1,1"));  // create a 3x3x12 hex mesh in the unit cube
  if (1)
    {
      stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part_with_topology("refine_new_nodes_part", stk::topology::NODE);
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      meta.declare_attribute_no_delete(part, &percept::auto_part);
    }
  eMesh.commit();
  eMesh.print_info();

  unsigned p_size = eMesh.get_parallel_size();
  unsigned p_rank = eMesh.get_rank();
  Util::setRank(eMesh.get_rank());

  if (p_size != 2) // FIXME
    return;

  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size == 2)
  {
    // pick an element on the processor boundary
    unsigned elem_num_local = (12/p_size)*3*3;
    unsigned elem_num_ghost = elem_num_local+(3*3);

    stk::mesh::Entity element_local_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local);
    stk::mesh::Entity element_ghost_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost);

    if (p_rank == 1)
    {
      element_local_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost);
      element_ghost_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local);
    }

    dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
    dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

    const stk::mesh::Entity element_local = element_local_p;
    const stk::mesh::Entity element_ghost = element_ghost_p;

    std::cout << "P["<<p_rank<<"] element_local = " << element_local << std::endl;
    std::cout << "P["<<p_rank<<"] element_ghost = " << element_ghost << std::endl;

    stk::mesh::EntityRank stk_mesh_Edge = stk::topology::EDGE_RANK;

    stk::mesh::EntityRank needed_entity_rank = stk_mesh_Edge;
    for (unsigned iSubDimOrd = 0; iSubDimOrd < 12; iSubDimOrd++)
    {
      SubDimCell_SDCEntityType subDimEntity(&eMesh);
      nodeRegistry.getSubDimEntity(subDimEntity, element_local, needed_entity_rank, iSubDimOrd);
      std::cout << "P[" << p_rank << "] element_local edge " << iSubDimOrd << " = " << subDimEntity << std::endl;
    }

  }
}


//=============================================================================
//=============================================================================
//=============================================================================

TEST(nodeRegistry, test_parallel_1_0)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  // start_demo_nodeRegistry_test_parallel_1

  percept::PerceptMesh eMesh(3u);

  // unsigned p_size = eMesh.get_parallel_size();
  // unsigned p_rank = eMesh.get_rank();
  // Util::setRank(eMesh.get_rank());
  unsigned p_size = stk::parallel_machine_size(eMesh.parallel());
  unsigned p_rank = stk::parallel_machine_rank(eMesh.parallel());
  Util::setRank(p_rank);

  eMesh.new_mesh(percept::GMeshSpec(std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,1")));

  // prepare for adding some quadratic elements
  stk::mesh::Part& block_hex_20 = eMesh.get_fem_meta_data()->declare_part("block_hex_20", stk::topology::ELEMENT_RANK);
  /// set cell topology for the part block_hex_20
  stk::mesh::set_cell_topology< shards::Hexahedron<20>  >( block_hex_20 );
  stk::io::put_io_part_attribute(block_hex_20);
  if (1)
    {
      stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part_with_topology("refine_new_nodes_part", stk::topology::NODE);
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      meta.declare_attribute_no_delete(part, &percept::auto_part);
    }

  eMesh.commit();
  eMesh.print_info();
  eMesh.save_as("./cube1x1x2_hex-20-orig.e");

  stk::mesh::Part* block_hex_8 = const_cast<stk::mesh::Part *>(eMesh.getPart("block_1"));

  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size <= 2)
  {
    // pick an element on the processor boundary
    unsigned elem_num_local = 1;
    unsigned elem_num_ghost = 2;
    unsigned elem_num_local_hex20 = 3;
    unsigned elem_num_ghost_hex20 = 4;
    if (p_size == 1)
      {
        elem_num_ghost = elem_num_local;
        elem_num_ghost_hex20 = elem_num_local_hex20;
      }

    eMesh.get_bulk_data()->modification_begin();
    if (p_size == 1)
      {
        stk::mesh::Entity elem0 = eMesh.get_bulk_data()->declare_element(elem_num_local_hex20, stk::mesh::ConstPartVector{&block_hex_20});
        for(unsigned node_ord = 0 ; node_ord < 20; ++node_ord)
        {
          stk::mesh::Entity new_node = eMesh.get_bulk_data()->declare_node(node_ord+1000);
          eMesh.get_bulk_data()->declare_relation( elem0 , new_node , node_ord);
        }
      }
    else
      {
        if (p_rank == 0)
        {
          stk::mesh::Entity elem1 = eMesh.get_bulk_data()->declare_element(elem_num_local_hex20, stk::mesh::ConstPartVector{&block_hex_20});
          for(unsigned node_ord = 0 ; node_ord < 20; ++node_ord)
          {
            stk::mesh::Entity new_node = eMesh.get_bulk_data()->declare_node(node_ord+2000);
            eMesh.get_bulk_data()->declare_relation( elem1 , new_node , node_ord);
          }
        }
        if (p_rank == 1)
        {
          stk::mesh::Entity elem2 = eMesh.get_bulk_data()->declare_element(elem_num_ghost_hex20, stk::mesh::ConstPartVector{&block_hex_20});
          for(unsigned node_ord = 0 ; node_ord < 20; ++node_ord)
          {
            stk::mesh::Entity new_node = eMesh.get_bulk_data()->declare_node(node_ord+3000);
            eMesh.get_bulk_data()->declare_relation( elem2 , new_node , node_ord);
          }
        }
      }
    stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
    eMesh.get_bulk_data()->modification_end();

    stk::mesh::Entity element_local = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local);
    stk::mesh::Entity element_ghost = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost);
    stk::mesh::Entity element_local_hex20 = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local_hex20);
    stk::mesh::Entity element_ghost_hex20 = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost_hex20);

    if (p_rank == 1)
    {
      element_local = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost);
      element_ghost = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local);
      element_local_hex20 = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost_hex20);
      element_ghost_hex20 = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local_hex20);
    }

    dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
    dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

    std::cout << "P["<<p_rank<<"] element_local = " << element_local << std::endl;
    std::cout << "P["<<p_rank<<"] element_ghost = " << element_ghost << std::endl;

    // choose edges to be used for new node locations (i.e., this would model a serendipity-like element with only edge Lagrange nodes)
    stk::mesh::EntityRank stk_mesh_Edge = stk::topology::EDGE_RANK;
    NeededEntityType needed_entity_rank( stk_mesh_Edge, 1u);
    std::vector<NeededEntityType> needed_entity_ranks(1, needed_entity_rank);

    /*
     * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
     *   from remote proc if necessary; finally, the local node database is ready to be queried
     *
     * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
     * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
     * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
     */
    nodeRegistry.beginRegistration();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_local, needed_entity_ranks,0);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_ghost, needed_entity_ranks,0);
    nodeRegistry.endRegistration();

    std::cout << "P["<<p_rank<<"] nodeRegistry size  = " << nodeRegistry.total_size() << std::endl;
    std::cout << "P["<<p_rank<<"] nodeRegistry lsize = " << nodeRegistry.local_size() << std::endl;

    dw() << "P["<<p_rank<<"] nodeRegistry size       = " << nodeRegistry.total_size() << DWENDL;
    dw() << "P["<<p_rank<<"] nodeRegistry lsize      = " << nodeRegistry.local_size() << DWENDL;

    // check if the newly requested nodes are local or remote
    nodeRegistry.beginCheckForRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_local, needed_entity_ranks,0);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_ghost, needed_entity_ranks,0);
    nodeRegistry.endCheckForRemote();

    // get the new nodes from other procs if they are nonlocal
    nodeRegistry.beginGetFromRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_local, needed_entity_ranks,0);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_ghost, needed_entity_ranks,0);
    nodeRegistry.endGetFromRemote();

    // now we can get the new node's id and entity
    unsigned iSubDimOrd = 4u;
    if (p_rank)
    {
      iSubDimOrd = 0u;
    }
    NodeIdsOnSubDimEntityType& nodeIds_onSE_0 = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_rank.first, iSubDimOrd));
    stk::mesh::Entity node_0   = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, eMesh.identifier(nodeIds_onSE_0[0]));

    // should be the same node on each proc
    std::cout << "P[" << p_rank << "] nodeId_0 = " << nodeIds_onSE_0 << " node_0= " << node_0 << std::endl;

    // end_demo

    if (p_size == 1)
      {
        if (1) {
          YAML::Emitter out;
          out << YAML::Anchor("NodeRegistry::map");
          out << YAML::BeginMap;
          out << YAML::Key << YAML::BeginSeq << 1 << 2 << YAML::EndSeq << YAML::Value << YAML::BeginSeq << -1 << -2 << YAML::EndSeq;
          out << YAML::Key << 1;
          out << YAML::Value << 2;
          out << YAML::Key << 3;
          out << YAML::Value << 4;
          out << YAML::EndMap;
          std::cout << "out=\n" << out.c_str() << "\n=out" << std::endl;
          //std::string expected_result = "&NodeRegistry::map\n?\n  - 1\n  - 2\n:\n  - -1\n  - -2\n1: 2\n3: 4";
          std::string expected_result = "&NodeRegistry::map\n? - 1\n  - 2\n: - -1\n  - -2\n1: 2\n3: 4";
          std::cout << "out2=\n" << expected_result << std::endl;
          EXPECT_TRUE(expected_result == std::string(out.c_str()));
        }

        YAML::Emitter yaml;
        std::cout << "\nnodeRegistry.serialize_write(yaml)" << std::endl;
        SerializeNodeRegistry::serialize_write(nodeRegistry, yaml, 0);
        //std::cout << yaml.c_str() << std::endl;
        if (!yaml.good())
          {
            std::cout << "Emitter error: " << yaml.good() << " " <<yaml.GetLastError() << "\n";
            EXPECT_TRUE(false);
          }
        std::ofstream file1("out.yaml");
        file1 << yaml.c_str();
        file1.close();
        std::ifstream file2("out.yaml");
        //YAML::Parser parser(file2);
        YAML::Node doc;

        try {
          //while(parser.GetNextDocument(doc)) {
          if (1) {
            doc = YAML::Load(file2);
            std::cout << "\n read doc.Type() = " << doc.Type() << " doc.Tag()= " << doc.Tag() << " doc.size= " << doc.size() << std::endl;
            if (doc.Type() == YAML::NodeType::Map)
              {
                for(YAML::const_iterator it=doc.begin();it!=doc.end();++it) {
                  int key, value;
                  std::cout << "read it.first().Type() = " << it->first.Type() << " it.first().Tag()= " << it->first.Tag() << std::endl;
                  std::cout << "read it.second().Type() = " << it->second.Type() << " it.second().Tag()= " << it->second.Tag() << std::endl;
                  const YAML::Node& keySeq = it->first;
                  for(YAML::const_iterator itk=keySeq.begin();itk!=keySeq.end();++itk) {
                    key = itk->as<int>();
                    std::cout << "read key= " << key << std::endl;
                  }

                  const YAML::Node& valSeq = it->second;
                  for(YAML::const_iterator itv=valSeq.begin();itv!=valSeq.end();++itv) {
                    value = itv->as<int>();
                    std::cout << "read value= " << value << std::endl;
                  }

                }
              }
          }
        }
        catch(YAML::ParserException& e) {
          std::cout << e.what() << "\n";
          EXPECT_TRUE(false);
        }

        file2.close();
        std::ifstream file3("out.yaml");
        NodeRegistry nrNew(eMesh);
        SerializeNodeRegistry::serialize_read(nrNew, file3);
        YAML::Emitter yaml3;
        std::cout << "\nnrNew.serialize_write(yaml3)" << std::endl;
        SerializeNodeRegistry::serialize_write(nrNew, yaml3, 0);
        std::cout << yaml3.c_str() << std::endl;

        //exit(1);
      }

    // start_demo_nodeRegistry_test_parallel_1_quadratic_elem

    // change element to be a serendipity quadratic element
    eMesh.get_bulk_data()->modification_begin();

    //getCellTopologyData< shards::Node  >()
    const CellTopologyData *const cell_topo_data = eMesh.get_cell_topology(block_hex_20);
    CellTopology cell_topo(cell_topo_data);

    MyPairIterRelation elem_nodes(eMesh, element_local, stk::topology::NODE_RANK);
    //MyPairIterRelation elem_nodes_ghost(eMesh, element_ghost, stk::topology::NODE_RANK);
    for (unsigned isd = 0; isd < 8; isd++)
    {
      eMesh.get_bulk_data()->declare_relation(element_local_hex20, elem_nodes[isd].entity(), isd);
      //eMesh.get_bulk_data()->declare_relation(element_ghost_hex20, elem_nodes_ghost[isd], stk::topology::NODE_RANK);
    }

    for (unsigned isd = 0; isd < 12; isd++)
    {
      nodeRegistry.prolongateCoords(element_local, needed_entity_rank.first, isd);
      NodeIdsOnSubDimEntityType& nodeIds_onSE_0_loc = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_rank.first, isd));

      stk::mesh::Entity node   = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, eMesh.identifier(nodeIds_onSE_0_loc[0]));

      unsigned edge_ord = 8u + isd;
      //unsigned n_edge_ord = cell_topo_data->edge[isd].topology->node_count;
      //std::cout << "n_edge_ord = " << n_edge_ord << std::endl;
      edge_ord = cell_topo_data->edge[isd].node[2];
      eMesh.get_bulk_data()->declare_relation(element_local_hex20, node, edge_ord);
    }

    if (0)
      {
        // stk_mesh no longer allows elements to change parts if topology is changed.
        std::vector<stk::mesh::Part*> add_parts(1, &block_hex_20);
        std::vector<stk::mesh::Part*> remove_parts(1, block_hex_8);
        eMesh.get_bulk_data()->change_entity_parts( element_local, add_parts, remove_parts );
      }

    stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
    eMesh.get_bulk_data()->modification_end();
    eMesh.print_info("After quadratic",2);

    eMesh.save_as("./cube1x1x2_hex-20.e");
    //exit(1);
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

#ifdef STK_BUILT_IN_SIERRA

TEST(nodeRegistry, test_serial_hex8_tet4_24_1)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

  // start_demo_nodeRegistry_test_serial_hex8_tet4_24_1

  percept::PerceptMesh eMesh(3u);

  unsigned p_size = eMesh.get_parallel_size();
  unsigned p_rank = eMesh.get_rank();
  Util::setRank(eMesh.get_rank());

  std::string gmesh_spec = std::string("1x1x")+toString(p_size)+std::string("|bbox:0,0,0,1,1,1");
  eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
  if (1)
    {
      stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part_with_topology("refine_new_nodes_part", stk::topology::NODE);
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      meta.declare_attribute_no_delete(part, &percept::auto_part);
    }

  // prepare for breaking into tet elements
  stk::mesh::Part& block_tet_4 = eMesh.get_fem_meta_data()->declare_part("block_tet_4", stk::topology::ELEMENT_RANK);
  /// set cell topology for the part block_tet_4
  stk::mesh::set_cell_topology< shards::Tetrahedron<4>  >( block_tet_4 );
  stk::io::put_io_part_attribute(block_tet_4);

  eMesh.commit();
  eMesh.print_info();
  eMesh.save_as(std::string("./")+std::string("cube1x1x")+toString(p_size)+std::string("-orig.e"));

  //mesh::Part* block_hex_8 = const_cast<mesh::Part *>(eMesh.getPart("block_1"));

  NodeRegistry nodeRegistry(eMesh);
  nodeRegistry.initialize();

  if (p_size <= 2)
  {
    // pick an element on the processor boundary
    unsigned elem_num_local = 1;
    unsigned elem_num_ghost = 2;
    if (p_size == 1)
      elem_num_ghost = 1;

    stk::mesh::Entity element_local_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local);
    stk::mesh::Entity element_ghost_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost);
    if (p_rank == 1)
    {
      element_local_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_ghost);
      element_ghost_p = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, elem_num_local);
    }

    dw() << "P["<<p_rank<<"] elem_num_local = " << elem_num_local << DWENDL;
    dw() << "P["<<p_rank<<"] elem_num_ghost = " << elem_num_ghost << DWENDL;

    stk::mesh::Entity element_local = element_local_p;
    stk::mesh::Entity element_ghost = element_ghost_p;

    std::cout << "P["<<p_rank<<"] element_local = " << element_local << std::endl;
    std::cout << "P["<<p_rank<<"] element_ghost = " << element_ghost << std::endl;

    // choose faces to be used for new node locations

    std::vector<NeededEntityType> needed_entity_ranks(2);
    needed_entity_ranks[0] = NeededEntityType(eMesh.face_rank(), 1u);
    needed_entity_ranks[1] = NeededEntityType(stk::topology::ELEMENT_RANK, 1u);

    /*
     * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
     *   from remote proc if necessary; finally, the local node database is ready to be queried
     *
     * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
     * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
     * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
     */
    nodeRegistry.beginRegistration();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_local, needed_entity_ranks,0);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element_ghost, needed_entity_ranks,0);
    nodeRegistry.endRegistration();

    std::cout << "P["<<p_rank<<"] nodeRegistry size  = " << nodeRegistry.total_size() << std::endl;
    std::cout << "P["<<p_rank<<"] nodeRegistry lsize = " << nodeRegistry.local_size() << std::endl;

    dw() << "P["<<p_rank<<"] nodeRegistry size       = " << nodeRegistry.total_size() << DWENDL;
    dw() << "P["<<p_rank<<"] nodeRegistry lsize      = " << nodeRegistry.local_size() << DWENDL;

    // check if the newly requested nodes are local or remote
    nodeRegistry.beginCheckForRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_local, needed_entity_ranks,0);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element_ghost, needed_entity_ranks,0);
    nodeRegistry.endCheckForRemote();

    // get the new nodes from other procs if they are nonlocal
    nodeRegistry.beginGetFromRemote();
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_local, needed_entity_ranks,0);
    nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element_ghost, needed_entity_ranks,0);
    nodeRegistry.endGetFromRemote();

    // now we can get the new node's id and entity
    unsigned iSubDimOrd = 5u;  // 6'th face of element 1 is shared with 5th of element 2 (top, bot resp.)
    if (p_rank)
    {
      iSubDimOrd = 4u;
    }
    NodeIdsOnSubDimEntityType& nodeIds_onSE_0 = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_ranks[0].first, iSubDimOrd));
    stk::mesh::Entity node_0   = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, eMesh.identifier(nodeIds_onSE_0[0]));

    // should be the same node on each proc
    std::cout << "P[" << p_rank << "] nodeId_0 = " << nodeIds_onSE_0 << " node_0= " << node_0 << std::endl;
    if (p_size == 2)
    {
      EXPECT_EQ(eMesh.identifier(nodeIds_onSE_0[0]), 15u);
    }

    // end_demo

    // start_demo_nodeRegistry_test_serial_hex8_tet4_24_1_a

    // check for center node
    if (p_rank == 0)
    {
      NodeIdsOnSubDimEntityType& nodeIds_onSE_1 = *( nodeRegistry.getNewNodesOnSubDimEntity(element_local, needed_entity_ranks[1].first, 0u));
      stk::mesh::Entity node_1   = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, eMesh.identifier(nodeIds_onSE_1[0]));
      std::cout << "P[" << p_rank << "] nodeId_1 = " << nodeIds_onSE_1 << " node_1= " << node_1 << std::endl;
      if (p_size == 2)
      {
        EXPECT_EQ(eMesh.identifier(nodeIds_onSE_1[0]), 13u);
      }
    }
    //exit(1);
    // end_demo
  }
}
#endif

} // namespace unit_tests
} // namespace percept

