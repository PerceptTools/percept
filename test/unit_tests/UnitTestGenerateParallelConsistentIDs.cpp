// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stddef.h>                     // for size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/parallel/GenerateParallelConsistentIDs.hpp>
#include <gtest/gtest.h>
#include <string>                       // for string
#include <stk_util/util/ReportHandler.hpp>
#include <stdint.h>

#include <percept/PerceptMesh.hpp>

uint64_t my_hash(uint64_t x1, uint64_t x2) {
  //return x1*x2+x1+x2;
  return 100*x1+x2;
}

TEST(DISABLED_UnitTestParallelConsistentIDs, createEdgeNodeIDs) {
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  uint64_t maxAllowableId = ~0U;

  std::vector<uint64_t> newIds;
  std::vector<uint64_t> localOrderArray;

  {
    // simple mesh with node IDs for 
    // original vertex nodes (1-6) and new edge nodes (7-15)
    //
    // 3 -11-2 -12-6
    // | \   | \   |
    // 13 8  7 15 14
    // |   \ |   \ |
    // 4 -9- 1 -10-5

    // generate new node IDs for edges

    std::vector<uint64_t> existingIds;
    switch (mpi_size) {
    case 4:
      {
        if      (mpi_rank==0) {
          existingIds = {1,3,4};
          localOrderArray = {my_hash(1,3),my_hash(1,4),my_hash(3,4)};
        }
        else if (mpi_rank==1) {
          existingIds = {1,2,3};
          localOrderArray = {my_hash(2,3),my_hash(1,2)};
        }
        else if (mpi_rank==2) {
          existingIds = {1,5,2};
          localOrderArray = {my_hash(1,5),my_hash(2,5)};
        }
        else if (mpi_rank==3) {
          existingIds = {5,6,2};
          localOrderArray = {my_hash(5,6),my_hash(2,6)};
        }
      }
      break;    
    case 3:
      {
        if      (mpi_rank==0) {
          existingIds = {1,3,4};
          localOrderArray = {my_hash(1,3),my_hash(1,4),my_hash(3,4)};
        }
        else if (mpi_rank==1) {
          existingIds = {1,2,3,5}; // take 2 inner elements
          localOrderArray = {my_hash(1,2),my_hash(2,3),my_hash(1,5),my_hash(2,5)};
        }
        else if (mpi_rank==2) {
          existingIds = {5,6,2};
          localOrderArray = {my_hash(5,6),my_hash(2,6)};
        }
      }
      break;    
    case 2:
      {
        if (mpi_rank==0) {
          existingIds = {1,2,3,4};
          localOrderArray = {my_hash(1,2),my_hash(1,3),my_hash(1,4),my_hash(2,3),my_hash(3,4)};
        }
        else {
          existingIds = {1,2,5,6};
          localOrderArray = {my_hash(1,5),my_hash(2,5),my_hash(2,6),my_hash(5,6)};
        }
      }
      break;    
    case 1:
      {
        existingIds = {1,2,3,4,5,6};
        localOrderArray = {my_hash(1,2),my_hash(1,3),my_hash(1,4),my_hash(1,5),
                           my_hash(2,3),my_hash(2,5),my_hash(2,6),
                           my_hash(3,4),
                           my_hash(5,6)};
      }
      break;
    }

    newIds = stk::generate_parallel_consistent_ids(maxAllowableId, existingIds, localOrderArray, MPI_COMM_WORLD);

    uint64_t localSize = newIds.size();
    uint64_t globalSize;
    MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    EXPECT_EQ(globalSize, (uint64_t)9);

    // DEBUG
    {
      std::vector<uint64_t> newIdsGlobal;
      stk::parallel_vector_concat(MPI_COMM_WORLD, newIds, newIdsGlobal);
      std::sort(newIdsGlobal.begin(), newIdsGlobal.end());
      //
      //  Verify properties of the returned ids
      //
      for(uint64_t i=0; i<newIdsGlobal.size(); ++i) {
        uint64_t curNewID = newIdsGlobal[i];
        
        EXPECT_EQ(curNewID, (uint64_t)(7+i));
      }
    }
  }
}

TEST(DISABLED_UnitTestParallelConsistentIDs, createEdgeNodeIDsSimpleMesh) {
  int mpi_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  //int mpi_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  uint64_t maxAllowableId = ~0U;

  std::vector<uint64_t> newIds;
  std::vector<uint64_t> localOrderArray;

  percept::PerceptMesh eMesh(3u);
  eMesh.set_ioss_read_options("auto-decomp:yes");
  eMesh.open("four_tri3.g");
  eMesh.commit();
  
  std::vector<uint64_t> existingIds;
  // populate existingIds from local nodes
  {
    // get nodes
    std::vector<stk::mesh::Entity> nodes;
    eMesh.get_bulk_data()->get_entities(stk::topology::NODE_RANK, 
                                        eMesh.get_fem_meta_data()->locally_owned_part(), nodes);

    std::vector<uint64_t> nodeIds;
    for (auto node : nodes) {

      if (eMesh.owned(node)) {
        //std::cout << "P" << mpi_rank << ": owns node = " << eMesh.identifier(node) << std::endl;
      }

      nodeIds.push_back(eMesh.identifier(node));
    }

    existingIds.assign(nodeIds.begin(), nodeIds.end());
  }

  //std::cout << "P" << mpi_rank << ": num node IDs = " << existingIds.size() << std::endl;

  // populate newIds from edges in parallel
  {
    // get elems
    std::vector<stk::mesh::Entity> elems;          
    eMesh.get_bulk_data()->get_entities(stk::topology::ELEMENT_RANK, 
                                        eMesh.get_fem_meta_data()->locally_owned_part(), elems);

    //std::cout << "P" << mpi_rank << ": num elems = " << elems.size() << std::endl;
  
    // loop elems
    std::vector<stk::mesh::Entity> side_node_entities;
    std::vector<uint64_t> side_node_Ids;
    std::set<uint64_t> localOrderSet;     
    for (auto elem : elems) {
      // loop edges - HACK only for tri3 now
      for (int iEdge=0; iEdge<3; iEdge++) {
        
        eMesh.element_side_nodes(elem, iEdge, stk::topology::EDGE_RANK, side_node_entities);
        
        EXPECT_EQ(side_node_entities.size(), (unsigned)2);

        side_node_Ids = {eMesh.identifier(side_node_entities[0]), 
                         eMesh.identifier(side_node_entities[1])};
        std::sort(side_node_Ids.begin(), side_node_Ids.end());
        
        stk::mesh::Entity neigh = eMesh.get_face_neighbor(elem, iEdge);

        // if we own the edge = own the node with lowest ID
        // or if there is no edge neighbor
        // add to localOrderArray
        /*if (eMesh.owned(eMesh.get_entity(stk::topology::NODE_RANK,side_node_Ids[0])) ||
          !eMesh.is_valid(neigh)) {*/
        if (!eMesh.is_valid(neigh) || eMesh.owner_rank(neigh) >= mpi_rank) {
          //std::cout << "P" << mpi_rank << ":   adding edge = " 
          //        << side_node_Ids[0] << " "
          //        << side_node_Ids[1] << std::endl;
          
          localOrderSet.insert(my_hash(side_node_Ids[0],side_node_Ids[1]));
        }
        else {
          //std::cout << "P" << mpi_rank << ":   skipping edge = " 
          //        << side_node_Ids[0] << " "
          //        << side_node_Ids[1] << std::endl;
          
        }
      }
    }
    localOrderArray.assign(localOrderSet.begin(),localOrderSet.end());
  }

  //std::cout << "P" << mpi_rank << ": num edges added = " << localOrderArray.size() << std::endl;
    
  newIds = stk::generate_parallel_consistent_ids(maxAllowableId, existingIds, localOrderArray, MPI_COMM_WORLD);
  
  uint64_t localSize = newIds.size();
  uint64_t globalSize;
  MPI_Allreduce(&localSize, &globalSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  EXPECT_EQ(globalSize, (uint64_t)9);
  
  {
    std::vector<uint64_t> newIdsGlobal;
    stk::parallel_vector_concat(MPI_COMM_WORLD, newIds, newIdsGlobal);
    std::sort(newIdsGlobal.begin(), newIdsGlobal.end());
    //
    //  Verify properties of the returned ids
    //
    for(uint64_t i=0; i<newIdsGlobal.size(); ++i) {
      uint64_t curNewID = newIdsGlobal[i];
      
      if (mpi_rank==0)
        std::cout << curNewID << std::endl;

      EXPECT_EQ(curNewID, (uint64_t)(7+i)); //
    }
  }
}
