// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/MemoryTracking.cpp>
#include <stk_mesh/base/MemoryUsage.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>   // for Entity
#include <stk_io/FillMesh.hpp>

#include <percept/PerceptMesh.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <climits>

namespace percept {
namespace unit_tests {

void my_get_memory_usage(size_t & now)
{
  std::string memtotal;
  std::string memfree;
  std::string line(128,'\0');

  size_t total = 0;
  size_t free = 0;

  std::ifstream proc_status("/proc/meminfo");
  if (!proc_status) return;

  while (memtotal.empty() || memfree.empty())
  {

    if(!std::getline(proc_status, line)) return;

    /* Find VmHWM */
    else if (line.substr(0, 9) == "MemTotal:")
    {
        memtotal = line.substr(10);
        std::istringstream iss(memtotal);
        iss >> total;
        total *= 1024;
    }
    /* Find VmRSS */
    else if (line.substr(0, 8) == "MemFree:")
    {
        memfree = line.substr(9);
        std::istringstream iss(memfree);
        iss >> free;
        free *= 1024;
    }
  }
  proc_status.close();

  now = total - free;

  //std::cout << "my_get_memory_usage (total/free/now): " << total << " " << free << " " << now << std::endl;
}

void print_memory(stk::ParallelMachine comm, size_t &current, size_t &hwm)
{
    stk::get_memory_usage(current, hwm);

    stk::all_reduce( comm, stk::ReduceSum<1>( &current ) );
    stk::all_reduce( comm, stk::ReduceSum<1>( &hwm ) );

    const size_t MB = 1024*1024;

    static int counter = 1;
    if (0==stk::parallel_machine_rank(comm))
        std::cerr << "Current(" << counter << ") : " << (double)current/MB
            << "\thwm: " << (double)hwm/MB << std::endl;
    ++counter;
}

TEST(memory, verify_get_memory_usage)
{
    if (1<stk::parallel_machine_size(MPI_COMM_WORLD)) return;

    size_t new0, hwm0;

    stk::get_memory_usage(new0, hwm0);

    size_t new1;
    my_get_memory_usage(new1);

    std::cout << std::setw(3) << "i"
            << std::setw(18) << "real mem increase"
            << std::setw(18) << "stk: mem increase"
            << std::setw(18) << "exp: mem increase"
            << std::endl;

    size_t array_size = 1;
    for (unsigned i=1; i<=9; i++) {
        double * x = new double[array_size];

        // if you don't touch the memory it does not seem to get allocated
        for (unsigned j=0; j<array_size; j++)
            x[j]=j*2-1;

        size_t new2, hwm2;
        stk::get_memory_usage(new2, hwm2);

        size_t new3;
        my_get_memory_usage(new3);

        std::cout << std::setw(3) << i
                << std::setw(18) << sizeof(double)*array_size
                << std::setw(18) << new2-new0
                << std::setw(18) << new3-new1
                << std::endl;

        delete [] x;
        array_size *= 10;
    }
}

TEST(memory, serial_stk_from_exodus)
{
    if (1<stk::parallel_machine_size(MPI_COMM_WORLD)) return;

    size_t new0, hwm0;
    stk::get_memory_usage(new0, hwm0);

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

    stk::io::fill_mesh("input.g", bulk);

    size_t new1, hwm1;
    stk::get_memory_usage(new1, hwm1);

    std::cout << "memory increase (bytes): current = " << (new1-new0)
              << " hwm = " << (hwm1-hwm0) << std::endl;
}

TEST(memory, vary_bucket_size)
{
    size_t new0, hwm0;
    print_memory(MPI_COMM_WORLD, new0, hwm0);

    for (unsigned i=7; i<=13; i++) {
        const unsigned bucket_size = pow(2,i);

        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA, false, NULL, bucket_size);

        stk::io::fill_mesh_with_auto_decomp("input.g", bulk);

        print_memory(MPI_COMM_WORLD, new0, hwm0);
    }
}

TEST(memory, delete_face_adjacent_element_graph)
{
    size_t new0, hwm0;
    print_memory(MPI_COMM_WORLD, new0, hwm0);

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

    stk::io::fill_mesh_with_auto_decomp("input.g", bulk);

    print_memory(MPI_COMM_WORLD, new0, hwm0);

    bulk.delete_face_adjacent_element_graph();

    print_memory(MPI_COMM_WORLD, new0, hwm0);
}

TEST(memory, Percept_generated_cube)
{
    size_t new0, hwm0;
    print_memory(MPI_COMM_WORLD, new0, hwm0);

    percept::PerceptMesh eMesh(3u);

    std::string gmesh_spec = std::string("10x10x10|bbox:0,0,0,1,1,1");
    eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
    eMesh.commit();

    print_memory(MPI_COMM_WORLD, new0, hwm0);
}

TEST(memory, generated_cube)
{
    size_t new0, hwm0;
    print_memory(MPI_COMM_WORLD, new0, hwm0);

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

    std::string gmesh_spec = std::string("10x10x10|bbox:0,0,0,1,1,1");

    stk::io::StkMeshIoBroker stkIo;

    stkIo.set_bulk_data(bulk);
    stkIo.add_mesh_database(gmesh_spec, "generated", stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
    stkIo.populate_bulk_data();

    print_memory(MPI_COMM_WORLD, new0, hwm0);
}

void createEntities(stk::mesh::BulkData& bulk,
        stk::mesh::EntityRank entityRank, int count, std::vector<stk::mesh::Entity>& requested_entities)
{
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    std::vector<size_t> requests(meta.entity_rank_count() , 0 );
    requests[entityRank] = count;
    bulk.generate_new_entities( requests, requested_entities );
    if (entityRank ==stk::topology::NODE_RANK)
    {
        stk::mesh::Part& nodePart = meta.get_topology_root_part(stk::topology::NODE);
        stk::mesh::PartVector nodeParts(1, &nodePart);
        for(size_t i=0; i < requested_entities.size(); ++i) {
            bulk.change_entity_parts(requested_entities[i], nodeParts);
        }
    }
}

void create_elements_and_nodes(stk::mesh::BulkData& bulk, size_t n_elements, size_t n_nodes)
{
    stk::mesh::MetaData& meta = bulk.mesh_meta_data();
    stk::mesh::Part& part = meta.declare_part_with_topology("elems", stk::topology::HEX_8);

    std::vector<stk::mesh::Entity> new_elements;
    std::vector<stk::mesh::Entity> new_nodes;

    meta.commit();

    bulk.modification_begin();

    createEntities(bulk, stk::topology::ELEMENT_RANK, n_elements, new_elements);

    createEntities(bulk, stk::topology::NODE_RANK, n_nodes, new_nodes);

    size_t i_node = 0;
    int n_node_per_element = 8;
    for (size_t i_element = 0; i_element < n_elements; i_element++)
    {
        stk::mesh::Entity element = new_elements[i_element];
        std::vector<stk::mesh::Part*> add_parts(1,&part);
        std::vector<stk::mesh::Part*> remove_parts;
        bulk.change_entity_parts( element, add_parts, remove_parts );

        for (int j_node = 0; j_node < n_node_per_element; j_node++)
        {
            stk::mesh::Entity node = new_nodes[i_node];

            bulk.declare_relation(element, node, j_node);

            i_node++;
            if (i_node >= n_nodes-1)
                i_node = 0;
        }
    }

    stk::mesh::fixup_ghosted_to_shared_nodes(bulk);
    bulk.modification_end();
}

TEST(memory, create_entities)
{
    size_t new0, hwm0;
    print_memory(MPI_COMM_WORLD, new0, hwm0);

    for (unsigned i=3; i<=4; i++) {

        const size_t size = pow(10,i);

        stk::mesh::MetaData meta(3);
        stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

        std::cout << "testing creation of " << size << " number of nodes and elements" << std::endl;

        create_elements_and_nodes(bulk, size, size);

        size_t new1, hwm1;
        print_memory(MPI_COMM_WORLD, new1, hwm1);
    }
}

}
}
