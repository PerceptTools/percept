// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_mesh/base/FieldParallel.hpp>

#include <gtest/gtest.h>

#include <percept/PerceptMesh.hpp>

#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept
  {
    namespace regression_tests
    {


#define DO_MEMORY_ACCOUNTING

typedef size_t MemorySizeType;

struct MemoryInfo
{
  MemorySizeType m_rss_current;
  MemorySizeType m_rss_high_water_mark;
  static const MemorySizeType MB = 1024*1024;

  MemoryInfo() { get_memory_usage(); }

  void get_memory_usage() 
  {
    m_rss_current = 0;
    m_rss_high_water_mark = 0;
    stk::get_memory_usage(m_rss_current, m_rss_high_water_mark);
  }
  void set_state() { get_memory_usage(); }
  void get_increment() {
    MemoryInfo old_state = *this;
    get_memory_usage();
    m_rss_current -= old_state.m_rss_current;
  }

};

inline double MegaByte(MemorySizeType x) { return  ((double)x/1024.0/1024.0); }

std::ostream& operator<<(std::ostream& os, const MemoryInfo& mem)
{
  char buf[1024];
  sprintf(buf, "\n%20s %20s\n%20g %20g\n", "current_rss [MB]", "max_rss [MB]",
          MegaByte(mem.m_rss_current), MegaByte(mem.m_rss_high_water_mark) );
  os << buf;
  return os;
}


TEST(adapt, count_memory)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD ;

  const unsigned p_size = stk::parallel_machine_size( pm );
  if (p_size == 1)
    {
      const unsigned n = 20;
      //const unsigned nx = n , ny = n , nz = p_size*n ;
      const unsigned nx = n , ny = n;

      percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, false);
      stk::mesh::Part& part = fixture.meta_data.declare_part_with_topology("elems", stk::topology::TRI_3);

      fixture.meta_data.commit();
      fixture.generate_mesh();

      percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
      //eMesh.print_info("quad mesh",2);

      //const size_t num_new_tris = 2000*2000;
      const size_t num_new_tris = 20*20;
      const size_t num_nodes_per_tri = 3;
      const size_t num_new_nodes = num_new_tris*num_nodes_per_tri;
      MemoryInfo mem_delta_node;
      double time = -percept::Util::cpu_time();

      std::vector<stk::mesh::Entity> new_nodes, new_elements;

      eMesh.get_bulk_data()->modification_begin();
      eMesh.createEntities(eMesh.node_rank(), num_new_nodes, new_nodes);
      stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
      eMesh.get_bulk_data()->modification_end();

      mem_delta_node.get_increment();
      double mem_per_node = double(mem_delta_node.m_rss_current)/double(num_new_nodes);
      std::cout << "\nstk_mesh count_memory mem_per_node = " << mem_per_node << "\n" << std::endl;

      MemoryInfo mem_delta_elem_0, mem_delta_elem_1;

      eMesh.get_bulk_data()->modification_begin();
      eMesh.createEntities(stk::topology::ELEMENT_RANK, num_new_tris, new_elements);
      //stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
      //eMesh.get_bulk_data()->modification_end();

      mem_delta_elem_0.get_increment();

      //eMesh.get_bulk_data()->modification_begin();
      size_t i_node=0;
      for (size_t i=0; i<num_new_tris; ++i) {

        std::vector<stk::mesh::Part*> add_parts(1,&part);
        std::vector<stk::mesh::Part*> remove_parts;
        eMesh.get_bulk_data()->change_entity_parts( new_elements[i], add_parts, remove_parts );

        unsigned ordinal = 0;
        for (size_t j=0; j < num_nodes_per_tri; ++j) {
          eMesh.get_bulk_data()->declare_relation(new_elements[i], new_nodes[i_node],ordinal);
          ++ordinal;
          ++i_node;
        }
      }
      stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
      eMesh.get_bulk_data()->modification_end();

      mem_delta_elem_1.get_increment();
      double mem_per_elem_0 = double(mem_delta_elem_0.m_rss_current)/double(num_new_tris);
      double mem_per_elem_1 = double(mem_delta_elem_1.m_rss_current)/double(num_new_tris);

      time += percept::Util::cpu_time();

      std::cout << "\nstk_mesh count_memory mem_per_elem (no connectivity) = " << mem_per_elem_0 << " with connectivity= " << mem_per_elem_1 << " cpu= " << time << std::endl;

    }
}


    }
  }

