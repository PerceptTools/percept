// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_PerceptUtils_hpp
#define percept_PerceptUtils_hpp

#include <string>
#include <mpi.h>

namespace stk {
struct topology;
namespace mesh {
  class Entity;
  class FieldBase;
}
namespace diag {
  class Timer;
}
}

class CellTopologyData;

namespace percept {

  void computeCentroid(stk::mesh::Entity entity, double centroid[3], const stk::mesh::FieldBase & coord_field);

  double volume(stk::mesh::Entity element, const stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data);

  stk::diag::Timer& rootTimerStructured();

  void printTimersTableStructured();

  void get_memory_high_water_mark_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg, size_t& hwm_sum);
  void get_memory_now_across_processors(MPI_Comm comm, size_t& hwm_max, size_t& hwm_min, size_t& hwm_avg, size_t& hwm_sum);
  std::string print_memory_high_water_mark(MPI_Comm comm);
  std::string print_memory_now(MPI_Comm comm);
  std::string print_memory_both(MPI_Comm comm);

  void convert_stk_topology_to_ioss_name(
          const stk::topology stk_topo,
          std::string& ioss_topo);
}

#endif
