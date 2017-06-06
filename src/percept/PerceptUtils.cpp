// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/PerceptUtils.hpp>
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/diag/Timer.hpp>

#include <Shards_CellTopologyData.h>

namespace percept {

void computeCentroid(stk::mesh::Entity elem, double centroid[3], const stk::mesh::FieldBase & coord_field)
{
  int spaceDim = coord_field.mesh_meta_data().spatial_dimension();
  
  const stk::mesh::BulkData& mesh = coord_field.get_mesh();
  
  const stk::mesh::Entity * nodes = mesh.begin(elem, stk::topology::NODE_RANK);
  const unsigned        num_nodes = mesh.num_connectivity(elem, stk::topology::NODE_RANK);

  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  double nnode = num_nodes;

  for (unsigned ii=0; ii < num_nodes; ++ii)
    {
      double * node_coord_data_0 = (double*)stk::mesh::field_data(coord_field, nodes[ii]);
      for (int iSpaceDimOrd = 0; iSpaceDimOrd < spaceDim; iSpaceDimOrd++)
        {
          centroid[iSpaceDimOrd] += node_coord_data_0[iSpaceDimOrd] / nnode;
        }
    }
}

double volume(stk::mesh::Entity element, const stk::mesh::FieldBase *coord_field, const CellTopologyData * cell_topo_data)
{
  VolumeUtil jacA;
  shards::CellTopology cell_topo(cell_topo_data);
  double volScale = jacA.getJacobianToVolumeScale(cell_topo);
  
  double jacobian = 0.0;
  jacA(jacobian, element, coord_field, cell_topo_data);
  double cellVol = jacobian*volScale;
  return cellVol;
}

stk::diag::Timer& rootTimerStructured() { 
  static stk::diag::TimerSet s_timerSet(sierra::Diag::TIMER_ALL);
  static stk::diag::Timer s_timer = stk::diag::createRootTimer("Structured", s_timerSet);
  return s_timer;
}

void printTimersTableStructured() {
  std::ostringstream str;
  rootTimerStructured().stop();
  stk::diag::printTimersTable(str, rootTimerStructured(), stk::diag::METRICS_ALL, false);
  
  std::cout << str.str() << std::endl;
}
}
