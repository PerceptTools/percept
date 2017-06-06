// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_PerceptUtils_hpp
#define percept_PerceptUtils_hpp

namespace stk {
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
}

#endif
