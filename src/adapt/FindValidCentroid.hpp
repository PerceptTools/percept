// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef FindValidCentroid_hpp
#define FindValidCentroid_hpp


#include <percept/PerceptMesh.hpp>
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>

namespace percept {

  class FindValidCentroid {
  public:
    PerceptMesh& m_eMesh;
    int ndiv;
    bool m_debug;
    bool m_use_finite_volume;
    FindValidCentroid(PerceptMesh& eM, int nd=5, bool deb=false, bool use_finite_volume=false) : m_eMesh(eM), ndiv(nd), m_debug(deb), m_use_finite_volume(use_finite_volume) {
      //std::cout << "m_debug= " << m_debug << std::endl;
      if (eM.getProperty("FindValidCentroid_use_finite_volume") == "true")
        m_use_finite_volume = true;
    }

    double getVolumes(std::vector<double>& volumes, stk::mesh::Entity element);
    double metric(std::vector<double>& volumes, bool& foundBad);
    // return if changed
    bool findCentroid(stk::mesh::Entity element, double *c_p, std::vector<stk::mesh::Entity>& nodes, stk::mesh::Entity c_node);
  };

}
#endif
