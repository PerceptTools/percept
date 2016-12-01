// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinerUtil_hpp
#define adapt_RefinerUtil_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <utility>
#include <math.h>
#include <map>
#include <set>
#include <vector>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <percept/stk_mesh.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/ProgressMeter.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>
#include <adapt/NodeRegistry.hpp>

#include <adapt/SubDimCell.hpp>

#include <adapt/RefinementInfoByType.hpp>


namespace percept {


class RefinerUtil
{
public:

  static BlockNamesType
  getBlockNames(const std::string& block_name, unsigned proc_rank, percept::PerceptMesh& eMesh, const std::string& geomFile="");

  static BlockNamesType
  correctBlockNamesForPartPartConsistency(percept::PerceptMesh& eMesh, BlockNamesType& blocks);

  static BlockNamesType
  correctBlockNamesForPartPartConsistency_1(percept::PerceptMesh& eMesh, BlockNamesType& blocks, const std::string& geomFile);

  /// For use with the Encore/Percept interface.
  /// for each element in elements_to_unref, add parents to the list, repeatedly,
  ///   for num_levels_to_add times, thus adding grand-parents, great-grandparents, ...
  /// If num_levels_to_add is < 0, add all up to root (ie, infinite number of levels)
  /// Throw an error if the parents are already in the list.

  static void
  addAncestorsToUnrefineList(percept::PerceptMesh& eMesh, int num_levels_to_add, ElementUnrefineCollection& elements_to_unref);

  /// create missing edges after adapt - for edge-based simulators

  static void
  create_missing_edges(percept::PerceptMesh& eMesh);

  static void
  rebuild_family_tree(PerceptMesh& eMesh, bool debug=false);

  /// NodeRegistry
  static void rebuild_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, bool initNR = true, PerceptMesh *eMeshNR = 0, NodeRegistry *compareNR=0, bool skipEmpty = true);
  static void save_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, const std::string& msg, bool doComm=true);
  static void compare(PerceptMesh& eMesh0, NodeRegistry& nr0, PerceptMesh& eMesh1, NodeRegistry& nr1, bool skipEmpty = true);
};

}

#endif
