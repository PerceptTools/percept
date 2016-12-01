// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_RefinementInfoByType_hpp
#define adapt_RefinementInfoByType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <stdint.h>
#include <limits>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>



  namespace percept {

    struct RefinementInfoByType
    {
      RefinementInfoByType() : 
        m_numOrigElems(0),
        m_numNewElems(0),
        m_numOrigNodes(0),
        m_numNewNodes(0),

        m_numOrigElemsLast(0),
        m_numNewElemsLast(0),
        m_rank(0)
      {}

      typedef uint64_t RefinementInfoCount ;
      //typedef unsigned RefinementInfoCount ;

      RefinementInfoCount m_numOrigElems;
      RefinementInfoCount m_numNewElems;
      shards::CellTopology m_topology;
      RefinementInfoCount m_numOrigNodes;
      RefinementInfoCount m_numNewNodes;

      RefinementInfoCount m_numOrigElemsLast;
      RefinementInfoCount m_numNewElemsLast;
      unsigned m_rank;

      static void estimateNew(std::vector< RefinementInfoByType >& refinementInfoByType, int iRefinePass);
      static void printTable(std::ostream& os, std::vector< RefinementInfoByType >& refinementInfoByType, int iRefinePass, bool printAll = false);
      static void countCurrentNodes(percept::PerceptMesh& eMesh, std::vector< RefinementInfoByType >& refinementInfoByType);

    };

  }

#endif
