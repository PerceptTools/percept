// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef BasisTable_hpp
#define BasisTable_hpp

#include "Teuchos_RCP.hpp"
#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>
#include <percept/function/MDArray.hpp>


namespace Intrepid {
	template<class Scalar, class ArrayScalar>
	class Basis;
}

  namespace percept {
    class BasisTable
    {
    public:
      typedef Intrepid::Basis<double, MDArray > BasisType;
      typedef Teuchos::RCP<BasisType>           BasisTypeRCP;
      typedef std::map<unsigned, BasisTypeRCP > BasisTableMap;
      static BasisTypeRCP getBasis(shards::CellTopology& topo);
      static void setupBasisTable();

    private:
      static BasisTableMap m_basisTable;


    };

  }

#endif
