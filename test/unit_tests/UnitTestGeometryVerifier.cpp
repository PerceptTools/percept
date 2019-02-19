// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/function/FunctionOperator.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/TopologyVerifier.hpp>
#include <percept/GeometryVerifier.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <gtest/gtest.h>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/fixtures/Fixture.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

namespace percept {
namespace unit_tests {

// on pathscale platform this doesn't work (something to do with static variables)

static int dw_enabled = 1;
static stk::diag::Writer s_diagWriter(std::cout.rdbuf(), dw_enabled);
static stk::diag::Writer &
dw()
{
  //static stk::diag::Writer s_diagWriter(dwout().rdbuf(), 0);

  s_diagWriter.setPrintMask(LOG_NORM+LOG_ALWAYS);

  return s_diagWriter;
}

#define EXTRA_PRINT 0

//=============================================================================
//=============================================================================
//=============================================================================

TEST(geom, volume)
{
  dw().m(LOG_GEOMETRY_VERIFIER) << "TEST::geom::volume " << stk::diag::dendl;

  const size_t num_x = 3;
  const size_t num_y = 3;
  const size_t num_z = 3;
  std::string config_mesh =
    Ioss::Utils::to_string(num_x) + "x" +
    Ioss::Utils::to_string(num_y) + "x" +
    Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
  PerceptMesh eMesh(3u);
  eMesh.new_mesh(GMeshSpec(config_mesh));
  eMesh.commit();
  // no need for this in create mode: eMesh.readBulkData();

  //MetaData& metaData = *eMesh.get_fem_meta_data();
  stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();

  eMesh.dump();
  GeometryVerifier geomVerifier(false);
  geomVerifier.isGeometryBad(bulkData, true);
  //setDoPause(true);
  //pause();
}

}
}
