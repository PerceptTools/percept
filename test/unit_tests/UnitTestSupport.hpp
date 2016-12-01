// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_unit_tests_UnitTestSupport_hpp
#define adapt_unit_tests_UnitTestSupport_hpp

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>

#include <gtest/gtest.h>

#include <stk_io/IossBridge.hpp>

#include <boost/lexical_cast.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <percept/Percept.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <percept/RunEnvironment.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#include <stk_util/parallel/Parallel.hpp>

  namespace percept {
    namespace unit_tests {


      class UnitTestSupport
      {
      public:
        /// The following defines where to put the input and output files created by this set of functions
        static const std::string input_files_loc;
        static const std::string output_files_loc;

        static bool always_do_regression_tests;

        /// This function either writes the given mesh to a file in Exodus format (option 0)
        ///   or, under option 1, checks if the file already exists, and if so, treats that
        ///   file as the "gold" copy and does a regression difference check.
        /// Overriden by always_do_regression_tests - if "true", then option 1 is always done.

        static void save_or_diff(PerceptMesh& eMesh, std::string filename, int option = 0);



      };
   }
  }

#endif
