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
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <gtest/gtest.h>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <gtest/gtest.h>
#include <boost/lexical_cast.hpp>
#include <stk_io/IossBridge.hpp>

#include <percept/Percept.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>
#include <unit_tests/TestLocalRefiner.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <percept/RunEnvironment.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/PyramidFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

// smoothing tests
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/ReferenceMeshSmoother1.hpp>


// this is for testing the local-refine refactoring
#define UNIFORM_REFINER UniformRefiner
//#define UNIFORM_REFINER TestLocalRefiner


  namespace percept
  {
    namespace heavy_tests
    {
      //using namespace regression_tests;

#if 1
      static const std::string path_sep = "._.";
      static const std::string input_files_loc="./input_files"+path_sep;
      static const std::string output_files_loc="./output_files"+path_sep;
#else
      static const std::string input_files_loc="./input_files/";
      static const std::string output_files_loc="./output_files/";
#endif

#define EXTRA_PRINT 0

      static std::string procs_string[9] = {"np0", "np1", "np2", "np3", "np4", "np5", "np6", "np7", "np8"};


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      //========= AREA for tests in progress of being debugged
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================



    }//    namespace heavy_tests
  }//  namespace percept


