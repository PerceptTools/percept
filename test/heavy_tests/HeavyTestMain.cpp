// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>

#if defined( STK_HAS_MPI )
#include <mpi.h>eavy_tests/HeavyTe
#endif
#if defined(WITH_KOKKOS)
#include <Kokkos_Core.hpp>
#endif
#include <percept/fixtures/Fixture.hpp>
#include <percept/RunEnvironment.hpp>
#include <percept/PerceptUtils.hpp>

#include <percept/pyencore.h>

#if !PY_PERCEPT 

#include <gtest/gtest.h>
#include <mpi.h>
#include <Kokkos_Core.hpp>

int gl_argc=0;
char** gl_argv=0;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);


#if defined(WITH_KOKKOS)
    Kokkos::initialize(argc, argv);
	std::cout << "INFO: Running tests WITH_KOKKOS set to true. Recompile without Kokkos if you don't want this." << std::endl;

#endif

    testing::InitGoogleTest(&argc, argv);

    gl_argc = argc;
    gl_argv = argv;

    int returnVal = RUN_ALL_TESTS();


#if defined(WITH_KOKKOS)


    Kokkos::finalize();
#endif

    MPI_Finalize();
     
    percept::printTimersTableStructured();

    return returnVal;
}

#else
  int main() {return 0;}
#endif
