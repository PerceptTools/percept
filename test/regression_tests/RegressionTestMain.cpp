// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <cstdlib>
#include <cstring>
#include <iostream>
#include <utility>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <percept/fixtures/Fixture.hpp>
#include <percept/RunEnvironment.hpp>
#include <percept/PerceptUtils.hpp>

#include <percept/pyencore.h>

#if !PY_PERCEPT 

#include <gtest/gtest.h>
#include <mpi.h>
#include <Kokkos_Core.hpp>

int main(int argc, char **argv)
{
  testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);

  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank != 0) {
    delete listeners.Release(listeners.default_result_printer());
  }

  int returnVal = RUN_ALL_TESTS();
  
  percept::printTimersTableStructured();
  
  Kokkos::finalize();
  MPI_Finalize();
  
  return returnVal;
}

#else
  int main() {return 0;}
#endif
