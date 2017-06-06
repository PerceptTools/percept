// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

int main(int argc, char **argv)
{
  Kokkos::initialize(argc,argv);

  testing::InitGoogleTest(&argc, argv);

  int returnVal = RUN_ALL_TESTS();

  Kokkos::finalize();

  return returnVal;
}
