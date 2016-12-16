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

#include <mpi.h>

#include <gtest/gtest.h>

int gl_argc=0;
char** gl_argv=0;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    testing::InitGoogleTest(&argc, argv);

    gl_argc = argc;
    gl_argv = argv;

    int returnVal = RUN_ALL_TESTS();

    MPI_Finalize();

    return returnVal;
}
