// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
#include <percept/Percept.hpp>

#include "MeshDifference.hpp"


using namespace percept;

int main(int argc,  char **argv)
{
  
  MeshDifference md;
  md.run(argc, argv);


  return 0;
}
