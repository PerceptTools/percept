// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.




#include <stdexcept>

#include <percept/stk_rebalance/Partition.hpp>

using namespace stk;
namespace stk {
  using namespace rebalance;
}

Partition::Partition(ParallelMachine comm) : comm_(comm) { }
Partition::~Partition() { }



