// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <stk_util/environment/memory_util.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/MemoryTracking.cpp>
#include <stk_mesh/base/MemoryUsage.hpp>
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>   // for Entity
#include <stk_io/FillMesh.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <climits>

namespace percept {
namespace unit_tests {

TEST(timer, verify_timer)
{
  //if (1<stk::parallel_machine_size(MPI_COMM_WORLD)) return;

  stk::diag::TimerSet timerSet(sierra::Diag::TIMER_ALL);
  stk::diag::Timer timer(stk::diag::createRootTimer("MeshAdapt", timerSet));

  timer.start();

  usleep(1e4);

  {
    stk::diag::Timer  timer_f1("Feature1", timer);
    stk::diag::TimeBlock tb_f1(timer_f1);

    usleep(2e4);

    for (int i=0; i<3; i++)
    {
      stk::diag::Timer  timer_f2("Feature2", timer_f1);
      stk::diag::TimeBlock tb_f2(timer_f2);
      
      usleep(1e4);
      
      for (int j=0; j<2; j++)
      {
        stk::diag::Timer  timer_f3("Feature3", timer_f2);
        stk::diag::TimeBlock tb_f3(timer_f3);
        
        usleep(5e4);
      }
    }
  }

  {
    stk::diag::Timer  timer_f4("Feature4", timer);
    stk::diag::TimeBlock tb_f4(timer_f4);

    usleep(4e4);
  }

  {
    std::ostringstream str;
    timer.stop();
    stk::diag::printTimersTable(str, timer,
                                stk::diag::METRICS_ALL, false,
                                MPI_COMM_WORLD);

    if (0==stk::parallel_machine_rank(MPI_COMM_WORLD))
      std::cout << str.str() << std::endl;
  }
}

}
}
