// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>

#include <stk_util/util/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>
#include <gtest/gtest.h>

#include <stk_search/CoarseSearch.hpp>

#include <iostream>

using namespace stk::diag;

static const int box_count = 100;

namespace percept {
namespace unit_tests {
/// PLATFORM_NOTE gcc4.3
/// If DIM here is set to 2, the code failes to compile due to dependence of Oct tree code on assuming Dim = 3 in bounding boxes
#define DIM 3

#ifndef REDS
TEST(search, test1)
{
  typedef stk::search::IdentProc<uint64_t,unsigned> IdentProc;
  typedef stk::search::Point<float> Point;
  typedef stk::search::Box<float> Box;
  typedef std::pair<Point,IdentProc> BoundingPoint;
  typedef std::pair<Box,IdentProc> BoundingBox;
  typedef std::vector<std::pair<IdentProc, IdentProc> > IdentProcRelation;

  stk::ParallelMachine  comm = MPI_COMM_WORLD;
  //stk::diag::WriterThrowSafe _write_throw_safe(dw());

  //!dw().m(LOG_SEARCH) << "Use case 1" << stk::diag::push << stk::diag::dendl;

  int parallel_rank = stk::parallel_machine_rank(comm);
  //  int parallel_size = stk::parallel_machine_size(comm);

  std::vector<BoundingBox> domain_vector;

  Point min_corner, max_corner;

  for (int i = parallel_rank*box_count; i < (parallel_rank + 1)*box_count; ++i) {
    min_corner[0] = i;
    min_corner[1] = 0;
    min_corner[2] = 0;

    max_corner[0] = i+1;
    max_corner[1] = 1;
    max_corner[2] = 1;

    BoundingBox   domain;
    domain.second.set_id(i);
    domain.first = Box(min_corner,max_corner);

    domain_vector.push_back(domain);
  }

  std::vector<BoundingPoint> range_vector;
  Point center;
  for (int i = parallel_rank*box_count; i < (parallel_rank + 1)*box_count; ++i) {
    center[0] = i+0.5f;
    center[1] = 0.5f;
    center[2] = 0.5f;

    BoundingPoint   p;
    p.second.set_id(i);
    p.first = center;

    range_vector.push_back(p);
  }

  //dw().m(LOG_SEARCH) << "range  " << range_vector << dendl;
  //dw().m(LOG_SEARCH) << "domain " << domain_vector << dendl;

  //dw().m(LOG_SEARCH) << "Search algorithm " << order.m_algorithm << dendl;

  IdentProcRelation relation;

  stk::search::coarse_search(range_vector, domain_vector, stk::search::BOOST_RTREE, comm, relation);

  if (0)
  {
    //for (unsigned i = 0; i < relation.size(); i++)
    for (unsigned i = 0; i < 10; i++)
    {
      std::cout << "relation[ " << i << "]= {" << relation[i].first << "} --> { " << relation[i].second << "}" << std::endl;
    }
  }
  //dw().m(LOG_SEARCH) << "relation " << relation << dendl;

  //dw().m(LOG_SEARCH) << stk::diag::pop;
}

#endif
}
}
