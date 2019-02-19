// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#define DO_TEST_TIME_SETS 0

#define NOMALLOC_ARRAY_CHECK_SIZES
#include <percept/NoMallocArray.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/fixtures/Fixture.hpp>

#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <boost/unordered_set.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>
#include <set>


namespace percept {
namespace unit_tests {

#define EXTRA_PRINT 0


//=============================================================================
//=============================================================================
//=============================================================================

//template<class K, class V, class MAP >
template<class SET >
static void setupSet(SET& set, unsigned N)
{
  for (unsigned i = 1; i <= N; i++)
  {
    set.insert( unsigned(i));
  }
}


template<class SET, class ITER >
static unsigned find1(SET& set,  ITER& i, unsigned key)
{
  i = set.find( key );
  return i != set.end() ? *i : 0 ;
}

template<class SET, class ITER >
struct FindSetItem2
{
  inline unsigned operator()(SET& set,  ITER& i, unsigned key)
  {
    i = set.find( key );
    return i != set.end() ? *i : 0 ;
    //return i->second;
  }
};



template<class SET, class ITER, class FUNC >
static double dot1(SET& set,  ITER& it, unsigned N, unsigned niter, FUNC& fm)
{
  //SET ::iterator ii = 0;
  double sum=0.0;
  unsigned ll=0u;
  for (unsigned iter = 0; iter < niter; iter++)
  {
    //std::cout << "iter= " << iter << std::endl;
    for (unsigned i = 0; i < N; i++)
    {
      //unsigned *pm = find1(set, it, i);
      unsigned pm = fm(set, it, i);
      ll+= pm;
      sum += (double)( pm)/(double(N));
      //set[i] = pm;
    }
  }
  return sum;
}


template<class SET, class ITER, class FUNC >
static void doTest(SET& set,  ITER& it, unsigned N, unsigned niter, FUNC& fm, std::string msg, double *times)
{
  EXCEPTWATCH;

  double t0 =  stk::cpu_time();
  setupSet(set, N);
  double t1 =  stk::cpu_time();
  std::cout << "settest:   setup time  = " << (t1-t0) << " [sec] for N = " << N << " for " << msg << std::endl;
  times[0] = t1-t0;

  double t2s =  stk::cpu_time();
  double dd= dot1(set, it, N, niter, fm);
  double t2e =  stk::cpu_time();
  times[1] = t2e-t2s;

  std::cout << "settest:  lookup time  = " << (t2e-t2s) << " [sec] for N = " << N << " for " << msg << " dd= " << dd << std::endl;

  double t3s =  stk::cpu_time();
  SET set3 = set;
  double t3e =  stk::cpu_time();
  times[2] = t3e-t3s;

  std::cout << "settest:  copy   time  = " << (t3e-t3s) << " [sec] for N = " << N << " for " << msg << std::endl;
}


TEST(time_sets, compare_different_sets)
{
#if DO_TEST_TIME_SETS

  EXCEPTWATCH;

  unsigned N = 10000000; // 10M
  //unsigned N = 100000; // 100K
  //unsigned N = 10; // 10M
  unsigned niter = 10;

  unsigned init_capacity = N;
  (void)init_capacity;

  typedef std::set<unsigned> std_set_type;
  std_set_type std_set1;
  std_set_type std_set2;
  std_set_type std_set3;
  //std::set<unsigned, unsigned*>::iterator std_set_it;
  std_set_type::iterator std_set_it1;
  std_set_type::iterator std_set_it2;

  typedef boost::unordered_set<unsigned> boost_set_type;
  boost_set_type boost_set1;
  boost_set_type boost_set2;
  boost_set_type boost_set3;
  boost_set_type::iterator boost_set_it1;
  boost_set_type::iterator boost_set_it2;

  //std::set<unsigned, unsigned *> std_set;

  std::ostringstream os;

  typedef double V6[6];
  os << "N : boost_setup : boost_lookup : boost_copy : std_setup : std_lookup : std_copy" << std::endl;
  for (unsigned iN = 1; iN <= N; iN *= 10)
    {
      V6 lt;
      {
        FindSetItem2< boost_set_type, boost_set_type::iterator > fm2;
        if (1) doTest(boost_set2, boost_set_it2, iN, niter, fm2, "boost_set, find(key)", &lt[0]);
      }

      {
        FindSetItem2< std_set_type, std_set_type::iterator > fm2;
        if (1) doTest(std_set2, std_set_it2, iN, niter, fm2, "std_set, find(key)", &lt[3]);
      }
      os << iN  ;
      for (int j=0; j < 6; j++)
        os << " : " << lt[j];
      os << std::endl;
    }
  std::cout << os.str() << std::endl;

  //doTest(boost_set, N, niter);

#endif

}

}
}
