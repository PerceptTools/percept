// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <iostream>
namespace testcpp11 {


  class U {
  public:
    static std::string s_opt;
  };
  std::string U::s_opt = {"test"};

  class A {
  public:
    std::string input{""};
    std::string output = {"out.e"};
    int option = -1;
    double hmesh_min_max_ave_factor[3] = {0,0,0};
    std::string abc {U::s_opt};
    std::string def {abc+"xx"};
  };
  void testcpp11()
  {
    A a;
    std::cout << "a= " << a.option << " " << a.input << std::endl;
  }

}
