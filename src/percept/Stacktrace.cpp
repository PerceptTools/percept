// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Stacktrace.hpp>

namespace percept {

#if DO_STACKTRACE_POS
#if DO_STACKTRACE_POS_USE_INT
  int Stacktrace::s_position = 0;
#else
  std::string Stacktrace::s_position = "";
#endif
#endif
}
