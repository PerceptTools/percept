// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/function/internal/GenericFunction.hpp>

  namespace percept
  {
    std::ostream &operator<<(std::ostream& out,  GenericFunction& func)
    {
      out <<  "GenericFunction:: domain dims: " << func.getDomainDimensions() << " codomain dims: " << func.getCodomainDimensions();
      return out;
    }
  }
