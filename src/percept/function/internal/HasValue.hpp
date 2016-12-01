// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_function_HasValue_hpp
#define percept_function_HasValue_hpp

  namespace percept
  {
    template<typename ValueType>
    class HasValue
    {
    public:
      virtual ValueType& getValue() = 0;
      virtual void setValue(ValueType& ) = 0;
      virtual ~HasValue() {}
    };

    template<typename ValueType>
    class HasConstValue
    {
    public:
      virtual ValueType& getValue() = 0;
      virtual ~HasConstValue() {}
    };


  }

#endif
