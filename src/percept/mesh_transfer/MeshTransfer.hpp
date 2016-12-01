// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_MeshTransfer_hpp
#define percept_MeshTransfer_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

namespace stk
{
  namespace io
  {
    class StkMeshIoBroker;
  }
}

  namespace percept
  {

    class MeshTransfer
    {
    public:
      MeshTransfer() : 
	src_mesh(),
	dst_mesh(),
	target_mesh(),
	dst_entity(),
	field_name(),
	thvec_name(),
	rzvec_name(),
	dst_field_name(),
	xrot(0),
	yrot(0),
	zrot(0),
	xtrans(0),
	ytrans(0),
	ztrans(0),
	clp(false)
	{}

      void run(int argc,  char** argv);

    private:
      void process_options();

      std::string src_mesh, dst_mesh, target_mesh;
      std::string dst_entity;
      std::string field_name, thvec_name, rzvec_name, dst_field_name;

      double xrot, yrot, zrot;
      double xtrans, ytrans, ztrans;

      Teuchos::CommandLineProcessor clp;
    };
      
  }//namespace percept

#endif
