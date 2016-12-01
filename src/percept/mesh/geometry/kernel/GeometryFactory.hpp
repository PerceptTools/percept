// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GEOMETRYFACTORY_HPP
#define GEOMETRYFACTORY_HPP

#include <string>
#include "GeometryKernel.hpp"
#include "MeshGeometry.hpp"
#include <percept/PerceptMesh.hpp>

namespace percept {

class GeometryFactory
{
public:
    GeometryFactory(GeometryKernel* kernel, MeshGeometry* geometry);
    virtual ~GeometryFactory();

    bool read_file(const std::string& filename, stk::mesh::MetaData* meta_data);
    bool read_file(const std::string& filename, PerceptMesh *mesh)
    {
      return read_file(filename, mesh->get_fem_meta_data());
    }

protected:
    GeometryKernel* geomKernel;
    MeshGeometry* geomDatabase;
};
}
#endif // GEOMETRYFACTORY_HPP
