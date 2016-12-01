// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GEOMETRYKERNEL_OPENNURBS_HPP
#define GEOMETRYKERNEL_OPENNURBS_HPP

#include <opennurbs.h>
#include "GeometryKernel.hpp"

namespace percept {

class GeometryKernelOpenNURBS : public GeometryKernel
{
public:
    GeometryKernelOpenNURBS();
    virtual ~GeometryKernelOpenNURBS();

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities);
    virtual bool debug_dump_file(const std::string& file_name);

    virtual std::string get_attribute(GeometryHandle geom) const;

    virtual void snap_to(KernelPoint& point, GeometryHandle geom,
                         double *converged_tolerance = NULL,
                         double *uvw_computed = NULL,
                         double *uvw_hint = NULL, void *extra_hint = NULL);

    virtual void normal_at(KernelPoint& point, GeometryHandle geom, std::vector<double>& normal, void *extra_hint = NULL);

    virtual bool is_curve(GeometryHandle geom) const;

    virtual bool is_surface(GeometryHandle geom) const;

private:
    ONX_Model onModel;
};
}
#endif // GEOMETRYKERNEL_OPENNURBS_HPP
