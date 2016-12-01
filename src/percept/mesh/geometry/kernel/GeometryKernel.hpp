// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GEOMETRYKERNEL_HPP
#define GEOMETRYKERNEL_HPP

#include <vector>
#include <string>
#include <percept/PerceptMesh.hpp>

namespace percept {

typedef double* KernelPoint;
typedef int GeometryHandle;

class GeometryKernel
{
public:
    GeometryKernel() : m_spatialDim(3), m_debug(false), m_useFoundFaceMap(false) {}
    virtual ~GeometryKernel() {};

    virtual bool read_file(const std::string& file_name,
                           std::vector<GeometryHandle>& geometry_entities ) = 0;

    virtual void pre_process(percept::PerceptMesh *eMeshWithNode, const stk::mesh::PartVector& parts)
    {
    }

    virtual bool debug_dump_file(const std::string& file_name) { return true; }
    virtual std::string get_attribute(GeometryHandle geom) const = 0;

    virtual void snap_to(KernelPoint& point, GeometryHandle geom,
                    double *converged_tolerance = NULL,
                    double *uvw_computed = NULL,
                    double *uvw_hint = NULL, void *extra_hint = NULL) = 0;

    virtual void normal_at(KernelPoint& point, GeometryHandle geom, std::vector<double>& normal, void *extra_hint = NULL) = 0;

    virtual bool is_curve(GeometryHandle geom) const = 0;

    virtual bool is_surface(GeometryHandle geom) const = 0;
    void set_spatial_dim(int dim) { m_spatialDim = dim; }
    int get_spatial_dim() { return m_spatialDim; }
    void set_debug(bool debug) { m_debug = debug; }
    bool get_debug() { return m_debug; }
    void set_property(const std::string& key, const std::string& val) { m_propertyMap[key] = val; }
    std::string get_property(const std::string& key) { return m_propertyMap[key]; }
private:
  int m_spatialDim;
  bool m_debug;
public:
    bool m_useFoundFaceMap;
    typedef std::map<stk::mesh::Entity, stk::mesh::Entity> NodeToFoundFaceMap;
    NodeToFoundFaceMap m_nodeToFoundFaceMap;
    std::map<std::string, std::string> m_propertyMap;
};

}
#endif // GEOMETRYKERNEL_HPP
