// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef MeshType_hpp_
#define MeshType_hpp_

#include <array>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>
//#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
//#include <percept/PerceptMesh.hpp>

#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>
#include <percept/function/MDArray.hpp>
#include <percept/StructuredCellIndex.hpp>

namespace percept {

  //using StructuredCellIndex = std::array<unsigned,4>;  // i,j,k, block
  class PerceptMesh;

  struct SGridSelector {
    virtual bool operator()(StructuredCellIndex& indx)
    {
      VERIFY_MSG("not impl");
    }
  };

  struct SGridMeshGeometry {
    void normal_at(PerceptMesh *eMesh, StructuredCellIndex node, std::vector<double>& norm)
    {
      VERIFY_MSG("not impl");
    }

    int classify_node(StructuredCellIndex node_ptr, size_t& curveOrSurfaceEvaluator)
    {
      // FIXME
      return 2;
    }
  };

  struct MTFieldImpl {
    using Array4D = MDArray;
    std::vector<std::shared_ptr<MDArray> > m_block_fields;
    std::string m_name;
    MTFieldImpl(const std::string& nm) : m_name(nm) {}
    //MTFieldImpl(const std::string& nm, StructuredBlock *sgrid) : m_name(mm) {}
    std::string name() { return m_name; }
  };

  struct StructuredGrid {
    typedef SGridSelector MTSelector;
    typedef SGridMeshGeometry MTMeshGeometry;
    typedef StructuredCellIndex MTNode;
    typedef StructuredCellIndex MTElement;
    typedef StructuredCellIndex MTBucket;
#if STK_PERCEPT_LITE
    typedef double MTField;
#else
    typedef MTFieldImpl MTField;
#endif
    typedef int MTCellTopology;
    enum {NELEM_TYPES = 1 };
  };

  class MeshGeometry;

  struct STKMesh {
    typedef stk::mesh::Selector MTSelector;
    typedef MeshGeometry MTMeshGeometry;
    typedef stk::mesh::Entity MTNode;
    typedef stk::mesh::Entity MTElement;
    typedef stk::mesh::Bucket MTBucket;
    typedef stk::mesh::FieldBase MTField;
    typedef CellTopologyData MTCellTopology;
    enum {NELEM_TYPES = 10 };
  };

  template<typename MeshType>
  unsigned get_num_nodes(PerceptMesh *eMesh, typename MeshType::MTElement element);

  template<typename MeshType>
  const typename MeshType::MTNode *get_nodes(PerceptMesh *eMesh, typename MeshType::MTElement element, std::vector<typename MeshType::MTNode> *nodes = 0);

  template<typename MeshType>
  bool MTisGhostNode(PerceptMesh *m_eMesh, typename MeshType::MTNode node);

  template<typename MeshType>
  bool MTnode_locally_owned(PerceptMesh *m_eMesh, typename MeshType::MTNode node);

  template<typename MeshType>
  void MTcommFields(std::vector<const typename MeshType::MTField*>& fields, PerceptMesh *m_eMesh);

  template<typename MeshType>
  void MTsum_fields(std::vector<typename MeshType::MTField*>& fields, PerceptMesh *m_eMesh);

  /// gets @param field data from @param node into @param fld
  template<typename MeshType>
  void get_field(double *fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  template<typename MeshType>
  void get_field_new(double *fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// gets @param field data from @param node into @param fld[@param index]
  template<typename MeshType>
  void get_field(double *fld, unsigned size, int index, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// sets @param field data from @param fld into @param node
  template<typename MeshType>
  void set_field(const double * fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// sets @param field data from @param fld[@param index] into @param node
  template<typename MeshType>
  void set_field(const double * fld, unsigned size, int index, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);


  template<typename T, size_t N>
  std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
  {
    for (unsigned i=0; i < N; ++i)
      os << arr[i] << (i == N-1 ? "" : " ");
    return os;
  }


}
#endif
