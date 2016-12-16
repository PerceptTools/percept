// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_StructuredGridRefiner_hpp
#define percept_StructuredGridRefiner_hpp

#include <percept/PerceptMesh.hpp>
#include <percept/BlockStructuredGrid.hpp>
#if !STK_PERCEPT_LITE
#  if defined(STK_BUILT_IN_SIERRA)
#    include <cgns/Iocgns_DatabaseIO.h>
#  else
#    include <Iocgns_DatabaseIO.h>
#  endif
#endif

namespace percept {

  class StructuredGridRefiner {

    PerceptMesh& m_eMesh;

    enum Algorithm {
      STRUCTURED_ONLY,
      UNSTRUCTURED_CONVERSION
    };
    Algorithm m_alg;

    // allowing for future 1-based indexing
    const unsigned m_index_base = 0;

    void do_refine_structured();

  public:
    std::shared_ptr<BlockStructuredGrid> m_input, m_output;
    int m_debug;

    StructuredGridRefiner(PerceptMesh& eMesh, int debug=0) :
      m_eMesh(eMesh), m_alg(STRUCTURED_ONLY), m_debug(debug)
    {
      m_input = eMesh.get_block_structured_grid();
      m_output.reset(new BlockStructuredGrid(&eMesh, 0));
    }

    // void read_cgns()
    // {
    //   m_input->read_cgns();
    // }

    void do_refine();

    void post_proc();

    // template<typename UInt, typename UInt64, typename Array4D>
    // void refine(const std::string& block_name, const UInt block_id_info[3],
    //             const SGridSizes<UInt, UInt64>& input_sizes,
    //             SGridSizes<UInt, UInt64>& output_sizes,
    //             const UInt loop_ordering[3], const UInt access_ordering[3],
    //             const Array4D& input_xyz,
    //             Array4D * output_xyz= 0);

    void print(std::ostream& out = std::cout, int level=0);

    // tmp - shows next steps for pulling in sideset info, etc.
    void double_structured_block(Ioss::Region *new_region,
                                 Ioss::StructuredBlock *oblock, size_t& num_node, size_t& num_cell);


  };


}

#include <percept/StructuredGridRefinerDef.hpp>


#endif
