// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/PerceptMesh.hpp>
#include <percept/BlockStructuredGrid.hpp>
#include <percept/MeshType.hpp>

#include <stk_util/parallel/CommSparse.hpp>

namespace percept {

  using Array4D = typename MTFieldImpl::Array4D;

  BlockStructuredGrid::BlockStructuredGrid(PerceptMesh* eMesh, Ioss::Region *region) : m_eMesh(eMesh), m_region(region)
  {
  }

  void BlockStructuredGrid::print(std::ostream& out, int level)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        m_sblocks[iblock]->print(out, level);
      }
  }

  void BlockStructuredGrid::read_cgns()
  {
    auto& blocks = m_region->get_structured_blocks();
    unsigned iblock=0;
    m_sblocks.resize(0);
    for (auto &block : blocks) {
      std::shared_ptr<StructuredBlock> ptr ( new StructuredBlock(m_eMesh, iblock, block, this) );
      ptr->read_cgns();
      m_sblocks.push_back(ptr);
      ++iblock;
    }
    std::cout << m_eMesh->rank() << " nblocks= " << m_sblocks.size() << std::endl;
  }

  void BlockStructuredGrid::register_field(const std::string& field, unsigned fourth_dim)
  {
    std::shared_ptr<MTFieldImpl> nf (new MTFieldImpl(field));
    m_fields[field] = nf;

    nf->m_block_fields.resize(m_sblocks.size());
    for (size_t ib=0; ib < m_sblocks.size(); ++ib)
      {
        nf->m_block_fields[ib] = m_sblocks[ib]->register_field(field, fourth_dim);
      }
  }

  void BlockStructuredGrid::create_pvd(const std::string& file_prefix, bool paraviewBroken)
  {
    if (m_eMesh->get_rank())
      return;

    std::ofstream out(file_prefix+".pvd");
    out << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    out << "  <Collection>\n";

    if (paraviewBroken)
      {
        int p_size = m_eMesh->get_parallel_size();
        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            for (int p_rank=0; p_rank < p_size; ++p_rank)
              {
                out << "      <DataSet part=\"" << iblock << "\" file=\"" << file_prefix << "." << iblock;
                if (p_size > 1 )
                  out << "." << p_size << "." << p_rank;
                out << ".vts\"/>" << std::endl;
              }
          }
      }
    else
      {
        std::string par_prefix = "";
        if (m_eMesh->get_parallel_size() > 1)
          par_prefix = "p";
        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            out << "      <DataSet part=\"" << iblock << "\" file=\"" << file_prefix << "." << iblock << "." + par_prefix + "vts\"/>" << std::endl;
            m_sblocks[iblock]->dump_vtk(file_prefix);
          }
      }

    out << "  </Collection>\n";
    out << "</VTKFile>\n";
  }

  /// creates a file, one per block, with pointers to all the parallel parts 
  ///   of the block in individual prefix_block.proc.vts files

  void BlockStructuredGrid::create_pvts(const std::string& file_prefix)
  {
    if (m_eMesh->get_parallel_size() == 1)
      return;

    stk::CommSparse comm_sparse(m_eMesh->parallel());
    //unsigned proc_rank = comm_sparse.parallel_rank();
    unsigned proc_size = comm_sparse.parallel_size();

    // send grid sizes to proc 0
    for (int stage=0; stage < 2; ++stage)
      {
        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            std::shared_ptr<StructuredBlock> sb = m_sblocks[iblock];

            for (int ii=0; ii < 3; ++ii)
              comm_sparse.send_buffer( 0 ).pack< unsigned > (sb->m_sizes.node_size[ii]);
          }

        if (stage == 0)
          {
            comm_sparse.allocate_buffers();
          }
        else
          {
            comm_sparse.communicate();
          }
      }

    if (m_eMesh->get_rank() == 0)
      {
        MDArrayUInt block_info(m_eMesh->get_parallel_size(), m_sblocks.size(), 3);

        for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
          {
            //if (from_proc == proc_rank)
            //  continue;
            stk::CommBuffer & recv_buffer = comm_sparse.recv_buffer( from_proc );

            for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
              {
                for (unsigned ii=0; ii < 3; ++ii)
                  recv_buffer.unpack< unsigned >( block_info(from_proc, iblock, ii) );
              }
            VERIFY_OP_ON(recv_buffer.remaining(), ==, false, "bad unpack");
          }

        for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
          {
            std::shared_ptr<StructuredBlock> sb = m_sblocks[iblock];

            std::ofstream out1(file_prefix+"."+toString(iblock)+".pvts");
            char buf0[1024];
            sprintf(buf0,
                    "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >\n"
                    "  <PStructuredGrid WholeExtent=\"0 %u 0 %u 0 %u\" GhostLevel=\"0\" >\n"
                    "      <PPoints>\n"
                    "         <PDataArray NumberOfComponents=\"3\" type=\"Float32\" />\n"
                    "      </PPoints>\n",
                    sb->m_sizes.node_size_global[0]-1, sb->m_sizes.node_size_global[1]-1, sb->m_sizes.node_size_global[2]-1);
            out1 << buf0;

            for(unsigned from_proc = 0; from_proc < proc_size; ++from_proc )
              {
                char buf1[1024];
                sprintf(buf1,
                        "      <PPiece Extent=\"0 %u 0 %u 0 %u\" Source=\"%s.%u.%u.vts\" />\n",
                        block_info(from_proc, iblock, 0) - 1,
                        block_info(from_proc, iblock, 1) - 1,
                        block_info(from_proc, iblock, 2) - 1,
                        file_prefix.c_str(), iblock, from_proc);
                out1 << buf1;
              }

            out1 <<  "  </PStructuredGrid>\n";
            out1 << "</VTKFile>​\n";
          }
      }
  }

  void BlockStructuredGrid::dump_vtk(const std::string& file_prefix)
  {
    double cpu = m_eMesh->start_cpu_timer();
    if (m_eMesh->get_rank() == 0)
      std::cout << "BlockStructuredGrid: saving VTK format file= " << file_prefix << " ... " << std::endl;

    bool paraviewBroken = true;
    create_pvd(file_prefix, paraviewBroken);
    // Paraview is broken, doesn't seem to like pvts, so, we just create more .vts pieces
    if (!paraviewBroken)
      create_pvts(file_prefix);

    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        m_sblocks[iblock]->dump_vtk(file_prefix);
      }

    cpu = m_eMesh->stop_cpu_timer(cpu);
    if (m_eMesh->get_rank() == 0)
      std::cout << "BlockStructuredGrid: ... cpu time= " << cpu << " sec." << std::endl;
  }

  void BlockStructuredGrid::get_nodes(std::vector<StructuredCellIndex>& nodes, unsigned offset)
  {
    nodes.resize(0);
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        const unsigned L0 = sgrid->m_loop_ordering[0], L1 = sgrid->m_loop_ordering[1], L2 = sgrid->m_loop_ordering[2];
        if (sgrid->is_empty())
          continue;

        unsigned sizes[3] = {sgrid->m_sizes.node_size[0], sgrid->m_sizes.node_size[1], sgrid->m_sizes.node_size[2]};

        unsigned indx[3]{0,0,0};

        for (indx[L2] = 0; indx[L2] < sizes[L2]-offset; ++indx[L2])
          {
            for (indx[L1] = 0; indx[L1] < sizes[L1]-offset; ++indx[L1])
              {
                for (indx[L0] = 0; indx[L0] < sizes[L0]-offset; ++indx[L0])
                  {
                    StructuredCellIndex node{{indx[0], indx[1], indx[2], iblock}};
                    nodes.push_back(node);
                  }
              }
          }
      }
  }

  void BlockStructuredGrid::get_elements(std::vector<StructuredCellIndex>& elements)
  {
    get_nodes(elements,  1);
  }

  struct SB_copy_field {
    Array4D& m_field_dest;
    const Array4D& m_field_src;

    SB_copy_field(Array4D& field_dest, const Array4D& field_src) : m_field_dest(field_dest), m_field_src(field_src) {}

    void run()
    {
      int64_t sz = m_field_dest.size();
      VERIFY_OP_ON(sz, ==, (int64_t)m_field_src.size(), "bad size");
      for (int64_t index = 0; index < sz; ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<SB_copy_field *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      m_field_dest[index] = m_field_src[index];
    }
  };

  void BlockStructuredGrid::copy_field(typename StructuredGrid::MTField* field_dest, typename StructuredGrid::MTField* field_src)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D& dest = (*field_dest->m_block_fields[iblock]);
        Array4D& src = (*field_src->m_block_fields[iblock]);
        SB_copy_field cf(dest, src);
        cf.run();
      }
  }


  /// axpbypgz calculates: z = alpha*x + beta*y + gamma*z
  struct SB_nodal_field_axpbypgz {
    const double m_alpha, m_beta, m_gamma;
    const Array4D& m_field_x;
    const Array4D& m_field_y;
    Array4D& m_field_z;

    SB_nodal_field_axpbypgz(double alp, const Array4D& field_x,
                            double bet, const Array4D& field_y,
                            double gam,  Array4D& field_z) :
      m_alpha(alp), m_beta(bet), m_gamma(gam),
      m_field_x(field_x), m_field_y(field_y), m_field_z(field_z) {}

    void run()
    {
      int64_t sz = m_field_x.size();
      VERIFY_OP_ON(sz, ==, (int64_t)m_field_y.size(), "bad size");
      VERIFY_OP_ON(sz, ==, (int64_t)m_field_z.size(), "bad size");
      for (int64_t index = 0; index < sz; ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<SB_nodal_field_axpbypgz *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      m_field_z[index] = m_alpha*m_field_x[index] + m_beta*m_field_y[index] + m_gamma*m_field_z[index];
    }
  };

  void BlockStructuredGrid::nodal_field_axpbypgz(double alpha, typename StructuredGrid::MTField* field_x,
                                                 double beta, typename StructuredGrid::MTField* field_y,
                                                 double gamma, typename StructuredGrid::MTField* field_z)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D& fx = (*field_x->m_block_fields[iblock]);
        Array4D& fy = (*field_y->m_block_fields[iblock]);
        Array4D& fz = (*field_z->m_block_fields[iblock]);
        SB_nodal_field_axpbypgz na(alpha, fx,
                                   beta,  fy,
                                   gamma, fz);
        na.run();
      }
  }

  /// axpby calculates: y = alpha*x + beta*y
  struct SB_nodal_field_axpby {
    const double m_alpha, m_beta;
    const Array4D& m_field_x;
    Array4D& m_field_y;

    SB_nodal_field_axpby(double alp, const Array4D& field_x,
                         double bet, Array4D& field_y) :
      m_alpha(alp), m_beta(bet),
      m_field_x(field_x), m_field_y(field_y) {}

    void run()
    {
      int64_t sz = m_field_x.size();
      VERIFY_OP_ON(sz, ==, (int64_t)m_field_y.size(), "bad size");
      for (int64_t index = 0; index < sz; ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<SB_nodal_field_axpby *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      m_field_y[index] = m_alpha*m_field_x[index] + m_beta*m_field_y[index];
    }
  };

  void BlockStructuredGrid::nodal_field_axpby(double alpha, typename StructuredGrid::MTField* field_x, double beta, typename StructuredGrid::MTField* field_y)
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D& fx = (*field_x->m_block_fields[iblock]);
        Array4D& fy = (*field_y->m_block_fields[iblock]);
        SB_nodal_field_axpby na(alpha, fx, beta, fy);
        na.run();
      }
  }


  struct SB_nodal_field_dot {
    long double m_sum;
    const Array4D& m_field_x;
    const Array4D& m_field_y;

    SB_nodal_field_dot(const Array4D& field_x,
                       const Array4D& field_y) :
      m_sum(0.0),
      m_field_x(field_x), m_field_y(field_y) {}

    void run()
    {
      int64_t sz = m_field_x.size();
      VERIFY_OP_ON(sz, ==, (int64_t)m_field_y.size(), "bad size");
      for (int64_t index = 0; index < sz; ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<SB_nodal_field_dot *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      m_sum += m_field_x[index]*m_field_y[index];
    }
  };

  long double BlockStructuredGrid::nodal_field_dot(typename StructuredGrid::MTField* field_x, typename StructuredGrid::MTField* field_y)
  {
    long double sum=0.0;
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D& fx = (*field_x->m_block_fields[iblock]);
        Array4D& fy = (*field_y->m_block_fields[iblock]);
        SB_nodal_field_dot na(fx, fy);
        na.run();
        sum += na.m_sum;
      }
    return sum;
  }

  /// set field to constant value
  struct SB_nodal_field_set_value {
    long double m_value;
    Array4D& m_field_x;

    SB_nodal_field_set_value(double val, Array4D& field_x) :
      m_value(val), m_field_x(field_x) {}

    void run()
    {
      int64_t sz = m_field_x.size();
      for (int64_t index = 0; index < sz; ++index)
        {
          (*this)(index);
        }
    }

    void operator()(int64_t& index) const
    {
      const_cast<SB_nodal_field_set_value *>(this)->operator()(index);
    }

    void operator()(int64_t& index)
    {
      m_field_x[index] = m_value;
    }
  };

  void BlockStructuredGrid::nodal_field_set_value(typename StructuredGrid::MTField* field_x, double value )
  {
    for (unsigned iblock=0; iblock < m_sblocks.size(); ++iblock)
      {
        std::shared_ptr<StructuredBlock> sgrid = m_sblocks[iblock];
        if (sgrid->is_empty())
          continue;
        Array4D& fx = (*field_x->m_block_fields[iblock]);
        SB_nodal_field_set_value na(value, fx);
        na.run();
      }
  }


  void BlockStructuredGrid::comm_fields(std::vector<const typename StructuredGrid::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    std::cout << "comm_fields not yet impl" << std::endl;
  }


  void BlockStructuredGrid::sum_fields(std::vector<typename StructuredGrid::MTField*>& fields, PerceptMesh *m_eMesh)
  {
    std::cout << "sum_fields not yet impl" << std::endl;
  }

  std::shared_ptr<BlockStructuredGrid>   BlockStructuredGrid::
  fixture_1(PerceptMesh* eMesh, std::array<unsigned,3> sizes)
  {
    std::shared_ptr<BlockStructuredGrid> bsg(new BlockStructuredGrid(eMesh,0));
    std::shared_ptr<StructuredBlock> sbi = StructuredBlock::fixture_1(eMesh, sizes, 0, 0, 0, bsg.get() );
    eMesh->set_block_structured_grid(bsg);
    bsg->m_sblocks.push_back(sbi);
    // FIXME - make this automatic, always have coordinates registered
    bsg->register_field("coordinates",3);
    return bsg;
  }

}

