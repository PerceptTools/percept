#ifndef percept_StructuredGridRefinerDef_hpp
#define percept_StructuredGridRefinerDef_hpp

namespace percept {

  /**
   * Refine a structured grid block
   *   Input the various sizing info @param input_sizes, @param output_sizes contains the new sizes.
   *   Input the coordinates @param input_xyz, output @param output_xyz
   *   Note: we assume the input/output_xyz are indexed as
   *     xyz(i0, i1, i2, {0,1,2})
   *   where loops are constructed as:
   *     // loop ordering
   *     const int L0 = loop_ordering[0], L1 = loop_ordering[1], L2 = loop_ordering[2];
   *     // access ordering
   *     const int A0 = access_ordering[0], A1 = access_ordering[1], A2 = access_ordering[2];
   *     int indx[3];
   *     for (indx[L2] = node_min[L2]; indx[L2] <= node_max[L2]; ++indx[L2])
   *       for (indx[L1] = node_min[L1]; indx[L1] <= node_max[L1]; ++indx[L1])
   *         for (indx[L0] = node_min[L0]; indx[L0] <= node_max[L0]; ++indx[L0])
   *           { do something with  xyz(indx[A0], indx[A1], indx[A2], {0,1,2}) }
   *
   *   e.g., loop_ordering = {0, 1, 2}, access_ordering = {0, 1, 2}:
   *      for (k = 0; k < nk; ++k)
   *        for (j = 0; j < nj; ++j)
   *          for (i = 0; i < ni; ++i)
   *            for (c = 0; c < 3; ++c)
   *               sum += xyz(i, j, k, c);
   *
   *   This way, the caller can decide the ordering of loops.
   *
   *   Note: node_min[], cell_max[] etc. are always interpreted as
   *      node_min[] = {node_imin, node_jmin, node_kmin}; etc.
   *
   *  If output_xyz == 0, only do the size computation step so the caller can
   *    construct output_xyz and call again.
   *
   *  TODO: boundary conditions interface
   */

  template<typename UInt, typename UInt64, typename Array4D>
  struct StructuredGridRefinerImpl {

    const std::string& block_name;
    const SGridSizes<UInt, UInt64>& input_sizes;
    SGridSizes<UInt, UInt64>& output_sizes;
    const std::array<UInt,3>& loop_ordering; // loop ordering
    const std::array<UInt,3>& access_ordering; // access ordering
    const Array4D& input_xyz;
    Array4D * output_xyz_p;
    int debug;

    // allowing for future 1-based indexing
    const unsigned m_index_base = 0;


    StructuredGridRefinerImpl(const std::string& block_name_in,
                              const SGridSizes<UInt, UInt64>& input_sizes_in,
                              SGridSizes<UInt, UInt64>& output_sizes_in,
                              const std::array<UInt,3>& loop_ordering_in, // loop ordering
                              const std::array<UInt,3>& access_ordering_in, // access ordering
                              const Array4D& input_xyz_in,
                              Array4D * output_xyz_p_in,
                              int debug_in)
      :
      block_name(block_name_in),
      input_sizes(input_sizes_in),
      output_sizes(output_sizes_in),
      loop_ordering(std::move(loop_ordering_in)),
      access_ordering(std::move(access_ordering_in)),
      input_xyz(input_xyz_in),
      output_xyz_p(output_xyz_p_in),
      debug(debug_in)
    {
    }

    void get_new_sizes() {
      // new sizes
      UInt sizes[3], new_sizes[3];
      for (UInt ii=0; ii < 3; ++ii)
        {
          sizes[ii] = input_sizes.node_max[ii] - input_sizes.node_min[ii] + 1;
          if (debug) std::cout << "StructuredGridRefiner::refine: sizes[" << ii << "]= " << sizes[ii] << std::endl;
          new_sizes[ii] = 2 * sizes[ii] - 1;
          output_sizes.node_min[ii] = 2*(input_sizes.node_min[ii] - m_index_base) + m_index_base;
          output_sizes.node_max[ii] = output_sizes.node_min[ii] + new_sizes[ii] - 1;
          output_sizes.cell_min[ii] = 2*(input_sizes.cell_min[ii] - m_index_base) + m_index_base;
          output_sizes.cell_max[ii] = output_sizes.cell_min[ii] + new_sizes[ii] - 2;
          output_sizes.node_size[ii] = new_sizes[ii];
          output_sizes.cell_size[ii] = new_sizes[ii] - 1;
          output_sizes.node_size_global[ii] = 2 * input_sizes.node_size_global[ii] - 1;
        }
    }

    // non-kokkos version
    void refine()
    {
      UInt indx[3]{0,0,0};
      const int L0 = loop_ordering[0], L1 = loop_ordering[1], L2 = loop_ordering[2];

      UInt index = 0;
      for (indx[L2] = input_sizes.node_min[L2]; indx[L2] <= input_sizes.node_max[L2]; ++indx[L2])
        {
          for (indx[L1] = input_sizes.node_min[L1]; indx[L1] <= input_sizes.node_max[L1]; ++indx[L1])
            {
              for (indx[L0] = input_sizes.node_min[L0]; indx[L0] <= input_sizes.node_max[L0]; ++indx[L0])
                {
                  (*this)(index);
                  ++index;
                }
            }
        }
    }

    // FIXME - need to find out if Kokkos requires including the last index of Array4D
    //    (i.e. the index over the xyz components)
    //    Also, need to generalize for general fields if we need to prolong fields as
    //    well as coordinates.
    inline
    void multi_dim_indices_from_index(const UInt& index, std::array<UInt,3>& indx) const
    {
      const int L0 = loop_ordering[0], L1 = loop_ordering[1], L2 = loop_ordering[2];
      const UInt sizes[3] = {
        1+ input_sizes.node_max[L0] - input_sizes.node_min[L0],
        1+ input_sizes.node_max[L1] - input_sizes.node_min[L1],
        1+ input_sizes.node_max[L2] - input_sizes.node_min[L2]
      };
      indx[L2] = input_sizes.node_min[L2] + (index / (sizes[L0]*sizes[L1] ));
      indx[L1] = input_sizes.node_min[L1] + ((index / sizes[L0]) % sizes[L1] );
      indx[L0] = input_sizes.node_min[L0] + (index % sizes[L0]);

      if (1)
        {
          UInt ii = indx[L0]-input_sizes.node_min[L0] + sizes[L0]*(indx[L1] - input_sizes.node_min[L1]) + sizes[L0]*sizes[L1]*(indx[L2] - input_sizes.node_min[L2]);
          VERIFY_OP_ON(ii, ==, index, "bad index");
        }
    }

    // KOKKOS_INLINE_FUNCTION
    inline
    void operator()(const UInt& index) const
    {
      bool ldebug = debug > 2;

      Array4D& output_xyz = *output_xyz_p;

      std::array<UInt,3> indx{{0,0,0}};
      //const int L0 = loop_ordering[0], L1 = loop_ordering[1], L2 = loop_ordering[2];
      const int A0 = access_ordering[0], A1 = access_ordering[1], A2 = access_ordering[2];

      multi_dim_indices_from_index(index, indx);

      UInt oindx[3]{2*(indx[0]-m_index_base)+m_index_base, 2*(indx[1]-m_index_base)+m_index_base, 2*(indx[2]-m_index_base)+m_index_base};

      // vertices
      for (unsigned ic = 0; ic < 3; ++ic)
        {
          output_xyz(oindx[A0], oindx[A1], oindx[A2], ic) = input_xyz(indx[A0], indx[A1], indx[A2], ic);
          if (ldebug) std::cout << "vertex(" << oindx[A0] << " " << oindx[A1] << " " << oindx[A2] << " " << ic << ") = " << output_xyz(oindx[A0], oindx[A1], oindx[A2], ic) << std::endl;
        }

      // edges
      const UInt delta[3][3]{{1,0,0},{0,1,0},{0,0,1}};
      for (unsigned ll = 0; ll < 3; ++ll)
        {
          if (indx[ll] < input_sizes.node_max[ll])
            {
              UInt oindxS[3]{oindx[0] + delta[ll][0], oindx[1] + delta[ll][1], oindx[2] + delta[ll][2]};
              UInt indxS[3]{indx[0], indx[1], indx[2]};
              UInt indxE[3]{indx[0] + delta[ll][0], indx[1] + delta[ll][1], indx[2] + delta[ll][2]};
              for (unsigned ic = 0; ic < 3; ++ic)
                {
                  output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) = 0.5*(input_xyz(indxS[A0], indxS[A1], indxS[A2], ic) + input_xyz(indxE[A0], indxE[A1], indxE[A2], ic));
                  if (ldebug) std::cout << "edge(" << oindxS[A0] << " " << oindxS[A1] << " " << oindxS[A2] << " " << ic << ") = " << output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) << std::endl;
                }
            }
        }

      // faces
      for (unsigned ll = 0; ll < 3; ++ll)
        {
          const unsigned addo[3][3] = { {0, 1, 1}, {1, 0, 1}, {1, 1, 0} };
          const unsigned lc[3] = {ll, (ll+1) % 3, (ll+2) % 3 };
          if (indx[lc[1]] < input_sizes.node_max[lc[1]] &&
              indx[lc[2]] < input_sizes.node_max[lc[2]])
            {

              UInt oindxS[3]{oindx[0]+addo[ll][0], oindx[1]+addo[ll][1], oindx[2]+addo[ll][2]};

              for (unsigned ic = 0; ic < 3; ++ic)
                {
                  output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) = 0;

                  for (unsigned ii = 0; ii < 2; ++ii)
                    {
                      for (unsigned jj = 0; jj < 2; ++jj)
                        {
                          const unsigned add[3][3] = { {0, ii, jj}, {ii, 0, jj}, {ii, jj, 0} };

                          UInt indxS[3]{indx[0] + add[ll][0], indx[1] + add[ll][1], indx[2] + add[ll][2]};
                          output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) += 0.25 * input_xyz(indxS[A0], indxS[A1], indxS[A2], ic);
                        }
                    }
                  if (ldebug) std::cout << "face(" << oindxS[A0] << " " << oindxS[A1] << " " << oindxS[A2] << " " << ic << ") = " << output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) << std::endl;
                }
            }
        }

      // centroid
      if (indx[0] < input_sizes.node_max[0]
          && indx[1] < input_sizes.node_max[1]
          && indx[2] < input_sizes.node_max[2])
        {
          UInt oindxS[3]{oindx[0] + 1, oindx[1] + 1, oindx[2] + 1};

          for (unsigned ic = 0; ic < 3; ++ic)
            {
              output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) = 0;
              for (unsigned i=0; i < 2; ++i)
                {
                  for (unsigned j = 0; j < 2; ++j)
                    {
                      for (unsigned k = 0; k < 2; ++k)
                        {
                          output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) += 0.125*input_xyz(indx[A0]+i, indx[A1]+j, indx[A2]+k, ic);
                        }
                    }
                }
              if (ldebug) std::cout << "centroid(" << oindxS[A0] << " " << oindxS[A1] << " " << oindxS[A2] << " " << ic << ") = " << output_xyz(oindxS[A0], oindxS[A1], oindxS[A2], ic) << std::endl;
            }
        }
    }

  };

}
#endif
