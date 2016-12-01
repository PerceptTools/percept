// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#include <iostream>
#include <adapt/RefinementInfoByType.hpp>

namespace percept {

  typedef size_t MemorySizeType;

  struct MemoryMultipliers
  {
    MemorySizeType num_hex8;
    MemorySizeType num_tet4;
    MemorySizeType num_nodes;

    typedef MemorySizeType MemMultType;
    MemMultType mult_hex8;
    MemMultType mult_tet4;
    MemMultType mult_nodes;

    MemoryMultipliers(MemMultType mult_hex8=1490, MemMultType mult_tet4=702, MemMultType mult_nodes=0):
      //MemoryMultipliers(MemMultType mult_hex8=381, MemMultType mult_tet4=538, MemMultType mult_nodes=1017):
      num_hex8(0ul),
      num_tet4(0ul),
      num_nodes(0ul),
      mult_hex8(mult_hex8),
      mult_tet4(mult_tet4),
      mult_nodes(mult_nodes)
    {
    }

    void read_simple(std::string file_name)
    {
      std::ifstream file(file_name.c_str());
      if (file.good())
        file >> mult_hex8 >> mult_tet4 >> mult_nodes;
    }

    MemorySizeType estimate_memory()
    {
      return mult_nodes*num_nodes + mult_hex8*num_hex8 + mult_tet4*num_tet4;
    }

    MemorySizeType estimate_memory(std::vector<RefinementInfoByType>& refInfo, bool use_new=true)
    {
      num_hex8=0ul;
      num_tet4=0ul;
      num_nodes=0ul;

      for (unsigned i = 0; i < refInfo.size(); i++)
        {
          num_nodes= refInfo[0].m_numNewNodes;
          //std::cout << "irank, rank, m_numNewNodes, m_numNewElems= " << i << " " << refInfo[i].m_rank << " " << refInfo[i].m_numNewNodes
          //<< " " << refInfo[i].m_numNewElemsLast
          //<< std::endl;

          //             if (refInfo[i].m_rank == 0)
          //               {
          //                 num_nodes += refInfo[i].m_numNewNodes;
          //               }
          //             else
          {
            switch(refInfo[i].m_topology.getKey())
              {
              case shards::Hexahedron<8>::key:
                if (use_new)
                  num_hex8 += refInfo[i].m_numNewElemsLast;
                else
                  num_hex8 += refInfo[i].m_numOrigElemsLast;

                break;
              case shards::Tetrahedron<4>::key:
                if (use_new)
                  num_tet4 += refInfo[i].m_numNewElemsLast;
                else
                  num_tet4 += refInfo[i].m_numOrigElemsLast;
                break;
              default:
                break;
              }
          }
        }

      return estimate_memory();
    }

    static double MegaByte(MemorySizeType x) { return  ((double)x/1024.0/1024.0); }

    static void process_estimate(MemorySizeType tot_mem, PerceptMesh& eMesh, std::vector<RefinementInfoByType>& refInfo, std::string memory_multipliers_file, std::string input_file, bool use_new=true)
    {
      //const stk::ParallelMachine& comm = eMesh.get_bulk_data()->parallel();

      // this is a data gather pass
      if (tot_mem)
        {
          /*
            stk::mesh::Selector sel_locally_owned(eMesh.get_fem_meta_data()->locally_owned_part());
            stk::mesh::Selector sel_globally_shared(eMesh.get_fem_meta_data()->globally_shared_part());
            stk::mesh::Selector sel_universal(eMesh.get_fem_meta_data()->universal_part());

            std::vector<unsigned> count ;
            stk::mesh::count_entities( sel_universal, *eMesh.get_bulk_data(), count );

            unsigned nnodes = count[0];

            stk::ParallelMachine pm = eMesh.get_bulk_data()->parallel();
            stk::all_reduce( pm, stk::ReduceSum<1>( &nnodes ) );
          */
          MemoryMultipliers memMults;
          // FIXME, here's where we would read in some values for memMults from memory_multipliers_file
          if (memory_multipliers_file.size())
            memMults.read_simple(memory_multipliers_file);
          RefinementInfoByType::countCurrentNodes(eMesh, refInfo);
          MemorySizeType estMem = memMults.estimate_memory(refInfo, use_new);
          //std::cout << "tmp srk tot_mem = " << MegaByte(tot_mem) << " estMem= " << MegaByte(estMem) << std::endl;
          if (eMesh.get_rank() == 0)
            {
              if (1)
                std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " actualMem[MB]= " << MegaByte(tot_mem)
                          << " estMem[MB]= " << MegaByte(estMem)
                          << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;

              std::cout << "(*MemEstMM: {nn,ntet,nhex,mem,estMem} " << input_file << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << ", " << memMults.num_hex8 << ", " << MegaByte(tot_mem)
                        << ", " << MegaByte(estMem) << "}" << std::endl;
            }

        }
      else
        {
          // this is an estimate multipliers pass (computes memory using current multipliers)
          MemoryMultipliers memMults;
          // FIXME, here's where we would read in some values for memMults from memory_multipliers_file
          if (memory_multipliers_file.size())
            memMults.read_simple(memory_multipliers_file);
          RefinementInfoByType::countCurrentNodes(eMesh, refInfo);
          MemorySizeType estMem = memMults.estimate_memory(refInfo);
          //std::cout << "tmp srk tot_mem = " << MegaByte(tot_mem) << " estMem= " << MegaByte(estMem) << std::endl;
          if (eMesh.get_rank() == 0)
            {
              std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8
                //<< " memory[MB]= " << MegaByte(tot_mem)
                        << " estimatedMem[MB]= " << MegaByte(estMem)
                        << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;

              if (0)
                std::cout << "(*MemEstMM: " << input_file << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << "," << memMults.num_hex8 << "," << MegaByte(tot_mem)
                          << ", " << MegaByte(estMem) << "}" << std::endl;
            }

        }
    }
  };


}
