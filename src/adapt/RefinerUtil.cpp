// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <adapt/SerializeNodeRegistry.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/Refiner.hpp>
#include <adapt/NodeRegistry.hpp>
#include <percept/PEnums.hpp>

#include <boost/unordered_set.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

namespace percept {

  using namespace std;
  using namespace percept;

#define EXTRA_PRINT_UR_GETBLOCKS 0

  static std::string strip_mult(const std::string& inp)
  {
    size_t pc = inp.find(":");
    if (pc == std::string::npos)
      return inp;
    else
      return inp.substr(0,pc);
  }

  BlockNamesType RefinerUtil::getBlockNames(const std::string& block_name_0, unsigned proc_rank, percept::PerceptMesh& eMesh, const std::string& geomFile)
  {
    std::string block_name = block_name_0;

    BlockNamesType blocks(percept::EntityRankEnd+1u);
    if (block_name.length() == 0)
      return blocks;

    if (block_name.substr(0, 5) == "file:")
      {
        //if (1) throw std::runtime_error("file: option Not implemented");
        std::string fileName = block_name.substr(5, block_name.length()-5);
        for (int iter=0; iter < 1000; ++iter)
          {
            if (fileName[0] == ' ')
              fileName = fileName.substr(1, fileName.length()-1);
            else
              break;
          }
        if (eMesh.get_rank() == 0)
          {
            std::cout << "block_names processing, found 'file:', will read from file = " << fileName << std::endl;
          }
        std::ifstream file(fileName.c_str());
        int line=0;
        block_name = "";
        std::string str_line;
        while(std::getline(file, str_line))
          {
            std::string block = str_line;
            if (eMesh.get_rank() == 0)
              std::cout << "file: line= " << line << " block= " << block << std::endl;
            if (block[0] != '#' && block.length())
              {
                if (line) block_name += ",";
                block_name += block;
                ++line;
              }
          }
        if (eMesh.get_rank() == 0)
          {
            std::cout << "block_names processing, after file read string is: " << block_name << std::endl;
          }
      }

    if (1)
      {
        std::string names = block_name;

        // pre-process to look for ".." range indicator

        std::string new_names = names;
        new_names = "";
        while(1)
          {
            if (!names.length())
              break;
            size_t ipos = names.find(',');
            bool last_one =  (ipos == std::string::npos);

            {
              std::string n1 = (last_one ? names : names.substr(0, ipos) );
              bool inc = true;
              //bool exc = false;
              if ('-' == n1[0])
                {
                  //exc = true;
                  inc = false;
                }
              else if ('+' == n1[0])
                {
                }
              else
                {
                  n1 = "+" + n1;
                }
              std::string plus_or_minus = (inc?"+":"-");
              std::string n2 = n1.substr(1, n1.length()-1);
              std::string id_string_start = "";
              std::string id_string_end = "";
              // leave open the possibility for other identifiers for range
              std::string dotdot = "..";
              int dotdot_len = dotdot.length();
              size_t pos_dotdot = n1.find(dotdot);
              if (pos_dotdot != std::string::npos)
                {
                  if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp with .., n1= " << n1 << " n2= " << n2 << std::endl;

                  if (n1.length() > 6 && n1.substr(1,6) == "block_")
                    {
                      // +block_123..block_125
                      // 0123456789^1234567890
                      id_string_start = n1.substr(7, pos_dotdot-7);
                      id_string_end = n1.substr(pos_dotdot+dotdot_len+6, n1.length()-(pos_dotdot+dotdot_len+6));
                    }
                  else if (n1.length() > 8 && n1.substr(1,8) == "surface_")
                    {
                      // error
                    }
                  else
                    {
                      // +12..45
                      // 012^456
                      //std::cout << "tmp pos_dotdot= " << pos_dotdot << std::endl;

                      id_string_start = n1.substr(1, pos_dotdot-1);
                      id_string_end = n1.substr(pos_dotdot+dotdot_len+0, n1.length()-(pos_dotdot+dotdot_len+0));
                    }

                  std::string mult = "";
                  size_t NPOS = std::string::npos;
                  size_t cpos = id_string_end.find(":");
                  if (cpos != std::string::npos)
                    {
                      size_t xpos = id_string_end.find("x");
                      if (xpos == NPOS)
                        xpos = id_string_end.find("X");
                      VERIFY_OP_ON(xpos, !=, NPOS, "missing x or X in blocks specifier");
                      mult = id_string_end.substr(cpos);
                      id_string_end = id_string_end.substr(0, cpos);
                    }
                  if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp with .., id_string_start= " << id_string_start << " id_string_end= " << id_string_end << " mult= " << mult << std::endl;

                  int id_start = 0;
                  int id_end = 0;
                  try {
                    id_start = boost::lexical_cast<int>(id_string_start);
                    id_end = boost::lexical_cast<int>(id_string_end);
                  }
                  catch (std::exception& X)
                    {
                      std::cout << "RefinerUtil::getBlockNames: exception: " << X.what() << std::endl;
                      std::cout << "RefinerUtil::getBlockNames: invalid range syntax in block_name: with .., id_string_start= "
                                << id_string_start << " id_string_end= " << id_string_end << std::endl;
                      throw std::runtime_error("invalid input syntax");
                    }
                  catch ( const std::exception * X )
                    {
                      std::cout << "RefinerUtil::getBlockNames: exception: " << X->what() << std::endl;
                      std::cout << "RefinerUtil::getBlockNames: invalid range syntax in block_name: with .., id_string_start= "
                                << id_string_start << " id_string_end= " << id_string_end << std::endl;
                      throw std::runtime_error("invalid input syntax");
                    }
                  catch( ... )
                    {
                      throw std::runtime_error("invalid input syntax");
                    }
                  if (EXTRA_PRINT_UR_GETBLOCKS)
                    {
                      std::cout << "tmp with .., id_string_start= " << id_string_start << " id_string_end= " << id_string_end << std::endl;
                      std::cout << "tmp with .., id_start= " << id_start << " id_end= " << id_end << std::endl;
                    }

                  for (int id=id_start; id <= id_end; id++)
                    {
                      new_names += plus_or_minus + (boost::lexical_cast<std::string>(id)) + mult + (id == id_end ? "" : ",");
                    }
                  if (!last_one)
                    new_names += ",";
                  if (EXTRA_PRINT_UR_GETBLOCKS)
                    std::cout << "tmp new_names with .. = " << new_names << std::endl;
                }
              else
                {
                  new_names += n1 + (last_one? "":",");
                  if (EXTRA_PRINT_UR_GETBLOCKS)
                    std::cout << "tmp new_names without .. = " << new_names << std::endl;
                }
              if (last_one)
                {
                  break;
                }
              else
                {
                  names = names.substr(ipos+1, names.length()-(ipos+1));
                }
            }
          }
        if (EXTRA_PRINT_UR_GETBLOCKS)
          std::cout << "tmp new_names after .. (range) processing = " << new_names << std::endl;

        names = new_names;
        std::string names_save = names;

        // post process to remove +name if -name exists
        bool has_plus = names.find('+') != std::string::npos;
        new_names = "";
        while(1)
          {
            if (!names.length())
              break;
            size_t pos_comma = names.find(',');
            bool last_one =  (pos_comma == std::string::npos);

            std::string n1 = (last_one ? names : names.substr(0, pos_comma) );

            bool inc = true;
            //bool exc = false;
            if ('-' == n1[0])
              {
                //exc = true;
                inc = false;
              }
            else if ('+' == n1[0])
              {
              }
            else
              {
                //error
              }
            std::string n2 = n1.substr(1, n1.length()-1);
            if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "n1= " << n1 << " n2= " << n2 << " strip_mult(n2)= " << strip_mult(n2) << std::endl;
            n2 = strip_mult(n2);
            // if (n2 == "block_910")
            //   {
            //     std::cout << "n1= " << n1 << " n2= " << n2 << std::endl;
            //   }

            if (inc)
              {

                size_t jpos = names_save.find("-"+n2);
                bool add_it = true;
                if (jpos != std::string::npos)
                  {
                    // don't add it
                    add_it = false;
                    // check for full match
                    std::string ns;
                    for (size_t ipos=jpos; ipos < names_save.length(); ++ipos)
                      {
                        if (names_save[ipos] == ',')
                          break;
                        ns += names_save[ipos];
                      }
                    std::string ns2 = ns.substr(1, ns.length()-1);
                    if (ns2 != n2)
                      add_it = true;
                    if (EXTRA_PRINT_UR_GETBLOCKS)
                      {
                        std::cout << "ns= " << ns << " ns2= " << ns2 << " n2= " << n2 << std::endl;
                      }
                  }
                if (add_it)
                  {
                    new_names += n1 + (last_one? "":",");
                  }
              }
            else
              {
                if (!has_plus)
                  new_names += n1 + (last_one? "":",");
              }

            if (last_one)
              {
                break;
              }
            else
              {
                names = names.substr(pos_comma+1, names.length()-(pos_comma+1));
              }
          }


        if (EXTRA_PRINT_UR_GETBLOCKS)
          std::cout << "tmp new_names after post-proc to remove +name if -name exists= " << new_names << std::endl;
        if (new_names.length() && !proc_rank)
          {
            std::cout << "RefinerUtil:: --block_name option after processing for removing -name= " << new_names << std::endl;
          }

        // final step
        names = new_names;
        while(1)
          {
            if (!names.length())
              break;
            size_t ipos = names.find(',');
            bool last_one =  (ipos == std::string::npos);

            {
              std::string n1 = (last_one ? names : names.substr(0, ipos) );

              bool inc = true;
              //bool exc = false;
              if ('-' == n1[0])
                {
                  //exc = true;
                  inc = false;
                }
              else if ('+' == n1[0])
                {
                }
              else
                {
                  n1 = "+" + n1;
                }
              std::string n2 = n1.substr(1, n1.length()-1);

              //std::cout << "n1= " << n1 << " n2= " << n2 << std::endl;
              if (n1.length() > 6 && n1.substr(1,6) == "block_")
                blocks[stk::topology::ELEMENT_RANK].push_back(n1);
              else if (n1.length() > 8 && n1.substr(1,8) == "surface_")
                blocks[eMesh.face_rank()].push_back(n1);
              else
                {
                  std::string pm = (inc?"+":"-");
                  blocks[stk::topology::ELEMENT_RANK].push_back(pm+"block_"+n2);
                }
              if (last_one)
                {
                  break;
                }
              else
                {
                  names = names.substr(ipos+1, names.length()-(ipos+1));
                }
            }
          }
        if (EXTRA_PRINT_UR_GETBLOCKS)
          std::cout << "tmp RefinerUtil::getBlockNames: blocks = " << blocks << std::endl;
      }

    // change so all -block get removed and replaced with + only
    for (unsigned irank=eMesh.side_rank(); irank <= eMesh.element_rank(); irank++)
      {
        bool found_minus = false, found_plus = false;
        for (unsigned ib=0; ib < blocks[irank].size(); ++ib)
          {
            if (blocks[irank][ib][0] == '-')
              {
                found_minus = true;
              }
            if (blocks[irank][ib][0] == '+')
              {
                found_plus = true;
              }
          }
        if (found_minus && found_plus)
          {
            std::ostringstream errmsg;
            errmsg << "You've specified an inconsistent combination of block names, mixing -block_n... with +block_m...\n";
            errmsg << "You can only specify -block_n if there are no +block_m's or if the +block_m's contain the -block_n\n";
            errmsg << "   to be deleted, e.g. +1,2,3.7,-4\n";
            errmsg << "Curent processed blocks list = " << blocks;
            errmsg << "\nSpecified on input: " << block_name << std::endl;
            throw std::runtime_error(errmsg.str());
          }

        if (found_minus)
          {
            std::vector<std::string> new_names;
            stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
            for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
              {
                stk::mesh::Part *  part = *i_part ;
                if (eMesh.is_auto_or_geom_part(geomFile, part))
                  continue;

                if (part->primary_entity_rank() == irank)
                  {
                    std::string partNamePlus = "+" + part->name();
                    std::string partNameMinus = "-" + part->name();
                    bool remove_it = false;
                    for (unsigned ib=0; ib < blocks[irank].size(); ++ib)
                      {
                        if (partNameMinus == blocks[irank][ib])
                          {
                            remove_it = true;
                            break;
                          }
                      }
                    if (!remove_it)
                      {
                        new_names.push_back(partNamePlus);
                      }
                  }
              }
            if (EXTRA_PRINT_UR_GETBLOCKS && !eMesh.get_rank())
              std::cout << "RefinerUtil::getBlockNames: irank= " << irank << " old blocks = \n" << blocks[irank] << "\n new= \n" << new_names << std::endl;
            blocks[irank] = new_names;
          }
      }

    return blocks;
  }

  /**
   * This method looks for surfaces that share nodes with the blocks specified in @param blocks and if it finds
   * any surfaces (sidesets), they are added to the blocks so they get refined properly.
   * TODO: If a surface is shared by more than one block, an error is thrown.
   * @deprecated - @see checkBreakPatternValidityAndBuildRanks_1() below
   */
  BlockNamesType RefinerUtil::correctBlockNamesForPartPartConsistency(percept::PerceptMesh& eMesh, BlockNamesType& blocks)
  {
    if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "RefinerUtil::correctBlockNamesForPartPartConsistency..." << std::endl;

    if (blocks[stk::topology::ELEMENT_RANK].size() == 0)
      return blocks;

    stk::mesh::EntityRank subDimRank = eMesh.side_rank();

    stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
    for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
      {
        stk::mesh::Part *  part = *i_part ;

        for (stk::mesh::PartVector::iterator i_surfacePart = all_parts.begin(); i_surfacePart != all_parts.end(); ++i_surfacePart)
          {
            stk::mesh::Part *  surfacePart = *i_surfacePart ;
            if ( stk::mesh::is_auto_declared_part(*surfacePart) )
              continue;

            const CellTopologyData * part_cell_topo_data = eMesh.get_cell_topology(*surfacePart);
            CellTopology surf_topo(part_cell_topo_data);
            //if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk surfacePart= " << surfacePart->name() << " topo= " << (part_cell_topo_data?surf_topo.getName() : "NULL") << std::endl;

            if (part_cell_topo_data && part->primary_entity_rank() == stk::topology::ELEMENT_RANK && surfacePart->primary_entity_rank() == subDimRank)
              {
                std::string partNamePlus = "+" + part->name();
                std::vector<std::string>::iterator partInBlocks = std::find(blocks[stk::topology::ELEMENT_RANK].begin(), blocks[stk::topology::ELEMENT_RANK].end(), partNamePlus);
                // if this part is not in the blocks list, skip it
                if (partInBlocks == blocks[stk::topology::ELEMENT_RANK].end())
                  {
                    //if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk skipping part= " << partNamePlus << std::endl;
                    continue;
                  }
                std::string surfacePartNamePlus = "+" + surfacePart->name();
                std::vector<std::string>::iterator surfacePartInBlocks = std::find(blocks[subDimRank].begin(), blocks[subDimRank].end(), surfacePartNamePlus);
                // if this surface is already in the list, skip it
                if (surfacePartInBlocks != blocks[subDimRank].end())
                  {
                    //if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk skipping surf= " << surfacePartNamePlus << std::endl;
                    continue;
                  }
                bool isBoundarySurface= eMesh.isBoundarySurface(*part, *surfacePart);

                if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk isBoundarySurface for part/surf= " << part->name() << " / " << surfacePart->name() << " = " << isBoundarySurface << std::endl;
                if (isBoundarySurface)
                  {
                    if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp part [" << part->name() << "] shares sideset [" << surfacePart->name() << "]" << std::endl;
                    blocks[subDimRank].push_back(std::string("+"+surfacePart->name()));
                  }
                else
                  {
                    //std::cout << "tmp part [" << part->name() << "] doesn't shares sideset [" << surfacePart->name() << "]" << std::endl;
                  }
              }
          }
      }
    if (0) std::cout << "tmp RefinerUtil::correctBlockNamesForPartPartConsistency: blocks = " << blocks << std::endl;
    return blocks;
  }

  static void make_parallel_consistent(PerceptMesh& eMesh, std::vector<std::string>& blocks)
  {
    char buf[1024];
    stk::CommAll comm_all(eMesh.parallel());

    unsigned proc_rank = comm_all.parallel_rank();
    unsigned proc_size = comm_all.parallel_size();

    // pack
    for (unsigned pr=0; pr < proc_size; pr++)
      {
        comm_all.send_buffer( pr ).pack< unsigned > (blocks.size());
        for (unsigned ib=0; ib < blocks.size(); ib++)
          {
            comm_all.send_buffer( pr ).pack< unsigned > ( blocks[ib].length());
            comm_all.send_buffer( pr ).pack< char > (blocks[ib].c_str(), blocks[ib].length());
          }
      }
    // allocateBuffers
    {
      bool local = true;
      unsigned num_msg_bounds = proc_size < 4 ? proc_size : proc_size/4 ;
      bool global = comm_all.allocate_buffers(num_msg_bounds , false, local );
      if ( not global )
        {
          std::cout << "P[" << proc_rank << "] : not global" << std::endl;
          throw std::runtime_error("not global");
        }
    }
    // pack
    for (unsigned pr=0; pr < proc_size; pr++)
      {
        comm_all.send_buffer( pr ).pack< unsigned > (blocks.size());
        for (unsigned ib=0; ib < blocks.size(); ib++)
          {
            comm_all.send_buffer( pr ).pack< unsigned > ( blocks[ib].length());
            comm_all.send_buffer( pr ).pack< char > (blocks[ib].c_str(), blocks[ib].length());
          }
      }
    // communicate
    comm_all.communicate();

    // unpack
    std::set<std::string> blocks_new(blocks.begin(), blocks.end());
    for (unsigned pr=0; pr < proc_size; pr++)
      {
        unsigned bsize=0;
        comm_all.recv_buffer( pr ).unpack< unsigned > (bsize);
        for (unsigned ib=0; ib < bsize; ib++)
          {
            unsigned len=0;
            comm_all.recv_buffer( pr ).unpack< unsigned > ( len );
            comm_all.recv_buffer( pr ).unpack< char > (buf, len);
            std::string str(buf, len);
            blocks_new.insert(str);
          }
      }
    blocks.resize(0);
    blocks.assign(blocks_new.begin(), blocks_new.end());
    std::sort(blocks.begin(), blocks.end());
    //std::cout << "make_parallel_consistent blocks= " << blocks << std::endl;
  }

  /**
   * This method looks for surfaces of blocks specified in @param blocks and if it finds
   * any surfaces (sidesets), they are added to the blocks so they get refined properly.
   * If a surface is shared by more than one block, and both blocks are not in the @param blocks list,
   *   an error is thrown.
   *
   * algorithm:
   *   for sides in side_part_map
   *     if all elem_part in side_part_map[side] are in blocks, ok, else throw
   *     if side not in blocks[side_rank], add it
   *
   */

#define DEBUG_CORRECT_BLOCKS_1 0

  BlockNamesType RefinerUtil::correctBlockNamesForPartPartConsistency_1(percept::PerceptMesh& eMesh, BlockNamesType& blocks, const std::string& geomFile)
  {
    if (blocks[stk::topology::ELEMENT_RANK].size() == 0)
      return blocks;

    SidePartMap side_part_map;
    Refiner::get_side_part_relations(eMesh, false, side_part_map);

    stk::mesh::EntityRank subDimRank = eMesh.side_rank();

    SidePartMap::iterator iter;
    for (iter = side_part_map.begin(); iter != side_part_map.end(); iter++)
      {
        stk::mesh::Part *side_part = iter->first;
        std::string side_part_name = side_part->name();
        const stk::mesh::PartVector *side_pv  = side_part->attribute<stk::mesh::PartVector>();
        stk::mesh::PartVector side_pv1;
        if (side_pv) side_pv1 = *side_pv;
        side_pv1.push_back(side_part);

        stk::mesh::PartVector elem_pv1;

        const stk::mesh::PartVector& epv = iter->second;
        for (unsigned iepv=0; iepv < epv.size(); iepv++)
          {
            stk::mesh::Part *elem_part = epv[iepv];
            std::string elem_part_name = elem_part->name();
            const stk::mesh::PartVector *elem_pv  = elem_part->attribute<stk::mesh::PartVector>();
            elem_pv1.push_back(elem_part);
            if (elem_pv) elem_pv1.insert(elem_pv1.end(), elem_pv->begin(), elem_pv->end());
          }

        for (unsigned iside_pv = 0; iside_pv < side_pv1.size(); iside_pv++)
          {
            stk::mesh::Part *surfacePart = side_pv1[iside_pv];
            const CellTopologyData * surfacePart_topo_data = eMesh.get_cell_topology(*surfacePart);
            CellTopology surf_topo(surfacePart_topo_data);
            if (!surfacePart_topo_data)
              continue;
            //if (DEBUG_CORRECT_BLOCKS_1) std::cout << "tmp srk surfacePart= " << surfacePart->name() << " topo= " << (surfacePart_topo_data?surf_topo.getName() : "NULL") << std::endl;

            //if (eMesh.is_auto_or_geom_part(geometry_file_name, *surfacePart))
            //  continue;

            bool at_least_one_elem_part_in_block_names = false;
            for (unsigned ielem_pv = 0; ielem_pv < elem_pv1.size(); ielem_pv++)
              {
                stk::mesh::Part *elementPart = elem_pv1[ielem_pv];
                //if (eMesh.is_auto_or_geom_part(geometry_file_name, *elementPart))
                //  continue;
                stk::mesh::Part *part = elementPart;
                std::string partNamePlus = "+" + part->name();

                std::vector<std::string>::iterator partInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), partNamePlus);
                if (partInBlocks != blocks[eMesh.element_rank()].end())
                  {
                    at_least_one_elem_part_in_block_names = true;
                    break;
                  }
              }
            if (DEBUG_CORRECT_BLOCKS_1) std::cout << "tmp srk surfacePart= " << surfacePart->name()
                                                  << " at_least_one_elem_part_in_block_names= " << at_least_one_elem_part_in_block_names
                                                  << " topo= " << (surfacePart_topo_data?surf_topo.getName() : "NULL")
                                                  << std::endl;

            if (!at_least_one_elem_part_in_block_names)
              continue;

            // if at_least_one_elem_part_in_block_names then all must be in there
            for (unsigned ielem_pv = 0; ielem_pv < elem_pv1.size(); ielem_pv++)
              {
                stk::mesh::Part *elementPart = elem_pv1[ielem_pv];
                //if (eMesh.is_auto_or_geom_part(geometry_file_name, *elementPart))
                //  continue;
                stk::mesh::Part *part = elementPart;
                std::string partNamePlus = "+" + part->name();

                std::vector<std::string>::iterator partInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), partNamePlus);
                bool found_it = partInBlocks != blocks[eMesh.element_rank()].end();

                if (DEBUG_CORRECT_BLOCKS_1) std::cout << "tmp srk surfacePart= " << surfacePart->name() << " elementPart= " << elementPart->name() << " found_it= " << found_it
                                                      << " part->primary_entity_rank()= " << part->primary_entity_rank()
                                                      << " surfacePart->primary_entity_rank()= " << surfacePart->primary_entity_rank()
                                                      << std::endl;

                if (!found_it)
                  {
                    std::ostringstream errmsg;
                    errmsg << "ERROR: correctBlockNamesForPartPartConsistency_1: found a surface (" + surfacePart->name() + ") that shares a block (" + part->name() + ") that is not in the specified\n"
                      " list of blocks to be refined.  Re-run by adding the missing block(s), which are = \n";
                    for (unsigned jelem_pv = 0; jelem_pv < elem_pv1.size(); jelem_pv++)
                      {
                        stk::mesh::Part *jelementPart = elem_pv1[jelem_pv];
                        //if (eMesh.is_auto_or_geom_part(geometry_file_name, *elementPart))
                        //  continue;
                        stk::mesh::Part *jpart = jelementPart;
                        std::string jpartNamePlus = "+" + jpart->name();
                        std::vector<std::string>::iterator jpartInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), jpartNamePlus);
                        if (jpartInBlocks == blocks[eMesh.element_rank()].end())
                          {
                            errmsg << " " << jpart->name();
                          }
                      }
                    errmsg << std::endl;
                    //std::cout << "\n\n" << errmsg.str() << "\n\n" << std::endl;
                    //throw std::runtime_error(errmsg.str());
                  }


                if (surfacePart_topo_data && part->primary_entity_rank() == eMesh.element_rank() && surfacePart->primary_entity_rank() == subDimRank)
                  {
                    std::string surfacePartNamePlus = "+" + surfacePart->name();
                    std::vector<std::string>::iterator surfacePartInBlocks = std::find(blocks[subDimRank].begin(), blocks[subDimRank].end(), surfacePartNamePlus);
                    // if this surface is already in the list, skip it
                    if (surfacePartInBlocks != blocks[subDimRank].end())
                      {
                        //if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk skipping surf= " << surfacePartNamePlus << std::endl;
                        continue;
                      }
                    bool isBoundarySurface= true; // by definition, the side part map is map of sides to shared element parts

                    if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp srk isBoundarySurface for part/surf= " << part->name() << " / " << surfacePart->name() << " = " << isBoundarySurface << std::endl;
                    if (isBoundarySurface)
                      {
                        if (EXTRA_PRINT_UR_GETBLOCKS) std::cout << "tmp part [" << part->name() << "] shares sideset [" << surfacePart->name() << "]" << std::endl;
                        blocks[subDimRank].push_back(std::string("+"+surfacePart->name()));
                      }
                    else
                      {
                        //std::cout << "tmp part [" << part->name() << "] doesn't shares sideset [" << surfacePart->name() << "]" << std::endl;
                      }
                  }
              }
          }
      }

    // add in geometry parts
    if (1)
      {
        stk::mesh::PartVector all_parts = eMesh.get_fem_meta_data()->get_parts();
        for (stk::mesh::PartVector::iterator i_part = all_parts.begin(); i_part != all_parts.end(); ++i_part)
          {
            stk::mesh::Part *  part = *i_part ;
            if (stk::mesh::is_auto_declared_part(*part))
              continue;

            if (part->primary_entity_rank() == eMesh.element_rank() )
              {
                if (part->name().find("tbs_") != std::string::npos || part->name().find("tbc_") != std::string::npos)
                  {
                    for (stk::mesh::PartVector::iterator j_part = all_parts.begin(); j_part != all_parts.end(); ++j_part)
                      {
                        stk::mesh::Part * block = *j_part ;
                        if (stk::mesh::is_auto_declared_part(*block))
                          continue;
                        if (block->name().find("tbs_") != std::string::npos || block->name().find("tbc_") != std::string::npos)
                          continue;

                        std::string jpartNamePlus = "+" + block->name();
                        std::vector<std::string>::iterator jpartInBlocks = std::find(blocks[eMesh.element_rank()].begin(), blocks[eMesh.element_rank()].end(), jpartNamePlus);

                        if (jpartInBlocks != blocks[eMesh.element_rank()].end())
                          {
                            bool allow_single_node_sharing = false;
                            bool isBoundarySurface= eMesh.isBoundarySurface(*block, *part, allow_single_node_sharing);

                            if (isBoundarySurface)
                              {
                                // FIXME check if already there?
                                blocks[eMesh.element_rank()].push_back(std::string("+"+part->name()));
                              }
                          }
                      }
                  }
              }
          }
        for (unsigned ibr=0; ibr < blocks.size(); ibr++)
          {
            make_parallel_consistent(eMesh,blocks[ibr]);
          }
      }

    if (0) std::cout << "tmp RefinerUtil::correctBlockNamesForPartPartConsistency_1: blocks = " << blocks << std::endl;
    return blocks;
  }


  //static
  void RefinerUtil::
  addAncestorsToUnrefineList(percept::PerceptMesh& eMesh, int num_levels_to_add, ElementUnrefineCollection& elements_to_unref)
  {
    int num_levels_to_add_1 = num_levels_to_add;
    if (num_levels_to_add < 0) num_levels_to_add_1 = 1000;
    for (int ilev=0; ilev < num_levels_to_add_1; ilev++)
      {
        ElementUnrefineCollection to_add(*eMesh.get_bulk_data());
        for (ElementUnrefineCollection::iterator iter = elements_to_unref.begin(); iter != elements_to_unref.end(); ++iter)
          {
            stk::mesh::Entity element = *iter;
            if (eMesh.hasFamilyTree(element))
              {
                stk::mesh::Entity parent = eMesh.getParent(element, false);
                if (elements_to_unref.find(parent) != elements_to_unref.end())
                  continue;
                std::vector<stk::mesh::Entity> children;
                bool hasChildren = eMesh.getChildren(parent, children, true, false);
                if (hasChildren && children.size())
                  {
                    bool allChildrenInUnrefSet = true;
                    for (unsigned ichild=0; ichild < children.size(); ichild++)
                      {
                        if (elements_to_unref.find(children[ichild]) == elements_to_unref.end())
                          {
                            allChildrenInUnrefSet = false;
                            break;
                          }
                      }
                    if (allChildrenInUnrefSet)
                      {
                        to_add.insert(parent);
                      }
                  }
              }
          }
        //std::cout << "RefinerUtil::addAncestorsToUnrefineList ilev= " << ilev << " to_add.size= " << to_add.size() << std::endl;
        if (to_add.size())
          {
            elements_to_unref.insert(to_add.begin(), to_add.end());
          }
        else
          {
            break;
          }
      }

  }

  /// create missing edges after adapt - for edge-based simulators
  void RefinerUtil::
  create_missing_edges(percept::PerceptMesh& eMesh)
  {
    typedef MySubDimCell<SDCEntityType, 2, CompareSDCEntityType> SubDimCell2;
    typedef boost::unordered_set<SubDimCell2, my_fast_hash<SDCEntityType, 2>, my_fast_equal_to<SDCEntityType, 2> > SubDimCellSet;
    SubDimCellSet edge_set, new_edge_set;

    // put existing edges in the set
    const stk::mesh::BucketVector & edge_buckets = eMesh.get_bulk_data()->buckets( eMesh.edge_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = edge_buckets.begin() ; k != edge_buckets.end() ; ++k )
      {
        stk::mesh::Bucket & bucket = **k ;
        // add all edges since we want to avoid creating them even if they are shared from another proc
        //if (bucket.owned())
        {
          const unsigned num_edges_in_bucket = bucket.size();

          for (unsigned i_edge = 0; i_edge < num_edges_in_bucket; i_edge++)
            {
              stk::mesh::Entity edge = bucket[i_edge];
              const percept::MyPairIterRelation edge_nodes (eMesh, edge, eMesh.node_rank());
              SubDimCell2 subDimEntity(eMesh);
              subDimEntity.clear();
              subDimEntity.insert(edge_nodes[0].entity());
              subDimEntity.insert(edge_nodes[1].entity());
              edge_set.insert(subDimEntity);
            }
        }
      }

    // visit elements (and sides if in 3D) and check if edges present; add as needed
    stk::mesh::EntityRank rank_start = eMesh.side_rank();
    if (eMesh.get_spatial_dim() == 2) rank_start = eMesh.element_rank();
    for (stk::mesh::EntityRank rank = rank_start; rank <= eMesh.element_rank(); ++rank)
      {
        std::string active_part_name = "refine_active_elements_part_"+toString(rank);
        //std::string inactive_part_name = "refine_inactive_elements_part_"+toString(rank);
        stk::mesh::Part* active_child_elements_part = eMesh.get_part(active_part_name);
        stk::mesh::Selector sel_part;
        if (active_child_elements_part)
          {
            sel_part = stk::mesh::Selector(*active_child_elements_part);
          }
        //stk::mesh::Part* inactive_parent_elements_part = eMesh.get_part(inactive_part_name);

        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            if (bucket.owned() && sel_part(bucket))
              {
                const CellTopologyData * cell_topo_data = eMesh.get_cell_topology(bucket);
                CellTopology topo(cell_topo_data);
                unsigned edge_count = cell_topo_data->edge_count;
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];
                    const percept::MyPairIterRelation elem_nodes (eMesh, element, eMesh.node_rank());
                    for (unsigned iedgeOrd=0; iedgeOrd < edge_count; iedgeOrd++)
                      {
                        unsigned in0 = cell_topo_data->edge[iedgeOrd].node[0];
                        unsigned in1 = cell_topo_data->edge[iedgeOrd].node[1];
                        SubDimCell2 subDimEntity(eMesh);
                        subDimEntity.clear();
                        subDimEntity.insert(elem_nodes[in0].entity());
                        subDimEntity.insert(elem_nodes[in1].entity());
                        SubDimCellSet::iterator fnd = edge_set.find(subDimEntity);
                        if (fnd == edge_set.end())
                          {
                            edge_set.insert(subDimEntity);
                            new_edge_set.insert(subDimEntity);
                          }
                      }
                  }
              }
          }
      }

    // now create new/missing edges
    eMesh.get_bulk_data()->modification_begin();
    std::vector<stk::mesh::Entity> new_edges;
    eMesh.createEntities(eMesh.edge_rank(), new_edge_set.size(), new_edges);

    unsigned count=0;
    for (SubDimCellSet::iterator edge_it = new_edge_set.begin(); edge_it != new_edge_set.end(); ++edge_it, ++count)
      {
        const SubDimCell2& edge = *edge_it;
        bool in_order = eMesh.identifier(edge[0]) < eMesh.identifier(edge[1]);
        if (in_order)
          {
            eMesh.get_bulk_data()->declare_relation(new_edges[count], edge[0], 0);
            eMesh.get_bulk_data()->declare_relation(new_edges[count], edge[1], 1);
          }
        else
          {
            eMesh.get_bulk_data()->declare_relation(new_edges[count], edge[1], 0);
            eMesh.get_bulk_data()->declare_relation(new_edges[count], edge[0], 1);
          }
      }
    stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
    eMesh.get_bulk_data()->modification_end();
  }

  void RefinerUtil::rebuild_family_tree(PerceptMesh& eMesh, bool debug)
  {
    if (!eMesh.m_parent_element_field_set)
      {
        eMesh.m_parent_element_field_set = true;
        eMesh.m_parent_element_field = eMesh.get_fem_meta_data()->get_field<ParentElementType>(eMesh.element_rank(), "parent_element");
        eMesh.m_parent_element_field_side = eMesh.get_fem_meta_data()->get_field<ParentElementType>(eMesh.side_rank(), "parent_element_side");
      }

    /// add nodes to new nodes part
    if (1)
      {
        if (!eMesh.m_new_nodes_field_set)
          {
            eMesh.m_new_nodes_field_set = true;
            eMesh.m_new_nodes_field = eMesh.get_fem_meta_data()->get_field<NewNodesType>(eMesh.node_rank(), "new_nodes");
          }

        eMesh.get_bulk_data()->modification_begin();
        stk::mesh::Part* new_nodes_part = eMesh.get_non_const_part("refine_new_nodes_part");
        VERIFY_OP_ON(new_nodes_part, != , 0, "new_nodes_part");
        if (new_nodes_part)
          {
            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.node_rank() );

            // first remove any existing nodes
            {
              std::vector<stk::mesh::Part*> remove_parts(1, new_nodes_part);
              std::vector<stk::mesh::Part*> add_parts;
              std::vector<stk::mesh::Entity> node_vec;

              stk::mesh::Selector removePartSelector(*new_nodes_part & eMesh.get_fem_meta_data()->locally_owned_part() );
              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;
                  if (removePartSelector(bucket))
                    {
                      const unsigned num_entity_in_bucket = bucket.size();
                      for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                        {
                          stk::mesh::Entity node = bucket[ientity];
                          node_vec.push_back(node);
                        }
                    }
                }
              for (unsigned ii=0; ii < node_vec.size(); ii++)
                {
                  eMesh.get_bulk_data()->change_entity_parts( node_vec[ii], add_parts, remove_parts );
                }
            }

            // now add new ones
            {
              std::vector<stk::mesh::Entity> new_nodes;
              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;
                  if (bucket.owned())
                    {
                      const unsigned num_entity_in_bucket = bucket.size();
                      for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                        {
                          stk::mesh::Entity node = bucket[ientity];
                          NewNodesType_type *ndata = stk::mesh::field_data(*eMesh.m_new_nodes_field, node);
                          if (ndata && ndata[0] != 0)
                            {
                              new_nodes.push_back(node);
                            }
                        }
                    }
                }

              std::vector<stk::mesh::Part*> add_parts(1, new_nodes_part);
              std::vector<stk::mesh::Part*> remove_parts;
              for (unsigned ind = 0; ind < new_nodes.size(); ind++)
                {
                  eMesh.get_bulk_data()->change_entity_parts( new_nodes[ind], add_parts, remove_parts );
                }
            }

            stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
            eMesh.get_bulk_data()->modification_end();
          }
      }

    //for (stk::mesh::EntityRank rank = eMesh.element_rank(); rank <= eMesh.element_rank(); ++rank)
    for (stk::mesh::EntityRank rank = eMesh.side_rank(); rank <= eMesh.element_rank(); ++rank)
      {
        std::ostringstream orank;
        orank << rank;
        std::string srank = orank.str();
        //int irank = (int)rank;
        std::vector<stk::mesh::Entity> ft_new_elements;
        const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);

        size_t num_elem_needed = 0;
        const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( rank );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();
            if (bucket.owned())
              {
                num_elem_needed += num_elements_in_bucket;
              }
          }

        eMesh.get_bulk_data()->modification_begin();
        if (debug)
          std::cout << "RefinerUtil::rebuild_family_tree_child: for rank= " << rank << " num_elem_needed= " << num_elem_needed << std::endl;
        eMesh.createEntities( FAMILY_TREE_RANK, 2*num_elem_needed, ft_new_elements);
        size_t i_ft = 0;

        std::vector<stk::mesh::Entity> elements;
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket = bucket.size();
            if (bucket.owned())
              {
                for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                  {
                    stk::mesh::Entity element = bucket[iElement];
                    elements.push_back(element);
                  }
              }
          }

        for (size_t ii=0; ii < elements.size(); ++ii)
          {
            stk::mesh::Entity element = elements[ii];
            ParentElementType_type *parent_element_id_from_field = NULL;

            if (is_matching_rank(*eMesh.m_parent_element_field, element))
              {
                parent_element_id_from_field = stk::mesh::field_data( *eMesh.m_parent_element_field , element );
              }
            else if (eMesh.m_parent_element_field_side && is_matching_rank(*eMesh.m_parent_element_field_side, element))
              {
                parent_element_id_from_field = stk::mesh::field_data( *eMesh.m_parent_element_field_side , element );
              }
            VERIFY_OP_ON(parent_element_id_from_field, !=, 0, "hmmm");
            stk::mesh::EntityId parent_element_id = static_cast<stk::mesh::EntityId>(parent_element_id_from_field[0]);
            stk::mesh::Entity parent_elem = stk::mesh::Entity();
            if (parent_element_id) parent_elem = eMesh.get_bulk_data()->get_entity(rank, parent_element_id);
            if (debug && rank == eMesh.side_rank())
              {
                std::cout << "RefinerUtil::rebuild_family_tree_child: side= " << eMesh.identifier(element) << " parent_element_id_from_field= "
                          << parent_element_id << " parent.valid= " << eMesh.is_valid(parent_elem);
                eMesh.print(element);
              }

            if (parent_element_id)
              {

                if (rank == eMesh.side_rank())
                  {
                    //VERIFY_OP_ON(eMesh.is_valid(parent_elem), ==, false, "root cause");
                    parent_elem = stk::mesh::Entity();

                    stk::mesh::EntityId id_new = 0;
                    stk::mesh::ConnectivityOrdinal ord;
                    eMesh.decipher_exodus_side_id(parent_element_id, id_new, ord);
                    VERIFY_OP_ON(id_new, !=, 0, "bad id_new");
                    if (id_new)
                      {
                        stk::mesh::Entity element_owner = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), id_new);
                        VERIFY_OP_ON(eMesh.is_valid(element_owner), ==, true, "bad parent side elem");
                        percept::MyPairIterRelation elem_to_side_rels (eMesh, element_owner, rank);
                        for (unsigned jj=0; jj < elem_to_side_rels.size(); ++jj)
                          {
                            if (elem_to_side_rels[jj].relation_ordinal() == ord)
                              {
                                parent_elem = elem_to_side_rels[jj].entity();
                                break;
                              }
                          }
                      }
                  }
                VERIFY_OP_ON(eMesh.is_valid(parent_elem), ==, true, "bad parent found, rank= "+srank);
                VERIFY_OP_ON(parent_elem, !=, element, "hmmm, parent=element");
                if (1)
                  {
                    stk::mesh::Entity family_tree = stk::mesh::Entity();
                    bool found = false;
                    unsigned ordinal = 0;
                    if (eMesh.hasFamilyTree(parent_elem))
                      {
                        //percept::MyPairIterRelation family_tree_relations (eMesh, element, FAMILY_TREE_RANK);
                        percept::MyPairIterRelation parent_to_family_tree_relations (eMesh, parent_elem, FAMILY_TREE_RANK);
                        VERIFY_OP_ON(parent_to_family_tree_relations.size(), > , 0, "hmmm");
                        VERIFY_OP_ON(parent_to_family_tree_relations.size(), <=, 2, "bad family tree relations");

                        bool isParentAlready = false;
                        if (parent_to_family_tree_relations.size() == 1)
                          {
                            unsigned parent_elem_ft_level_0 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, parent_elem);
                            family_tree = parent_to_family_tree_relations[parent_elem_ft_level_0].entity();
                          }
                        else if (parent_to_family_tree_relations.size() == 2)
                          {
                            unsigned parent_elem_ft_level_1 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, parent_elem);
                            family_tree = parent_to_family_tree_relations[parent_elem_ft_level_1].entity();
                            isParentAlready = true;
                          }

                        VERIFY_OP_ON(eMesh.is_valid(family_tree), ==, true, "hmmm");
                        percept::MyPairIterRelation family_tree_relations (eMesh, family_tree, eMesh.entity_rank(parent_elem));
                        VERIFY_OP_ON(family_tree_relations.size(), > , 0, "hmmm");
                        if (family_tree_relations[0].entity() == parent_elem)
                          {
                            isParentAlready = true;
                            VERIFY_OP_ON(eMesh.isParentElement(parent_elem), ==, true, "hmmm");
                          }
                        for (unsigned i = 1; i < family_tree_relations.size(); i++)
                          {
                            //if (family_tree_relations[i].relation_ordinal() == (ordinal + 1))
                            if (family_tree_relations[i].entity() == element)
                              {
                                found = true;
                                break;
                              }
                          }
                        if (!found && isParentAlready)
                          {
                            VERIFY_OP_ON(family_tree_relations.size(), > , 0, "hmmm");
                            ordinal = family_tree_relations.size() - 1;
                          }
                      }
                    if (!found)
                      {
                        if (i_ft >= ft_new_elements.size())
                          {
                            throw std::runtime_error("ran out of ft_new_elements");
                          }
                        stk::mesh::Entity familyTreeNewElement = ft_new_elements[i_ft++];
                        UniformRefinerPatternBase::set_parent_child_relations(eMesh, parent_elem, element, familyTreeNewElement, ordinal);
                      }
                  }
              }
          }

        std::vector<stk::mesh::Entity> ft_delete;
        for (unsigned ii=i_ft; ii < ft_new_elements.size(); ++ii)
          {
            if (eMesh.get_bulk_data()->has_no_relations(ft_new_elements[ii]))
              ft_delete.push_back(ft_new_elements[ii]);
          }
        for (unsigned ii=0; ii < ft_delete.size(); ++ii)
          {
            if (!eMesh.get_bulk_data()->destroy_entity(ft_delete[ii]))
              throw std::runtime_error("bad destroy of family_tree");
          }
        stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
        eMesh.get_bulk_data()->modification_end();

      }

    // get parts correct
    eMesh.get_bulk_data()->modification_begin();
    Refiner::set_active_part(eMesh);
    stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
    eMesh.get_bulk_data()->modification_end();

  }

  std::string printSubDim(PerceptMesh& eMesh, const SubDimCell_SDCEntityType& c)
  {
    std::ostringstream out;
    out << "SDC sz= " << c.size() << " {";
    for (unsigned i = 0; i < c.size(); i++)
      {
        out << eMesh.id(c[i]) << " ";
      }
    out << "}";
    return out.str();
  }

  void RefinerUtil::save_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, const std::string& msg, bool doComm)
  {
    if (eMesh.m_node_registry_field == 0)
      return;

    // if (eMesh.get_rank() == 0)
    //   std::cout << "save_node_registry start eMesh.file= " << eMesh.getProperty("in_filename") << std::endl;

    bool debug = false;

    if (1)
      {
        std::vector<stk::mesh::Entity> vec;
        //stk::mesh::Selector sel = eMesh.get_fem_meta_data()->locally_owned_part() | eMesh.get_fem_meta_data()->globally_shared_part();
        stk::mesh::Selector sel = eMesh.get_fem_meta_data()->universal_part();
        stk::mesh::get_selected_entities(sel , eMesh.get_bulk_data()->buckets(eMesh.node_rank()), vec);

        SubDimCell_SDCEntityType subDimEntity(eMesh);

        // FIXME
        for (size_t ii=0; ii < vec.size(); ++ii)
          {
            stk::mesh::Entity node = vec[ii];

            //debug = eMesh.id(node) == 148 || eMesh.id(node) == 2568;
            double *node_data = stk::mesh::field_data(*eMesh.m_node_registry_field, node);
            for (unsigned kk=0; kk < NUM_NR_FIELD_SLOTS; ++kk)
              {
                node_data[kk] = 0.0;
              }
          }
      }

    SubDimCellToDataMap& map = nodeRegistry.getMap();
    SubDimCellToDataMap::iterator iter;

    for (iter = map.begin(); iter != map.end(); ++iter)
      {
        const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
        SubDimCellData& nodeId_elementOwnderId = (*iter).second;

        // tuple storage: SDC_DATA_GLOBAL_NODE_IDS, SDC_DATA_OWNING_ELEMENT_KEY,  SDC_DATA_OWNING_SUBDIM_RANK, SDC_DATA_OWNING_SUBDIM_ORDINAL, SDC_DATA_SPACING
        NodeIdsOnSubDimEntityType& nodeIds_onSE      = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        //debug = eMesh.is_valid(nodeIds_onSE[0]) && eMesh.id(nodeIds_onSE[0]) == 38870;
        if (debug && nodeIds_onSE.size() && eMesh.is_valid(nodeIds_onSE[0]))
          std::cout << eMesh.rank() << " node= " << eMesh.print_entity_compact(nodeIds_onSE[0]) << std::endl;

        stk::mesh::EntityKey       owningElementKey  = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>();
        stk::mesh::EntityId        owningElementId   = owningElementKey.id();
        (void)owningElementId;
        stk::mesh::EntityRank      owningElementRank = owningElementKey.rank();
        //VERIFY_OP_ON(owningElementRank, ==, eMesh.element_rank(), "bad owningElementRank");
        unsigned               owningSubDimRank  = nodeId_elementOwnderId.get<SDC_DATA_OWNING_SUBDIM_RANK>();
        unsigned               owningSubDimOrd   = nodeId_elementOwnderId.get<SDC_DATA_OWNING_SUBDIM_ORDINAL>();
        VERIFY_OP_ON(owningSubDimOrd, >, 0, "hmm 2");
        --owningSubDimOrd;
        unsigned owningSubDimSize = subDimEntity.size();
        (void)owningSubDimSize;

        //Double3                    sdcSpacing        = nodeId_elementOwnderId.get<SDC_DATA_SPACING>();

        stk::mesh::Entity owningElement = eMesh.get_bulk_data()->get_entity(owningElementKey.rank(), owningElementKey.id());
        if (!eMesh.is_valid(owningElement))
          {
            if (debug) std::cout << eMesh.rank() << "save_node_registry invalid owningElement= " << owningElementKey << std::endl;
            continue;
          }

        if (debug)
          {
            SubDimCell_SDCEntityType subDimEntity1(eMesh);
            nodeRegistry.getSubDimEntity(subDimEntity1, owningElement, static_cast<stk::mesh::EntityRank>(owningSubDimRank), static_cast<unsigned>(owningSubDimOrd));
            for (unsigned ii=0; ii < subDimEntity.size(); ++ii)
              {
                std::cout << "save_node_registry ii= " << ii << " ndi= " << eMesh.id(subDimEntity1[ii]) << " sdi= " << eMesh.id(subDimEntity[ii]) << std::endl;
                if (0) VERIFY_OP_ON(eMesh.id(subDimEntity1[ii]), ==,  eMesh.id(subDimEntity[ii]) , "save_node_registry bad ndi");
              }
          }

        for (unsigned jj = 0; jj < nodeIds_onSE.size(); ++jj)
          {
            stk::mesh::Entity node = nodeIds_onSE[jj];

            //debug = eMesh.id(node) == 148 || eMesh.id(node) == 2568;

            if (!eMesh.is_valid(node))
              {
                if (debug) std::cout << "save_node_registry:: node invalid" << std::endl;
                continue;
              }
            NodeRegistryFieldType_type *node_data = stk::mesh::field_data(*eMesh.m_node_registry_field, node);
            if (!node_data)
              {
                if (debug) std::cout << "save_node_registry:: node_data invalid, id= " << eMesh.id(node) << std::endl;
                continue;
              }

            // FIXME
#if 0
            if (node_data[NR_FIELD_OWNING_ELEMENT_ID] != 0.0)
              {
#define PP(x) " " << #x << "= " << x
                std::cout << eMesh.rank() << "node data duplicated, existing/new: msg= " << msg
                          << "\n" << PP(node_data[NR_FIELD_OWNING_ELEMENT_ID]) << PP(owningElementKey.id())
                          << "\n" << PP(node_data[NR_FIELD_OWNING_ELEMENT_RANK]) << PP(owningElementRank)
                          << "\n" << PP(node_data[NR_FIELD_MARK]) << PP(nodeIds_onSE.m_mark)
                          << "\n" << PP(node_data[NR_FIELD_OWNING_SUBDIM_RANK]) << PP(owningSubDimRank)
                          << "\n" << PP(node_data[NR_FIELD_OWNING_SUBDIM_ORDINAL]) << PP(owningSubDimOrd + 1)
                          << "\n" << PP(node_data[NR_FIELD_OWNING_SUBDIM_SIZE]) << PP(owningSubDimSize)
                          << std::endl;
                std::cout << eMesh.demangled_stacktrace(20)
                          << std::endl;
                if (1)
                  {
                    SubDimCell_SDCEntityType subDimEntity1(eMesh);

                    stk::mesh::EntityId  owningElementId1 (  static_cast<stk::mesh::EntityId>(node_data[NR_FIELD_OWNING_ELEMENT_ID]) );
                    stk::mesh::EntityRank  owningElementRank1 (  static_cast<stk::mesh::EntityRank>(node_data[NR_FIELD_OWNING_ELEMENT_RANK]) );
                    stk::mesh::EntityKey owningElementKey1(owningElementRank1, owningElementId1);
                    stk::mesh::Entity owningElement1 = eMesh.get_bulk_data()->get_entity(owningElementRank1, owningElementId1);
                    //VERIFY_OP_ON(eMesh.is_valid(owningElement), ==, true, "bad owningElement");
                    stk::mesh::EntityRank owningSubDimRank1 = static_cast<stk::mesh::EntityRank>(node_data[NR_FIELD_OWNING_SUBDIM_RANK]);
                    unsigned              owningSubDimOrd1  = static_cast<unsigned>(node_data[NR_FIELD_OWNING_SUBDIM_ORDINAL]);
                    VERIFY_OP_ON(owningSubDimOrd1, >, 0, "hmmm 3");
                    --owningSubDimOrd1;
                    unsigned              owningSubDimSize1 = static_cast<unsigned>(node_data[NR_FIELD_OWNING_SUBDIM_SIZE]);

                    bool foundGhostNode1 = nodeRegistry.getSubDimEntity(subDimEntity1, owningElement1, owningSubDimRank1, owningSubDimOrd1);
                    if (foundGhostNode1)
                      std::cout << eMesh.rank() << " rebuild_node_registry foundGhostNode1= " << foundGhostNode1 << std::endl;
                    std::cout << eMesh.rank() << " msg= " << msg
                              << "\n" << PP(owningElementId) << " " << PP(owningElementId1)
                              << "\n" << PP(owningElementRank) << " " << PP(owningElementRank1)
                              << "\n" << PP(owningSubDimRank) << " " << PP(owningSubDimRank1)
                              << "\n" << PP(owningSubDimOrd) << " " << PP(owningSubDimOrd1)
                              << "\n" << PP(owningSubDimSize) << " " << PP(owningSubDimSize1)
                              << "\n" << "subDimEntity= " << printSubDim(eMesh, subDimEntity)
                              << "\n" << "subDimEntity1= " << printSubDim(eMesh, subDimEntity1)
                              << std::endl;
                  }
              }

            VERIFY_OP_ON(node_data[NR_FIELD_OWNING_ELEMENT_ID], ==, 0.0, "node data is duplicated");
#endif

            node_data[NR_FIELD_OWNING_ELEMENT_ID]     = static_cast<NodeRegistryFieldType_type>(owningElementKey.id());
            node_data[NR_FIELD_OWNING_ELEMENT_RANK]   = static_cast<NodeRegistryFieldType_type>(owningElementRank);
            node_data[NR_FIELD_MARK]                  = static_cast<NodeRegistryFieldType_type>(nodeIds_onSE.m_mark);
            node_data[NR_FIELD_OWNING_SUBDIM_RANK]    = static_cast<NodeRegistryFieldType_type>(owningSubDimRank);
            node_data[NR_FIELD_OWNING_SUBDIM_ORDINAL] = static_cast<NodeRegistryFieldType_type>(owningSubDimOrd + 1);
            node_data[NR_FIELD_OWNING_SUBDIM_SIZE]    = static_cast<NodeRegistryFieldType_type>(owningSubDimSize);

            if (debug) {
              std::cout << eMesh.rank() << " SNR:: [OE: " << owningElementKey.rank() << " " << owningElementKey.id() << " SDR/I: " << (int)owningSubDimRank << " " << (int)owningSubDimOrd << " NDS: " << eMesh.id(node)
                        << "\n" << eMesh.rank() << " owningElement= " << eMesh.print_entity_compact(owningElement)
                        << std::endl;
            }

          }
      }

    if (doComm)
      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(eMesh.m_node_registry_field);
        stk::mesh::copy_owned_to_shared(*eMesh.get_bulk_data(), fields);
        //stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);
      }
  }

  // rebuild
  void RefinerUtil::rebuild_node_registry(PerceptMesh& eMesh, NodeRegistry& nodeRegistry, bool initNR, PerceptMesh* eMeshNR, NodeRegistry *compareNR, bool skipEmpty)
  {
    if (eMesh.m_node_registry_field == 0)
      return;

    // if (eMesh.get_rank() == 0)
    //   std::cout << "rebuild_node_registry start eMesh.file= " << eMesh.getProperty("in_filename") << std::endl;

    if (initNR)
      {
        nodeRegistry.initialize();
        nodeRegistry.getMap().clear();
        nodeRegistry.init_comm_all();
      }

    if (1)
      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(eMesh.m_node_registry_field);
        stk::mesh::copy_owned_to_shared(*eMesh.get_bulk_data(), fields);
        //stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), fields);
      }

    std::vector<stk::mesh::Entity> vec;
    //stk::mesh::Selector sel = eMesh.get_fem_meta_data()->locally_owned_part() | eMesh.get_fem_meta_data()->globally_shared_part();
    stk::mesh::Selector sel = eMesh.get_fem_meta_data()->universal_part();
    stk::mesh::get_selected_entities(sel , eMesh.get_bulk_data()->buckets(eMesh.node_rank()), vec);

    SubDimCell_SDCEntityType subDimEntity(eMesh);

    for (size_t inode=0; inode < vec.size(); ++inode)
      {
        stk::mesh::Entity node = vec[inode];

        bool debug = false; // eMesh.id(node) == 38870; // || eMesh.id(node) == 2568;
        // if (eMesh.id(node) == 69776 && eMesh.get_rank() == 5)
        //   {
        //     debug = true;
        //   }
        // if (eMesh.aura(node))
        //   continue;
        double *node_data = stk::mesh::field_data(*eMesh.m_node_registry_field, node);

        stk::mesh::EntityId  owningElementId (  static_cast<stk::mesh::EntityId>(node_data[NR_FIELD_OWNING_ELEMENT_ID]) );
        stk::mesh::EntityRank  owningElementRank (  static_cast<stk::mesh::EntityRank>(node_data[NR_FIELD_OWNING_ELEMENT_RANK]) );
        stk::mesh::EntityKey owningElementKey(owningElementRank, owningElementId);
        stk::mesh::Entity owningElement = eMesh.get_bulk_data()->get_entity(owningElementRank, owningElementId);
        if (!eMesh.is_valid(owningElement))
          {
            if (debug) std::cout << eMesh.rank() << " rebuild_node_registry:: bad owningElement node= " << eMesh.id(node) << " owningElement= " << owningElementKey << std::endl;
            continue;
          }
        VERIFY_OP_ON(eMesh.is_valid(owningElement), ==, true, "bad owningElement");
        stk::mesh::EntityRank owningSubDimRank = static_cast<stk::mesh::EntityRank>(node_data[NR_FIELD_OWNING_SUBDIM_RANK]);
        unsigned              owningSubDimOrd  = static_cast<unsigned>(node_data[NR_FIELD_OWNING_SUBDIM_ORDINAL]);
        VERIFY_OP_ON(owningSubDimOrd, >, 0, "hmmm 3");
        --owningSubDimOrd;
        unsigned              owningSubDimSize = static_cast<unsigned>(node_data[NR_FIELD_OWNING_SUBDIM_SIZE]);

        bool foundGhostNode = nodeRegistry.getSubDimEntity(subDimEntity, owningElement, owningSubDimRank, owningSubDimOrd);
        if (foundGhostNode && nodeRegistry.getCheckForGhostedNodes())
          {
            //std::cout << eMesh.rank() << " rebuild_node_registry foundGhostNode= " << foundGhostNode << std::endl;
            continue;
          }
        VERIFY_OP_ON(subDimEntity.size(), ==, owningSubDimSize, "bad owningSubDimSize");

        if (0)
        for (unsigned ii=0; ii < subDimEntity.size(); ++ii)
          {
            stk::mesh::EntityId ndi = static_cast<stk::mesh::EntityId>(node_data[ii]);
            std::cout << "ii= " << ii << " ndi= " << ndi << " sdi= " << eMesh.id(subDimEntity[ii]) << std::endl;
          }

        for (unsigned ii=0; ii < subDimEntity.size(); ++ii)
          {
            VERIFY_OP_ON(eMesh.is_valid(subDimEntity[ii]), ==, true, "bad node in rebuild_node_registry");
            stk::mesh::EntityId ndi = static_cast<stk::mesh::EntityId>(node_data[ii]);
            if (0) VERIFY_OP_ON(eMesh.id(subDimEntity[ii]), ==, ndi, "bad node data");
          }
        static SubDimCellData empty_SubDimCellData;
        SubDimCellData* nodeId_elementOwnderId_ptr = nodeRegistry.getFromMapPtr(subDimEntity);
        SubDimCellData& nodeId_elementOwnderId = (nodeId_elementOwnderId_ptr ? *nodeId_elementOwnderId_ptr : empty_SubDimCellData);
        bool is_empty = nodeId_elementOwnderId_ptr == 0;
        NodeIdsOnSubDimEntityType& nodeIds_onSE1 = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();

        if (debug) {
          std::cout << eMesh.rank() << "tmp srk  is_empty= " << is_empty << " nodeIds_onSE1.size= " << nodeIds_onSE1.size()  << " " << (nodeIds_onSE1.size() ? nodeIds_onSE1[0].m_value : 1234567ll)
                    << std::endl;
          std::cout << "nid: " << ((nodeIds_onSE1.size() && eMesh.is_valid(nodeIds_onSE1[0]) )? eMesh.print_entity_compact(nodeIds_onSE1[0]) : " none") << std::endl;
        }
        if (is_empty)
          {
            SubDimCellData data( NodeIdsOnSubDimEntityType(1, stk::mesh::Entity(), nodeRegistry.NR_MARK_NONE),
                                 owningElementKey, owningSubDimRank, owningSubDimOrd + 1);
            NodeIdsOnSubDimEntityType& nid_new = data.get<SDC_DATA_GLOBAL_NODE_IDS>();
            nid_new.resize(1);
            nid_new[0] = node;
            nid_new.m_entity_id_vector[0] = eMesh.id(node);
            nid_new.m_mark = static_cast<unsigned>(node_data[NR_FIELD_MARK]);

            if (debug) std::cout << eMesh.rank() << "tmp srk  is_empty= " << is_empty << " node= " << eMesh.print_entity_compact(node) << std::endl;

            // Double3& spc = data.get<SDC_DATA_SPACING>();
            // spc[0] = node_data[NR_FIELD_OWN_SPACING_0];
            // spc[1] = node_data[NR_FIELD_OWN_SPACING_1];
            nodeRegistry.putInMap(subDimEntity,  data);
          }
        else
          {
            NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
            nodeIds_onSE.push_back(node);
            nodeIds_onSE.m_entity_id_vector.push_back(eMesh.id(node));
            nodeIds_onSE.m_mark = static_cast<unsigned>(node_data[NR_FIELD_MARK]);
            // Double2& spc = nodeId_elementOwnderId.get<SDC_DATA_SPACING>();
            // spc[0] = node_data[NR_FIELD_OWN_SPACING_0];
            // spc[1] = node_data[NR_FIELD_OWN_SPACING_1];
            std::ostringstream out;
            if (debug)
              {
                for (unsigned ij=0; ij < nodeIds_onSE.size(); ++ij)
                  {
                    stk::mesh::Entity n1 = nodeIds_onSE[ij];
                    double *node_data_1 = stk::mesh::field_data(*eMesh.m_node_registry_field, n1);

                    out << eMesh.rank() << " ij= " << ij << " nid= " << eMesh.print_entity_compact(nodeIds_onSE[ij])
                        << " nd= ";
                    for (unsigned ik=0; ik < NUM_NR_FIELD_SLOTS; ++ik)
                      {
                        out << " " << (int64_t)node_data_1[ik];
                      }
                    out << std::endl;
                  }

                out << eMesh.rank()
                    << "subDimEntity= " << printSubDim(eMesh, subDimEntity)
                    << std::endl;
                std::cout << out.str();
              }
            VERIFY_OP_ON(nodeIds_onSE.size(), ==, 1, "bad size on proc: "+eMesh.rank());
          }
      }

    //nodeRegistry.checkDB("after rebuild_node_registry");
    if (compareNR && eMeshNR)
      {
        compare(eMesh, nodeRegistry, *eMeshNR, *compareNR, skipEmpty);
      }
  }

  static void write(std::ofstream& out , NodeRegistry& nodeRegistry, bool skipEmpty = true)
  {
    SubDimCellToDataMap::iterator iter;
    SubDimCellToDataMap& map = nodeRegistry.getMap();
    for (iter = map.begin(); iter != map.end(); ++iter)
      {
        const SubDimCell_SDCEntityType& subDimEntity = (*iter).first;
        SubDimCellData& nodeId_elementOwnderId = (*iter).second;
        NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
        if (skipEmpty && nodeIds_onSE.size() == 0)
          continue;

        out << "[SDE: ";
        for (unsigned k=0; k < subDimEntity.size(); k++)
          {
            out << " " << nodeRegistry.getMesh().id(subDimEntity[k]) << " ";
            //emitter << nodeRegistry.getMesh().identifier(subDimEntity[k]);          YAML_ERRCHECK;
          }

        stk::mesh::EntityKey& owningElementKey = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>();
        unsigned char         owningSubDimRank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_SUBDIM_RANK>();
        unsigned char         owningSubDimOrd  = nodeId_elementOwnderId.get<SDC_DATA_OWNING_SUBDIM_ORDINAL>();
        unsigned              owningSubDimSize = subDimEntity.size();
        out << "]: [OE: " << owningElementKey.rank() << " " << owningElementKey.id() << " SDR/I: " << (int)owningSubDimRank << " " << (int)owningSubDimOrd << " NDS: " << " sz= " << owningSubDimSize;
        for (unsigned ii = 0; ii < nodeIds_onSE.size(); ii++)
          {
            //emitter << (int)nodeIds_onSE[ii]->identifier();      YAML_ERRCHECK;
            stk::mesh::EntityId id = nodeIds_onSE.m_entity_id_vector[ii];
            stk::mesh::EntityId id_check = nodeRegistry.getMesh().identifier(nodeIds_onSE[ii]);
            VERIFY_OP_ON(id_check, ==, id, "write id");

            out << " " << (int)nodeIds_onSE.m_entity_id_vector[ii];
          }
        out << "]\n";
      }
  }

  void RefinerUtil::compare(PerceptMesh& eMesh0, NodeRegistry& nr0, PerceptMesh& eMesh1, NodeRegistry& nr1, bool skipEmpty)
  {
    std::ostringstream msg;
    // if (nr0.getMap().size() != nr1.getMap().size())
    //   msg << "sizes different: nr0= " << nr0.getMap().size() << " nr1= " << nr1.getMap().size() << std::endl;

    for (SubDimCellToDataMap::iterator iter0 = nr0.getMap().begin(); iter0 != nr0.getMap().end(); ++iter0)
      {
        const SubDimCell_SDCEntityType& subDimEntity0 = (*iter0).first;
        SubDimCellData&            subDimCellData0      = (*iter0).second;
        NodeIdsOnSubDimEntityType& nodeIds_onSE0 = subDimCellData0.get<SDC_DATA_GLOBAL_NODE_IDS>();
        if (skipEmpty && nodeIds_onSE0.size() == 0)
          continue;

        SubDimCellToDataMap::iterator iter1 = nr1.getMap().find(subDimEntity0);
        if (iter1 == nr1.getMap().end())
          {
            msg << " couldn't find subDimEntity of 0 in nr1: " << subDimEntity0 << std::endl;
            continue;
          }

        {
          SubDimCellData&            subDimCellData1       = (*iter1).second;
          const SubDimCell_SDCEntityType& subDimEntity1 = iter1->first;
          if (subDimEntity0.size() != subDimEntity1.size())
            {
              msg << "\nsubDimEntity sizes not equal: " << subDimEntity0.size() << " " << subDimEntity1.size();
            }
          NodeIdsOnSubDimEntityType& nodeIds_onSE1 = subDimCellData1.get<SDC_DATA_GLOBAL_NODE_IDS>();
          bool nidEq = nodeIds_onSE1.size() == nodeIds_onSE0.size();
          if (!nidEq)
            msg << " sizes not equal: nr0.size,1 ="  << nodeIds_onSE0.size() << " " << nodeIds_onSE1.size() << std::endl;

          if (nodeIds_onSE0.m_mark != nodeIds_onSE1.m_mark)
            msg << " marks not equal: " << nodeIds_onSE0.m_mark << " " << nodeIds_onSE1.m_mark;

          if (nidEq)
            for (unsigned ii=0; ii < nodeIds_onSE0.size(); ++ii)
              {
                if (nodeIds_onSE0[ii] != nodeIds_onSE1[ii]
                    || nodeIds_onSE0.m_entity_id_vector[ii] != nodeIds_onSE1.m_entity_id_vector[ii])
                  {
                    msg << "nid not equal = " << eMesh0.id(nodeIds_onSE0[ii]) << " " << eMesh1.id(nodeIds_onSE1[ii])
                        << " " << nodeIds_onSE0.m_entity_id_vector[ii] << " " << nodeIds_onSE1.m_entity_id_vector[ii]
                        << std::endl;
                    nidEq = false;
                  }
              }

          if (!nidEq) {
            msg << " data not equal: nr0,1 =\n"  << subDimCellData0
                << "\n" << subDimCellData1
                << "\nkey0= sz= " << subDimEntity0.size();
            for (unsigned i0=0; i0 < subDimEntity0.size(); ++i0)
              msg << " " << subDimEntity0[i0];

            msg << "\nkey1= sz= " << subDimEntity1.size();
            for (unsigned i1=0; i1 < subDimEntity1.size(); ++i1)
              msg << " " << subDimEntity1[i1];
            msg << std::endl;
          }
          bool success = true;

#define COMP(a) do { if (subDimCellData0.get<a>() != subDimCellData1.get<a>()) {\
              if (0) std::cout << eMesh0.rank() << "\nmsg=" << msg.str() << "\nbad " << std::string( #a ) << " node= " << eMesh0.id(nodeIds_onSE0[0]) << std::endl; } \
            VERIFY_OP_ON_BOOL_NO_FAIL((int)subDimCellData0.get<a>(), ==, (int)subDimCellData1.get<a>(), "bad: " + std::string( #a ), success); \
          } while (0)
#define COMP1(a,x) do { if (subDimCellData0.get<a>().x() != subDimCellData1.get<a>().x()) { \
              if (0) std::cout << eMesh0.rank() << "\nmsg=" << msg.str() << "\nbad " << std::string( #a )+std::string( #x ) << " node= " << eMesh0.id(nodeIds_onSE0[0]) << std::endl; } \
            VERIFY_OP_ON_BOOL_NO_FAIL(subDimCellData0.get<a>().x(), ==, subDimCellData1.get<a>().x(), "bad: " + std::string( #a ), success); \
          } while (0)

          COMP1(SDC_DATA_OWNING_ELEMENT_KEY,rank);
          COMP1(SDC_DATA_OWNING_ELEMENT_KEY,id);
          COMP(SDC_DATA_OWNING_SUBDIM_RANK);
          COMP(SDC_DATA_OWNING_SUBDIM_ORDINAL);
          //COMP(SDC_DATA_SPACING);
#undef COMP
          if (!success)
            {
              if (subDimEntity0.size() == subDimEntity1.size())
                {
                  for (unsigned i0=0; i0 < subDimEntity0.size(); ++i0)
                    {
                      if (subDimEntity0[i0] != subDimEntity1[i0])
                        {
                          msg << " failed comparing data and subDimEntity's not equal: subDimEntity0 = ";
                          for (unsigned j0=0; j0 < subDimEntity0.size(); ++j0)
                            msg << " " << subDimEntity0[j0] << " " << subDimEntity1[j0] << std::endl;
                        }
                    }
                }
              else
                {
                  msg << " failed comparing data";
                }
            }
        }
      }

    if (msg.str().size())
      {
        std::cout << "compare is false: " << msg.str() << std::endl;

        if (1)
          {
            YAML::Emitter yaml0, yaml1;
            SerializeNodeRegistry::serialize_write(nr0, yaml0, 0);
            SerializeNodeRegistry::serialize_write(nr1, yaml1, 0);
            static int iter = 0;
            std::ofstream y0("cy0_"+toString(iter)+".yaml");
            std::ofstream y1("cy1_"+toString(iter)+".yaml");
            ++iter;
            y0 << "# " << eMesh0.getProperty("in_filename") << std::endl;
            y1 << "# " << eMesh1.getProperty("in_filename") << std::endl;

            y0 << yaml0.c_str() << std::endl;
            y1 << yaml1.c_str() << std::endl;

          }
        if (1)
          {
            static int iter = 0;
            std::ofstream y0("dy0_"+toString(iter)+".yaml");
            std::ofstream y1("dy1_"+toString(iter)+".yaml");
            ++iter;
            y0 << "# " << eMesh0.getProperty("in_filename") << std::endl;
            y1 << "# " << eMesh1.getProperty("in_filename") << std::endl;

            write(y0, nr0, skipEmpty);
            write(y1, nr1, skipEmpty);

          }

        VERIFY_MSG("compare");
      }
  }

}

