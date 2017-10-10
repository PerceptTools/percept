// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#include <adapt/main/MeshAdapt.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include <adapt/AdaptHelperFunctions.hpp>
#include <adapt/FixSideSets.hpp>

#include <stk_util/registry/ProductRegistry.hpp>

#if defined(STK_BUILT_IN_SIERRA)
#include <AuditLog.h>
#endif

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <percept/mesh/geometry/kernel/GeometryKernelGregoryPatch.hpp>
#include <percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>

#ifdef HAVE_CUBIT
#include "PGeom.hpp"
#include "PGeomAssoc.hpp"
#ifdef HAVE_ACIS
#include "PGeomACIS.hpp"
#endif

#include <percept/mesh/geometry/kernel/GeometryKernelPGEOM.hpp>
#endif

#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <percept/mesh/geometry/kernel/GeometryFactory.hpp>

#endif


namespace percept {

  MeshAdapt::MeshAdapt() 
  {



#if !MESH_ADAPT_CPP11_MVI
#define CPP11MVITYPE(x)
#define CPP11MVI(x) = x
#define CPP11MVI_3(x,y,z)
#define CPP11MVIEQ(x) x
#define CPP11MV1_ARRAY_INIT(type, var, sz, tmpvar, arrInit) do {type tmpvar[sz] = arrInit; for (unsigned itmp=0; itmp < sz; itmp++) { var[itmp] = tmpvar[itmp]; } } while(0)

#include "MeshAdaptMemberVarInit.hpp"

#undef CPP11MVITYPE
#undef CPP11MVI
#undef CPP11MVI_3
#undef CPP11MVIEQ
#undef CPP11MV1_ARRAY_INIT



#endif

  }

  MemorySizeType MeshAdapt::memory_dump(int dump_level, const stk::ParallelMachine& comm, stk::mesh::BulkData& bulkData, NodeRegistry* node_reg, std::string msg)
  {
    MemorySizeType returned_total_memory = 0;
    MemorySizeType rss_current_0 = 0;
    MemorySizeType rss_high_water_mark_0 = 0;
    stk::get_memory_usage(rss_current_0, rss_high_water_mark_0);

    const MemorySizeType MB = 1024*1024;

    stk::all_reduce( comm, stk::ReduceSum<1>( &rss_current_0 ) );
    stk::all_reduce( comm, stk::ReduceSum<1>( &rss_high_water_mark_0 ) );


    MemorySizeType node_reg_mem = 0;
    if (node_reg)
      node_reg_mem = node_reg->get_memory_usage();
    MemorySizeType node_reg_mem_sum = node_reg_mem;
    MemorySizeType node_reg_mem_max = node_reg_mem;
    stk::all_reduce( comm, stk::ReduceSum<1>( &node_reg_mem_sum ) );
    stk::all_reduce( comm, stk::ReduceMax<1>( &node_reg_mem_max ) );

    // TODO: 2014-09-11 compute_memory_usage is deprecated in recent versions of stk, replace
    //stk::mesh::MemoryUsage mem_usage;
    //stk::mesh::compute_memory_usage(bulkData, mem_usage);
    //MemorySizeType mem_total_bytes_sum = mem_usage.total_bytes;
    //MemorySizeType mem_total_bytes_max = mem_usage.total_bytes;
    //stk::all_reduce( comm, stk::ReduceSum<1>( &mem_total_bytes_sum ) );
    //stk::all_reduce( comm, stk::ReduceMax<1>( &mem_total_bytes_max ) );
    if (bulkData.parallel_rank() == 0)
      {
        if (dump_level > 1)

          {
            /*
              std::cout << "P[" << bulkData.parallel_rank() << "] AdaptMain::memory_dump stk_mesh counted memory usage at stage [" << msg << "] "
              " parallel sum, max memory [MB]= " << ((double)mem_total_bytes_sum)/MB << " , " << ((double)mem_total_bytes_max)/MB << std::endl;
              if (dump_level > 2)
              stk::mesh::print_memory_usage(mem_usage, std::cout);
            */

            std::cout << "P[" << bulkData.parallel_rank() << "] AdaptMain::memory_dump rss total (sum all proc) at stage [" << msg << "] = "
                      << " current high-water-mark [MB]= " << ((double)rss_current_0)/MB
                      << " , " << ((double)rss_high_water_mark_0)/MB
                      << std::endl;
          }

        {
          std::cout << "AdaptMain::memory_dump summary for " << msg << ": "
            //<< "stk_mesh [sum]"     << ((double)mem_total_bytes_sum)/MB << "MB , "
                    << "NodeRegistry [sum]" << ((double)node_reg_mem_sum)/MB << "MB , "
                    << "rss[sum] [MB]"      << ((double)rss_current_0)/MB << "MB"
                    << std::endl;
        }

      }
    //if (rss_current_0 )
    returned_total_memory = rss_current_0;
    //else
    //  returned_total_memory = mem_total_bytes_sum+node_reg_mem_sum;

    return returned_total_memory;
  }

  //extern void test_memory(int, int);
  void MeshAdapt::test_memory(MemorySizeType n_elements, MemorySizeType n_nodes)
  {
    vector<stk::mesh::Entity> new_elements;
    vector<stk::mesh::Entity> new_nodes;

    eMeshP->reopen();
    //stk::mesh::Part& part = eMeshP->get_fem_meta_data()->declare_part("elems", stk::topology::ELEMENT_RANK);
    stk::mesh::Part& part = eMeshP->get_fem_meta_data()->declare_part_with_topology("elems", stk::topology::HEX_8);

    eMeshP->commit();

    eMeshP->get_bulk_data()->modification_begin();

    std::cout << "creating " << n_elements << " elements..." <<std::endl;
    eMeshP->createEntities( stk::topology::ELEMENT_RANK, n_elements, new_elements);
    std::cout << "... done creating " << n_elements << " elements" << std::endl;

    std::cout << "creating " << n_nodes << " nodes..." <<std::endl;
    eMeshP->createEntities( stk::topology::NODE_RANK, n_nodes, new_nodes);
    std::cout << "... done creating " << n_nodes << " nodes" << std::endl;

    MemorySizeType num_prints = std::min(static_cast<MemorySizeType>(100ul), n_elements);
    MemorySizeType print_mod = n_elements/num_prints;
    MemorySizeType i_node = 0;
    int n_node_per_element = 8;
    for (MemorySizeType i_element = 0; i_element < n_elements; i_element++)
      {
        if (!i_element || (i_element % print_mod == 0))
          {
            std::cout << "declare_relation for i_element = " << i_element << " [" << n_elements << "] = " << ((double)i_element)/((double)n_elements)*100 << "%"
                      << std::endl;
          }

        stk::mesh::Entity element = new_elements[i_element];
        std::vector<stk::mesh::Part*> add_parts(1,&part);
        std::vector<stk::mesh::Part*> remove_parts;
        eMeshP->get_bulk_data()->change_entity_parts( element, add_parts, remove_parts );

        for (int j_node = 0; j_node < n_node_per_element; j_node++)
          {
            stk::mesh::Entity node = new_nodes[i_node];

            eMeshP->get_bulk_data()->declare_relation(element, node, j_node);

            i_node++;
            if (i_node >= n_nodes-1)
              i_node = 0;
          }
      }

    std::cout << " doing modification_end ... " << std::endl;
    stk::mesh::fixup_ghosted_to_shared_nodes(*eMeshP->get_bulk_data());
    eMeshP->get_bulk_data()->modification_end();
    std::cout << " done modification_end ... " << std::endl;
  }

  // Determine the basename and refinement level from a given an ExodusII file name.
  // This is the last integer before the final '.e' or '.exo' in the file name.
  // If no level was found, returns 0 as the level.
  //
  // Return:
  //    string - the basename which is everything before the number
  //    int    - the refinement level
  std::pair<std::string,int> MeshAdapt::get_basename_and_level(const std::string& input_mesh)
  {
    // find the last number before the first dot character
    size_t dot_index = input_mesh.find_last_of(".");
    std::string prefix(input_mesh, 0, dot_index);
    size_t last_index = prefix.find_last_not_of("0123456789");

    int level = 0; // level will be one if we don't find a level
    if (last_index + 1 < dot_index) {
      std::string level_str(prefix, last_index + 1, std::string::npos);
      level = std::stoi(level_str);
    }
    std::string basename(prefix, 0, last_index + 1);

    return std::make_pair(basename, level);
  }

  void MeshAdapt::checkInput(std::string option, std::string value, std::string allowed_values, Teuchos::CommandLineProcessor& clp)
  {
    //if (value.length() == 0) return;
    std::vector<std::string> vals = Util::split(allowed_values, ", ");
    for (unsigned i = 0; i < vals.size(); i++)
      {
        if (vals[i] == value)
          return;
      }
    std::ostringstream oss;
    oss << "\nCommand line syntax error: bad option for " << option << " (= " << value << ") \n allowed values = " << allowed_values;
    //throw std::runtime_error(oss.str());
    {
      std::cout << oss.str() << std::endl;
      printHelp(clp);
      std::cout << oss.str() << std::endl;
      exit(1);
    }

  }

  void MeshAdapt::print_simple_usage(int argc, char **argv)
  {
    std::string exe_name(argv[0]);
    std::cout << "usage: " << argv[0] <<
      " convert input_mesh_filename [output_mesh_filename]" << std::endl;
    std::cout << "       " << argv[0] <<
      " [refine|enrich] input_mesh_filename  [output_mesh_filename] [number_levels]" << std::endl;
    std::cout << "       " << argv[0] <<
      " adapt input_mesh_filename [analysis_results_filename]" << std::endl;
  }

  int MeshAdapt::check_for_simple_options(int argc, char **argv)
  {
    first_extra_option_index = -1;
    int simple = 0;
    for (int i = 1; i < argc; i++)
      {
        std::string sargv = std::string(argv[i]);
        if (sargv == "adapt" ||
            sargv  == "refine"  ||
            sargv == "enrich"  ||
            sargv == "convert")
          {
            for (int j = i+1; j < argc; j++)
              {
                if (Util::startsWith(sargv, "--"))
                  {
                    first_extra_option_index = j;
                    break;
                  }
              }
            return i;
          }
      }
    return simple;
  }

  int MeshAdapt::adapt_main_simple_options(int argc_in, char **argv_in)
  {
    first_extra_option_index = -1;
    int simple_options_index = check_for_simple_options(argc_in, argv_in);

    int argc = argc_in;
    if (argc < 2 + simple_options_index ||
        argc > 4 + simple_options_index ) {
      if (first_extra_option_index == -1 || first_extra_option_index == simple_options_index + 1)
        print_simple_usage(argc_in, argv_in);
      return 1;
    }

    std::vector<std::string> argv(argc);
    for (int i = 0; i < argc; i++)
      argv[i] = (const char *)argv_in[i];

    std::string exe_name        = argv[0];
    std::string option          = argv[0+simple_options_index];
    if (option != "adapt" && option != "refine" && option != "enrich" && option != "convert") {
      print_simple_usage(argc_in, argv_in);
      return 1;
    }
    std::string input_mesh = argv[1+simple_options_index];
    std::string number_refines = (argc == 3+simple_options_index? argv[2+simple_options_index] : "1");
    int nref=0;

    bool isInt = false;
    try {
      nref = boost::lexical_cast<int>(number_refines);
      (void)nref;
      isInt = true;
    }
    catch( ... )
      {
      }

    std::string output_mesh = input_mesh;
    std::string extension = input_mesh.substr(input_mesh.length()-2,input_mesh.length());
    if (debug) std::cout << " extension= " << extension << std::endl;
    std::string new_ext = "";
    new_ext += "_";
    if (option == "refine")
      new_ext += "refined_"+number_refines+extension;
    else
      new_ext += option+"ed_"+extension;

    std::string previous_adapted_mesh("");
    std::string next_adapted_mesh("");

    //
    // determine name of output mesh, and other auxiliary meshes
    //
    if (option == "adapt") {
      std::pair<std::string,int> basename_and_level = get_basename_and_level(input_mesh);
      std::string basename(basename_and_level.first);
      int input_level = basename_and_level.second;

      size_t dot_index = input_mesh.find_last_of(".");
      std::string suffix = input_mesh.substr(dot_index+1);

      if (input_level > 0) {
        previous_adapted_mesh = basename + std::to_string(input_level) + "_ft." + suffix;
      }
      std::string next_level(std::to_string(input_level + 1));
      next_adapted_mesh = basename + next_level + "_ft." + suffix;
      output_mesh = basename + next_level + "." + suffix;

      // set to name of analysis mesh, if it exists
      if (argc >= 2 + simple_options_index)
        input_mesh = argv[2+simple_options_index];

      number_refines = "1";
    }
    else {
      if (!isInt && (argc == 3+simple_options_index)) {
        output_mesh = number_refines;
        number_refines = (argc == 4+simple_options_index? argv[3+simple_options_index] : "1");
        if (Util::startsWith(number_refines, "--"))
          number_refines = "";
        std::cout << "tmp output_mesh= " << output_mesh << std::endl;
        std::cout << "tmp number_refines= " << number_refines << std::endl;
      }
      else {
        Util::replace(output_mesh, extension, new_ext);
      }
    }
    if (debug) std::cout << " output_mesh= " << output_mesh << std::endl;

    std::vector<std::string> argv_new;
    for (int i = 0; i < simple_options_index; i++)
      argv_new.push_back(argv[i]);
    argv_new.push_back("--input_mesh="+input_mesh);
    argv_new.push_back("--output_mesh="+output_mesh);
    if (option == "refine")
      argv_new.push_back("--refine=DEFAULT");
    else if (option == "adapt") {
      argv_new.push_back("--refine=DEFAULT");
      argv_new.push_back("--respect_spacing=0");
      if (previous_adapted_mesh != "")
        argv_new.push_back("--previous_adapted_mesh=" + previous_adapted_mesh);
      argv_new.push_back("--next_adapted_mesh=" + next_adapted_mesh);
      argv_new.push_back("--RAR_info=\"adapt.yaml\"");
    }
    else if (option == "enrich")
      argv_new.push_back("--enrich=DEFAULT");
    else
      argv_new.push_back("--convert=Hex8_Tet4_24");
    argv_new.push_back("--load_balance=0");
    argv_new.push_back("--remove_original_elements=1");
    argv_new.push_back("--number_refines="+number_refines);
    argv_new.push_back("--ioss_read_options=\"large\"");
    argv_new.push_back("--ioss_write_options=\"large\"");

    if (first_extra_option_index > 0)
      {
        for (int i = first_extra_option_index; i < argc; ++i)
          {
            argv_new.push_back(argv[i]);
          }
      }
    if ( debug) std::cout << "new argv = \n" << argv_new << std::endl;
    int argc_new = argv_new.size();
    char **argv_new_cstr = new char*[argc_new];
    for (int i = 0; i < argc_new; i++)
      {
        argv_new_cstr[i] = (char *)argv_new[i].c_str();
      }
    int ret_val = adapt_main_full_options(argc_new, argv_new_cstr);
    delete[] argv_new_cstr;
    return ret_val;
  }

  void MeshAdapt::dump_args(int argc, char **argv)
  {
    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
      {
        std::cout << "argv[" << i << "]= " << argv[i] << std::endl;
      }
  }

  void MeshAdapt::check_args(int argc, char **argv)
  {
    std::string errors="";
    for (int i = 1; i < argc; i++)
      {
        std::string av(argv[i]);
        if (av == "--help" || av == "--Help" || av == "-h" || av == "--h" ||
            av == "-H" || av == "--H" || av == "--version" || av == "-v") continue;
        size_t equal_pos = av.find("=",0);
        if (equal_pos == std::string::npos)
          {
            errors += "ERROR in options: no '=' found in option: "+av+"\n";
          }
        if (av.length() < equal_pos+2)
          {
            errors += "ERROR in options: no value given after '=', found in option: "+av+"\n";
          }
      }
    if (errors.length())
      {
        std::cout << "ERRORS found in options: debug dump of options= " << std::endl;
        dump_args(argc, argv);
        std::cout << errors << std::endl;
        throw std::runtime_error(errors);
      }
  }

  int MeshAdapt::adapt_main(int argc, char **argv)
  {
    // allow positional arguments, etc.
    if (check_for_simple_options(argc, argv))
      return adapt_main_simple_options(argc, argv);
    else
      {
        check_args(argc, argv);
        return adapt_main_full_options(argc, argv);
      }
  }

  int MeshAdapt::setup_options(Teuchos::CommandLineProcessor& clp, const stk::ParallelMachine& comm, int argc, char **argv)
  {
    bool found_Help = false;
    bool found_print_version = false;
    for (int i = 0; i < argc; ++i) {
      const std::string s(argv[i]);
      if ( s == "-H" || s == "-Help" || s == "--Help" || s == "--H") {
        found_Help = true;
      }
      if (s == "--version" || s == "-v")
        found_print_version = true;
    }

    clp.setOption("help"                     , &help                     , "print this usage message");

    // version
    clp.setOption("version"                  , &print_version            , "print version and exit");

    // files
    clp.setOption("input_mesh"               , &input_mesh               , "input mesh name");
    clp.setOption("output_mesh"              , &output_mesh              , "output mesh name");
    clp.setOption("load_balance"             , &load_balance             , "load balance (decomp/slice/spread) input mesh file");
    clp.setOption("previous_adapted_mesh"    , &previous_adapted_mesh    , "previous adapted mesh name for offline adaptivity (INTERNAL USE)");
    clp.setOption("next_adapted_mesh"        , &next_adapted_mesh        , "next adapted mesh name for offline adaptivity (INTERNAL USE)");
    clp.setOption("save_intermediate_meshes" , &save_intermediate_meshes ,
                                  "save meshes in refinement sequence to 'refined_mesh_I.e' for each I <= number_refines\n"
                                  "  Note: if < 0, save in animation format: refined_mesh.e, refined_mesh.e-s0001, ...");

    // operation type
    clp.setOption("refine"                   , &refine                   , refine_options.c_str());
    clp.setOption("number_refines"           , &number_refines           ,
                                  "number of refinement passes"
                                  "\n  Note: must be set to be >= max of any 'N' in the --blocks=<block-name>:Nx,... option."
                                  "\n  Note: if the :Nx optional specification is not used for a block, it is equivalent"
                                  "\n    to specifying :<number_refines>x for that block.");
    clp.setOption("convert"                  , &convert                  , convert_options.c_str());
    clp.setOption("enrich"                   , &enrich                   , enrich_options.c_str());

    // run-adapt-run
    clp.setOption("RAR_info"                 , &adapt_input              , "name of input file for advanced usage");

    // spacing
#if !defined(NO_GEOM_SUPPORT)
    clp.setOption("respect_spacing"          , &respect_spacing          , "respect the initial mesh spacing during refinement");
#endif

    // query/control
    clp.setOption("verify_meshes"            , &verify_meshes            , "verify positive volumes for input and output meshes, set to 1 for finite element volume checks, FV for finite volume checks");
    clp.setOption("query_only"               , &query_only               , "query only, no refinement done");
    clp.setOption("progress_meter"           , &progress_meter           , "progress meter on or off");
    clp.setOption("print_timers"             , &print_timers             , "print more detailed timing info");
    clp.setOption("print_info"               , &print_info               , ">= 0  (higher values print more info)");
    clp.setOption("skin_mesh"                , &skin_mesh                , "produce a boundary sideset for any exposed boundaries");
    clp.setOption("dihedral_angle_check"     , &dihedral_angle_check     ,
                                  ">= 0  check dihedral angles of a tetrahedral mesh, print those that are larger than 90 deg. (higher values print more info)\n");
    clp.setOption("print_memory_usage"       , &print_memory_usage       , "print memory usage");
    clp.setOption("precheck_memory_usage"    , &precheck_memory_usage    ,
                                  "before refinement, check there's enough memory by using an estimated number of bytes per element (1000)\n"
                                  "\nIf set to 1, we internally estimate how many cores/node using a list of SNL hostnames, but it can be set\n"
                                  "explicitly to the number of cores/node, e.g. on chama, cores/node is 16");
    clp.setOption("estimate_memory_usage"    , &estimate_memory_usage    ,
                                  " use internal or memory_multipliers_file values to estimate memory needed.\n"
                                  "   if query_only=1, use multipliers from memory_multipliers_file to estimate memory to be used in refinements, if memory_multipliers_file is set.\n"
                                  "   if query_only=1, and no memory_multipliers_file is set, use internal values for memory_multipliers.\n"
                                  "   If query_only=0, print actual memory data and estimates.");

    // properties
    clp.setOption("property_map"             , &property_map             , "YAML-based list of string: string pairs, e.g. {smoother_type: Newton, smoother_niter: 200, smoother_tol: 1.e-3}");

    // geometry
    clp.setOption("input_geometry"           , &input_geometry           , "input geometry name");
    clp.setOption("smooth_geometry"          , &smooth_geometry          , "smooth geometry - moves nodes after geometry projection to try to avoid bad meshes");
    clp.setOption("smoother_niter"           , &smoother_niter           , "mesh smoother number of iterations");
    clp.setOption("smoother_tol"             , &smoother_tol             , "mesh smoother convergence tolerance");
    clp.setOption("smooth_surfaces"          , &smooth_surfaces          , "allow nodes to move on surfaces when smoothing");
    clp.setOption("dump_geometry_file"       , &dump_geometry_file       , "debug print geometry (OpenNURBS 3dm) file contents");
    clp.setOption("fit_geometry_file"        , &fit_geometry_file        , "for 2D meshes, create an OpenNURBS 3dm file fitting cubic splines to all boundaries");
    clp.setOption("fit_angle"                , &fit_angle                ,
                                  "for 2D meshes, specify angle criterion for determining corners - defined by the\n"
                                  "   angle between the two adjacent faces, so a flat surface has angle 180.\n"
                                  "   Specified in degrees - if the included angle is less than fit_angle, the\n"
                                  "   node is considered a corner and a hard corner will be enforced.");
    clp.setOption("fit_3d_file"              , &fit_3d_file              ,
                                  "for 3D meshes, fit bi-cubics to mesh surface geometry, store in fields for evaluation\n"
                                  "  specify a YAML file to read with control information on which surfaces\n"
                                  "  to fit and how to treat seams in the surface using an angle criterion\n"
                                  "  a sample YAML file will be printed if fit_3d_file is set to sample.yaml.\n"
                                  "Note: if the YAML file has 'QA: ... data ...' set, then the surfaces will be\n"
                                  "  converted to shells, and all seams will be detected and put in a\n"
                                  "  special edge part and saved to the file specified - this is\n"
                                  "  for QA'ing of the input data before doing the actual fit.");

    clp.setOption("convert_geometry_parts_OpenNURBS"  , &convert_geometry_parts_OpenNURBS , "remove beam and shell elements, place nodes in similarly-named nodesets, save to specified new file in this option");

    // smoothing
    clp.setOption("smooth_use_reference_mesh", &smooth_use_reference_mesh, "for most cases, set to 1 (default) - can be used for smoothing with no reference mesh");
    clp.setOption("fix_all_block_boundaries" , &fix_all_block_boundaries , "when smoothing without geometry, fix all inner and outer block boundaries");

    // mesh query
    clp.setOption("compute_hmesh"            , &compute_hmesh            , "compute mesh parameter using method eigens|edges");
    clp.setOption("print_hmesh_surface_normal"  , &print_hmesh_surface_normal            , "prints a table of normal mesh spacing at each surface");

    // subsetting
    clp.setOption("blocks"                   , &block_name_inc           , block_name_desc_inc.c_str());
    clp.setOption("use_transition_elements"  , &use_transition_elements  ,
                                  "when specifying --blocks option, set this to 1 to enable conformal meshes\n"
                                  "by introducing so-called transition elements in the border region between\n"
                                  "a refined and unrefined block");
    // histograms
#if STK_ADAPT_HAVE_YAML_CPP
    clp.setOption("histogram"  , &histogram_options  , "options for histograms:"
                                  "\n  either a single filename, which reads further commands in YAML format (yaml.org) from that file, or\n"
                                  "  a string of the form \"{ fields: [field_1,...,field_n], file_root: my_histograms,\n"
                                  "     mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume], time: 0.1, step: 2 }\"\n"
                                  "  where field_i are field names to get stats for, file_root is a root filename,\n"
                                  "  and mesh: gives options for mesh quality histograms.\n"
                                  "  time: or step: options allow specifying which timestep in the database should be used.\n"
                                  "  If read from a file, file format is like this:\n"
                                  "    fields:\n"
                                  "      - pressure\n"
                                  "      - velocity\n"
                                  "      - temperature\n"
                                  "    file_root: my_histograms\n"
                                  "    mesh:\n"
                                  "      - edge_length\n"
                                  "      - quality_edge\n"
                                  "      - quality_vol_edge_ratio\n"
                                  "      - volume");
#endif
    clp.setOption("histogram_file_root"      , &histograms_root          , " if cout, use screen, else use this as the root name of histogram files.");

    // ioss options
    clp.setOption("ioss_read_options"        , &ioss_read_options        ,
                                  "options to IOSS/Exodus for e.g. large files | auto-decomp | auto-join\n"
                                  "to use, set the string to a combination of\n"
                                  "{\"large\", \"auto-decomp:yes\",  \"auto-decomp:no\",\n"
                                  "   \"auto-join:yes\", \"auto-join:no\" }\n"
                                  "   e.g. \"large,auto-decomp:yes\"\n"
                                  " Note: set options for read and/or write (ioss_write_options)");
    clp.setOption("ioss_write_options"       , &ioss_write_options       , "see ioss_read_options");
    clp.setOption("generated_mesh"           , &generated_mesh           ,
                  "use internal fixture to generate an NxMxP hexahedral mesh and write it out to --output_mesh, --input_mesh should contain the following\n"
                  "  syntax: NxMxP|bbox:xmin,xmax,ymin,ymax,zmin,zmax|sideset:xXyYzZ (bbox and sideset are optional; sideset args are where to put sidesets\n"
                  "     e.g: '--generated_mesh=1 --input_mesh=\"10x10x10|bbox:-1,-1,-1,2,2,1|sideset:xXyZ\"' would generate a 10x10x10 mesh with sidesets on minX,maxX,minY and maxZ sides"
                  );

    // debugging/advanced
    clp.setOption("debug"                    , &debug                    , " turn on debug printing");

    clp.setOption("use_side_map"             , &use_side_map             , " experimental - for shell-based meshes");
    clp.setOption("proc_rank_field"          , &proc_rank_field          , " add an element field to show processor rank");
    clp.setOption("remove_original_elements" , &remove_original_elements , " remove original elements (default=true)");

    clp.setOption("remove_geometry_blocks"   , &remove_geometry_blocks   , "remove geometry blocks from output Exodus file after refinement/geometry projection");

    // internal
    clp.setOption("delete_parents"           , &delete_parents           , "DEBUG: delete parents from a nested, multi-refine mesh - used for debugging");
    clp.setOption("snap_geometry"            , &snap_geometry            , "project nodes to geometry - used for internal testing only");
    clp.setOption("internal_test"            , &internal_test            , "run the specified internal test");
    clp.setOption("sync_io_regions"          , &sync_io_regions          ,
                                  "synchronize input/output region's Exodus id's (default=0)\n"
                                  "   use this option if you want to ensure output mesh has\n"
                                  "   the same block ids and names as the input mesh, which\n"
                                  "   only makes sense if in refine (not enrich or convert) mode");

#if defined(STK_BUILT_IN_SIERRA)
    // Salinas
    clp.setOption("rbar_blocks"              , &rbar_blocks              , "list of blocks to treat in special Salinas fashion for RBARs - see block_name description for format.");
#endif
    //clp.setOption("exclude"                  , &block_name_exc           , block_name_desc_exc.c_str());

    clp.setOption("memory_multipliers_file"  , &memory_multipliers_file  ,
                                  "[experimental]  filename with 3 space-separated entries, with estimate for bytes-per-hex8 tet4 and nodes, e.g. 300 280 200\n"
                                  "  If not set, use internal estimates for memory multipliers.");

#if ALLOW_MEM_TEST
    clp.setOption("test_memory_elements"     , &test_memory_elements     , " give a number of elements");
    clp.setOption("test_memory_nodes"        , &test_memory_nodes        , " give a number of nodes");
#endif
    clp.setOption("serialized_io_group_size" , &serialized_io_group_size , "[experimental] set to non-zero to use this many i/o groups to minimize disk contention");

    // Detailed doc
    std::string docString = s_docString;

#ifdef STK_BUILT_IN_SIERRA
    std::string docStringSalinas =
      "Salinas Special Command --rbar_blocks\n"
      "\n"
      "Percept can refine a mesh that contains Salinas RBAR elements (beams connecting nodes of one surface\n"
      "to another, e.g. to model a joint or spring).  The new mesh will contain RBARs connecting new nodes\n"
      "on one surface to new nodes on the other.  Specify a list of block names that contain the RBAR\n"
      "elements (see the --blocks command for format).\n";

    docString = docString + docStringSalinas;
#endif

    if (found_print_version)
      {
        std::string v;
        int vb = get_version(&v);
        if(p_rank == 0) std::cout << "Version: " << v << std::endl;
#if defined( STK_HAS_MPI )
        stk::parallel_machine_finalize();
#endif
#if defined(WITH_KOKKOS)
        Kokkos::finalize();
#endif
        std::exit(vb);
      }

    if (!found_Help) docString = "";
    clp.setDocString(docString.c_str());

    if (found_Help) {
      try {
        printHelp(clp);
      }
      catch (...) {}
#if defined( STK_HAS_MPI )
      stk::parallel_machine_finalize();
#endif
#if defined(WITH_KOKKOS)
      Kokkos::finalize();
#endif
      std::exit(0);
    }

    int err_clp = processCommandLine(clp, argc, argv);
    if (err_clp) return err_clp;

    if (serialized_io_group_size)
      {
        std::cout << "Info: found non-zero serialized_io_group_size on command-line= "
                  << serialized_io_group_size << std::endl;
        if (serialized_io_group_size < 0 || serialized_io_group_size > (int)p_size || (int)p_size % serialized_io_group_size != 0)
          {
            if (p_rank==0)
              std::cout << "Error: Job requested serialized_io_group_size of " << serialized_io_group_size
                        << " which is incompatible with MPI size= " << p_size
                        << "... shutting down." << std::endl;
            throw std::runtime_error("bad value for serialized_io_group_size");
          }
        Ioss::SerializeIO::setGroupFactor(serialized_io_group_size);
      }

    if (convert.length()+enrich.length()+refine.length() == 0)
      {
        std::cout << "\nCommand line syntax error: you must give a value for one (and only one) of the options: refine, enrich, or convert.\n" << std::endl;
        printHelp(clp);
        std::cout << "\nCommand line syntax error: you must give a value for one (and only one) of the options: refine, enrich, or convert.\n" << std::endl;
        exit(1);
      }

    if (convert.length())
      checkInput("convert", convert, convert_options, clp);

    if (enrich.length())
      checkInput("enrich", enrich, enrich_options, clp);

    if (refine.length())
      checkInput("refine", refine, refine_options, clp);

    if (print_info || convert_geometry_parts_OpenNURBS.size())
      {
        doRefineMesh = false;
      }

    if (help
        || input_mesh.length() == 0
        || output_mesh.length() == 0
        || ((convert.length() == 0 && refine.length()==0 && enrich.length()==0) && number_refines)
        //||  not (convert == "Hex8_Tet4_24" || convert == "Quad4_Quad4_4" || convert == "Quad4_Tri3_6")
        )
      {
        printHelp(clp);
        return 1;
      }

#if defined( STK_HAS_MPI )
    MPI_Barrier( MPI_COMM_WORLD );
#endif

    // check geometry file extensions
    if (input_geometry.length())
      {
        if (p_rank == 0) std::cout << "mesh_adapt: input_geometry type  file= " << input_geometry << std::endl;
        if (input_geometry.find(".3dm") != std::string::npos )
          {
            std::string m2gFile = input_geometry.substr(0,input_geometry.length()-3) + "m2g";

            struct stat s;
            if (0 == stat(m2gFile.c_str(), &s))
            {
#ifdef HAVE_CUBIT
              input_geometry_type = PGEOM_OPENNURBS;
              if (p_rank == 0) std::cout << "mesh_adapt: input_geometry type is OpenNURBS (PGeom) for file= " << input_geometry << std::endl;
#else
              throw std::runtime_error("CUBIT not supported on this platform");
#endif              
            }
            else
            {
              input_geometry_type = OPENNURBS;
              if (p_rank == 0) std::cout << "mesh_adapt: input_geometry type is OpenNURBS for file= " << input_geometry << std::endl;
            }
          }
        else if (input_geometry.find(".e") != std::string::npos || input_geometry.find(".g") != std::string::npos || input_geometry.find(".exo") != std::string::npos)
          {
            input_geometry_type = MESH_BASED_GREGORY_PATCH;
            if (p_rank == 0) std::cout << "mesh_adapt: input_geometry type is mesh-based (Gregory patch) for file= " << input_geometry << std::endl;
          }
        else if(input_geometry.find(".sat") != std::string::npos)
        {
#ifdef HAVE_ACIS
          input_geometry_type = PGEOM_ACIS;
          if (p_rank == 0) std::cout << "mesh_adapt: input_geometry type is ACIS for file= " << input_geometry << std::endl;
#else
          throw std::runtime_error("ACIS is not available in this platform.");
#endif
        }
        else
          {
            VERIFY_MSG("invalid file extension on --input_geometry file\n   "
                       "-- valid extensions are .3dm (OpenNURBS) or .e,.g,.exo for GregoryPatch Exodus files - file= " + input_geometry);
          }
      }


    return 0;
  }

  void MeshAdapt::do_precheck_memory_usage()
  {
    size_t avail = 0;
    std::string hostname;
    int nnodes = 0;

    if (precheck_memory_usage)
      {
        char host[256];
        gethostname(host, sizeof(host));
        hostname = host;
        std::string names[] = {"chama","redsky","glory","skybridge","muzia","curie","uno"};
        int cores_per_node[] = {16,     8,       16,     16,         8,      8,      8};
        int nnames = sizeof(names)/sizeof(std::string);
        int cpn = 0;
        if (precheck_memory_usage == 1)
          {
            for (int in = 0; in < nnames; ++in)
              {
                if (hostname.compare(0, names[in].length(), names[in]) == 0)
                  {
                    cpn = cores_per_node[in];
                    if (!eMeshP->get_rank())
                      {
                        std::cout << "NOTE: found hostname in our internal list, hostname= " << hostname << std::endl;
                      }
                    break;
                  }
              }
          }
        else
          {
            cpn = precheck_memory_usage;
          }
        bool unknown_host = false;
        if (!cpn)
          {
            unknown_host = true;
            cpn = 8;
            if (!eMeshP->get_rank())
              std::cout << "WARNING: unknown host = " << hostname << " using 8 for cores_per_node - note: if this is incorrect, re-run with --precheck_memory_usage=<your cores_per_node> or 0 to turn off memory checking."
                        << std::endl;
          }

        uint64_t nelLocal = static_cast<uint64_t> (eMeshP->get_ioss_mesh_data()->get_input_io_region()->get_property("element_count").get_int());
        uint64_t nelGlobal = static_cast<uint64_t> (eMeshP->get_ioss_mesh_data()->get_input_io_region()->get_property("global_element_count").get_int());
        uint64_t nelMin, nelMax, nelAve;
        stk::get_max_min_avg(eMeshP->parallel(), nelLocal, nelMax, nelMin, nelAve);
        uint64_t bytes_per_elem = 1300; // our estimate
        uint64_t mem = nelGlobal*bytes_per_elem;

        uint64_t fac = std::pow(uint64_t(8), uint64_t(std::max(number_refines,1)));
        if (number_refines == 0)
          fac = 0;
        mem *= (fac + 1 + 1);  // uniform refinement estimate + memory for initial mesh + 1 more for luck

        stk::get_memory_available(avail);
        if (avail == 0) avail = 1;
        nnodes = eMeshP->get_parallel_size() / cpn;
        nnodes = std::max(nnodes, 1);
        uint64_t totAvailableMem = avail*nnodes;
        double Gb = 1024.0*1024.0*1024.0;
        if (eMeshP->get_rank() == 0)
          {
            std::cout << "P[" << eMeshP->get_rank() << "] precheck_memory_usage is set, before reading the bulk data of the mesh, the computed available memory= "
                      << double(totAvailableMem)/Gb << " Gb" << " while the estimated needed memory to refine the mesh is " << double(mem)/Gb << " Gb"
                      << "\nThe memory available on each node is = " << double(avail)/Gb
                      << "\nThe local element count = " << nelLocal << " global element count= " << nelGlobal
                      << " min = " << nelMin << " max= " << nelMax << " ave= " << nelAve
                      << "\nThe number of compute-nodes = " << nnodes << " using cores_per_node= " << cpn
                      << "\nThe number of cores = " << eMeshP->get_parallel_size()
                      << std::endl;
          }
        if (mem > totAvailableMem)
          {
            if (unknown_host)
              {
                if (eMeshP->get_rank() == 0)
                  std::cout << "WARNING: not enough memory to read and refine the mesh, use more procs/cores, we estimate needed nodes (minimum) = " + toString(mem/avail) +
                    "\nHowever, since this host is unknown, execution will continue" << std::endl;
              }
            else
              {
                throw std::runtime_error("not enough memory to read and refine the mesh, use more procs/cores, we estimate needed nodes (minimum) = " + toString(mem/avail));
              }
          }
      }
  }

  void MeshAdapt::mesh_based_geometry_setup()
  {
    FitGregoryPatches::SurfaceSets surfaceSets;
    FitGregoryPatches::AngleMap angleMap;
    fitter.reset( new FitGregoryPatches(*eMeshP, surfaceSets, angleMap, 15.0));
    if (fit_3d_file.length())
      {
        fitter->parse(fit_3d_file);
        bool doRegister=true;
        fitter->register_or_set_fields(doRegister);
      }
    else if (input_geometry_type == MESH_BASED_GREGORY_PATCH)
      {
        bool doRegister=true;
        fitter->register_or_set_fields(doRegister);
        eMeshP->add_input_field(eMeshP->m_gregory_control_points_field);
        eMeshP->add_input_field(eMeshP->m_gregory_control_points_field_shell);
        eMeshP->add_input_field(eMeshP->m_node_normals);
      }

    grsf.reset(new GeometryRecoverySplineFit(*eMeshP, fit_angle));
#if defined(STK_PERCEPT_HAS_GEOMETRY)
    if (fit_geometry_file != "")
      {
        grsf->fit_geometry_create_parts_meta();
      }
#elif defined(NO_GEOM_SUPPORT)
    std::ostringstream oss;
    oss << "\nERROR: Geometry and/or smoothing is not currently supported on this platform. Try running with geometry turned off.";
    throw std::runtime_error(oss.str());
#endif
  }
  void MeshAdapt::mesh_based_geometry_fitting()
  {
    if (fit_3d_file.length())
      {
        if (!eMeshP->get_rank()) std::cout << "Fitting 3D geometry using control data from file: " << fit_3d_file
                                           << " QA.activate= " << fitter->m_QA.m_activate << " QA.file= " << fitter->m_QA.m_file
                                           << std::endl;
        if (fitter->m_QA.m_activate && fitter->m_QA.m_file.length())
          {
            fitter->fit();
            fitter->createQA(fitter->m_QA.m_file+"_"+input_mesh);
            if (fitter->m_QA.m_visualizer_command_prefix.length() && p_rank == 0)
              {
                std::ostringstream sz, rank;
                if (p_size > 1)
                  {
                    int width = int(std::log10(double(p_size))) + 1;
                    sz << std::setfill('0') << std::setw(width) << p_size;
                    rank << std::setfill('0') << std::setw(width) << p_rank;
                  }
                std::string ext =  (p_size == 1 ? "" : "."+sz.str()+"."+rank.str());
                for (unsigned ii=0; ii < fitter->m_QA.m_filenames.size(); ++ii)
                  {
                    std::string command = fitter->m_QA.m_visualizer_command_prefix + fitter->m_QA.m_filenames[ii] + ext;
                    std::cout << "\n\n======\n=====\nrunning visualization command: " << command << "\n\n" << std::endl;
                    runCommand(command);
                    if (fitter->m_QA.m_num_divisions)
                      {
                        command = fitter->m_QA.m_visualizer_command_prefix + "div_" + fitter->m_QA.m_filenames[ii] + ext;
                        std::cout << "\n\n======\n=====\nrunning visualization command: " << command << "\n\n" << std::endl;
                        runCommand(command);
                      }
                  }
              }
          }
        else
          {
            fitter->fit();
          }
      }

#if defined(STK_PERCEPT_HAS_GEOMETRY)
    if (fit_geometry_file != "")
      {
        grsf->fit_geometry(fit_geometry_file);
      }
#endif
  }

  void MeshAdapt::verify_mesh_util(bool isInit)
  {
    if (verify_meshes == "1" || verify_meshes == "2" || verify_meshes == "FV")
      {
        bool print_table=true;
        double badJac=0.0;
        bool use_finite_volume = false;
        if (verify_meshes == "FV")
          use_finite_volume = true;
        int dump_all_elements = 0;
        if (verify_meshes == "2")
          dump_all_elements = 1;
        std::string type = (isInit ? "input":"refined");
        if (!eMeshP->get_rank()) std::cout << "Verify " << type << " mesh..." << std::endl;
        if (eMeshP->check_mesh_volumes(print_table, badJac, dump_all_elements, use_finite_volume))
          {
            throw std::runtime_error("ERROR: verify_meshes shows a bad "+type+" mesh");
          }
        if (isInit) adaptedMeshVerifier.reset( new percept::AdaptedMeshVerifier(true));
        if (!adaptedMeshVerifier->isValid(*eMeshP, isInit))
          throw std::runtime_error("ERROR: AdaptedMeshVerifier found invalid "+type+" mesh");
      }
    else
      {
        if (verify_meshes != "0")
          VERIFY_MSG("--verify_meshes option unrecognized, use 0, 1, 2 or FV, your optoin: "+verify_meshes);
      }
  }

  void MeshAdapt::do_histograms_init()
  {
#if STK_ADAPT_HAVE_YAML_CPP
    if (histogram_options.size() != 0)
      {
        Histograms<double> histograms(histograms_root);
        HistogramsParser<double> hparser(histogram_options);
        hparser.create(histograms);

        double hopt_time = histograms.m_database_time;
        int hopt_step = histograms.m_database_step;
        int current_step = eMeshP->get_current_database_step();
        if (hopt_time >= 0.0)
          {
            eMeshP->read_database_at_time(hopt_time);
          }
        else if (hopt_step >= 0)
          {
            eMeshP->read_database_at_step(hopt_step);
          }
        else
          {
            // get last step
            int step = eMeshP->get_database_time_step_count();
            eMeshP->read_database_at_step(step);
            //std::cout << "step= " << step << " current_step= " << current_step << std::endl;
            //eMeshP->read_database_at_step(step?step:1);
          }

        eMeshP->mesh_field_stats(&histograms);
        histograms.compute_uniform_bins(10);

        if (!p_rank) {
          std::cout << "Before refine, user-requested histograms= " << std::endl;
          histograms.print(true);
        }
        // reset
        eMeshP->read_database_at_step(current_step);
      }
#endif
  }

  void MeshAdapt::do_histograms_final()
  {
#if STK_ADAPT_HAVE_YAML_CPP
    if (histogram_options.size() != 0)
      {
        Histograms<double> histograms(histograms_root);
        HistogramsParser<double> hparser(histogram_options);
        hparser.create(histograms);

        eMeshP->mesh_field_stats(&histograms);
        histograms.compute_uniform_bins(10);

        if (!p_rank) {
          std::cout << "After refine, user-requested histograms= " << std::endl;
          histograms.print(true);
        }
      }
#endif
  }

  void MeshAdapt::compute_hmesh_sizes_init()
  {
    if (compute_hmesh.size() != 0)
      {
        double hmesh=0.0;
        Histograms<double> histograms(histograms_root);
        //Histograms<double> histograms;
#if STK_ADAPT_HAVE_YAML_CPP
        HistogramsParser<double> hparser(histogram_basic_options);
        hparser.create(histograms);
#endif
        if (compute_hmesh == "eigens")
          {
            hmesh = eMeshP->hmesh_stretch_eigens(hmesh_min_max_ave_factor, &histograms["stretch_eigens"].m_data, &histograms["quality_edge_eigens"].m_data);
            histograms["stretch_eigens"].set_titles("Stretch Eigens Histogram");
            histograms["quality_edge_eigens"].set_titles("Stretch Eigens Max/Min Quality Histogram");
          }
        else if (compute_hmesh == "edges")
          {
            hmesh = eMeshP->hmesh_edge_lengths(hmesh_min_max_ave_factor, &histograms["edge_length"].m_data, &histograms["quality_edge"].m_data);
            histograms["edge_length"].set_titles("Edge Length Histogram");
            histograms["quality_edge_eigens"].set_titles("Edge Max/Min Quality Histogram");
          }
        else
          {
            throw std::runtime_error("unknown option for compute_hmesh: "+compute_hmesh);
          }
        histograms.compute_uniform_bins(10);

        if (!p_rank) {
          std::cout << "Before refine, Mesh size (h-parameter) = " << hmesh
                    << " (min = " << hmesh_min_max_ave_factor[0]
                    << " max = " << hmesh_min_max_ave_factor[1]
                    << " ave = " << hmesh_min_max_ave_factor[2]
                    << ") "
                    << std::endl;
          histograms.print(true);
        }
        hmesh_factor = hmesh;
      }
  }

  void MeshAdapt::compute_hmesh_sizes_final()
  {
    if (compute_hmesh.size() != 0)
      {
        double hmesh=0.0;
        double min_max_ave[3];
        Histograms<double> histograms(histograms_root);
        //Histograms<double> histograms;
#if STK_ADAPT_HAVE_YAML_CPP
        HistogramsParser<double> hparser(histogram_basic_options);
        hparser.create(histograms);
#endif

        if (compute_hmesh == "eigens")
          {
            hmesh = eMeshP->hmesh_stretch_eigens(min_max_ave, &histograms["stretch_eigens"].m_data, &histograms["quality_edge_eigens"].m_data);
            histograms["stretch_eigens"].set_titles("Stretch Eigens Histogram");
            histograms["quality_edge_eigens"].set_titles("Stretch Eigens Max/Min Quality Histogram");
          }
        else if (compute_hmesh == "edges")
          {
            hmesh = eMeshP->hmesh_edge_lengths(min_max_ave, &histograms["edge_length"].m_data, &histograms["quality_edge"].m_data);
            histograms["edge_length"].set_titles("Edge Length Histogram");
            histograms["quality_edge_eigens"].set_titles("Edge Max/Min Quality Histogram");
          }
        else
          {
            throw std::runtime_error("unknown option for compute_hmesh: "+compute_hmesh);
          }
        histograms.compute_uniform_bins(10);

        hmesh_factor /= hmesh;
        hmesh_min_max_ave_factor[0] /= min_max_ave[0];
        hmesh_min_max_ave_factor[1] /= min_max_ave[1];
        hmesh_min_max_ave_factor[2] /= min_max_ave[2];
        if (!p_rank) {
          std::cout << "After refine, Mesh size (h-parameter) = " << hmesh << " oldH/newH factor= " << hmesh_factor
                    << "\n (new min = " << min_max_ave[0]
                    << " max = " << min_max_ave[1]
                    << " ave = " << min_max_ave[2]
                    << ") "
                    << "\n (old/new min = " << hmesh_min_max_ave_factor[0]
                    << " max = " << hmesh_min_max_ave_factor[1]
                    << " ave = " << hmesh_min_max_ave_factor[2]
                    << ") "
                    << std::endl;
          histograms.print(true);
        }
      }
  }

  void MeshAdapt::do_test_memory()
  {
    std::cout << "test_memory_elements and nodes are nonzero, will not refine exodus files." << std::endl;

    test_memory(test_memory_elements, test_memory_nodes);

    if (print_memory_usage)
      memory_dump(print_memory_usage, eMeshP->get_bulk_data()->parallel(), *eMeshP->get_bulk_data(), 0, "after test memory");

    if (estimate_memory_usage && !query_only)
      {
        MemorySizeType tot_mem = memory_dump(false, eMeshP->get_bulk_data()->parallel(), *eMeshP->get_bulk_data(), 0, "after test memory");

        //std::cout << "MemEst: num_nodes= " << test_memory_nodes << " num_tet4=0 hum_hex8= " << test_memory_elements << " memory= " << MegaByte(tot_mem) << std::endl;
        //MemoryMultipliers::process_estimate(tot_mem, eMesh, refiner.getRefinementInfoByType(), memory_multipliers_file);
        MemoryMultipliers memMults;
        if (memory_multipliers_file.size())
          memMults.read_simple(memory_multipliers_file);
        memMults.num_hex8=test_memory_elements;
        memMults.num_nodes=test_memory_nodes;
        MemorySizeType estMem = memMults.estimate_memory();
        //                 std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory= " << MegaByte(tot_mem)
        //                           << " estMem= " << MegaByte(estMem) << std::endl;
        std::cout << "MemEst: num_nodes= " << memMults.num_nodes << " num_tet4= " << memMults.num_tet4 << " num_hex8= " << memMults.num_hex8 << " memory[MB]= " << MegaByte(tot_mem)
                  << " estMem[MB]= " << MegaByte(estMem)
                  << " mult_hex8= " << memMults.mult_hex8 << " mult_tet4= " << memMults.mult_tet4 << " mult_nodes=" << memMults.mult_nodes << std::endl;
        std::cout << "(*MemEstMM: " << input_mesh << " *) ,{" << memMults.num_nodes << ", " << memMults.num_tet4 << "," << memMults.num_hex8 << "," << MegaByte(tot_mem)
                  << ", " << MegaByte(estMem) << "}" << std::endl;

      }
    if (estimate_memory_usage && query_only)
      {
        //MemoryMultipliers::process_estimate(0, eMesh, refiner.getRefinementInfoByType(), memory_multipliers_file, input_mesh);
      }
  }

  void MeshAdapt::create_refine_pattern()
  {
    if (use_transition_elements)
      {
        eMeshP->register_and_set_refine_fields();
        std::set<stk::mesh::Part *> pl;
        if (1)
          {
            for (unsigned ii=0; ii < m_block_names[eMeshP->element_rank()].size(); ++ii)
              {
                std::string bn = m_block_names[eMeshP->element_rank()][ii];
                VERIFY_OP_ON ((bn[0] == '+' || bn[0] == '-'), ==, true, "bad block name: "+bn);
                std::string bname = bn.substr(1);
                stk::mesh::Part *part = eMeshP->get_fem_meta_data()->get_part(bname);

                if (!part) {
                    const std::string alias = eMeshP->get_ioss_mesh_data()->get_input_io_region()->get_alias(bname);

                    if (debug && !eMeshP->get_rank())
                      std::cout << "block= " << bname << " replaced with alias=" << alias << std::endl;

                    part = eMeshP->get_fem_meta_data()->get_part(alias);

                    const int mult = block_names_x_map[bname];
                    block_names_x_map[alias] = mult;
                }

                VERIFY_OP_ON(part, !=, 0, "couldn't find part: "+bname);

                if (bn[0] == '+')
                  {
                    pl.insert(part);
                  }
              }
            for (unsigned ii=0; ii < m_block_names[eMeshP->element_rank()].size(); ++ii)
              {
                std::string bn = m_block_names[eMeshP->element_rank()][ii];
                std::string bname = bn.substr(1);
                stk::mesh::Part *part = eMeshP->get_fem_meta_data()->get_part(bname);

                if (bn[0] == '-')
                  {
                    if (pl.find(part) != pl.end())
                      pl.erase(part);
                  }
              }
          }
        stk::mesh::PartVector pv(pl.begin(), pl.end());
        block_selector = stk::mesh::selectUnion( pv );

        if (debug && !eMeshP->get_rank())
          std::cout << "block_names= " << m_block_names << "\nfrom parts = " << eMeshP->print_part_vector_string(pv)
                    << "\nblock_selector= " << block_selector << std::endl;

        pattern = Teuchos::get_shared_ptr( make_local_break_pattern(*eMeshP) );
      }
    else
      {
        pattern = Teuchos::get_shared_ptr(UniformRefinerPatternBase::createPattern(refine, enrich, convert, *eMeshP, m_block_names));
      }
  }

  void MeshAdapt::create_refiner()
  {
    if (use_transition_elements)
      {
        univ_selector.reset(new stk::mesh::Selector(eMeshP->get_fem_meta_data()->universal_part()));
        element_refine_predicate.reset(new ElementRefinePredicate(*eMeshP, univ_selector.get(), eMeshP->m_refine_field, 0.0));
        refiner.reset(new TransitionElementAdapter<ElementRefinePredicate>(*element_refine_predicate, *eMeshP, *pattern, 0));
        refiner->setRemoveOldElements(false);
        refiner->setAlwaysInitializeNodeRegistry(false);
      }
    else
      {
        refiner.reset(new UniformRefiner(*eMeshP, *pattern, proc_rank_field_ptr));
        refiner->setRemoveOldElements(remove_original_elements);
        if (block_name_inc.size())
          {
            if (!eMeshP->get_rank()) std::cout << "WARNING: setting remove_original_elements=0 and output_active_elements_only=1 since --blocks was specified" << std::endl;
            refiner->setRemoveOldElements(false);
            refiner->setAlwaysInitializeNodeRegistry(false);
            eMeshP->output_active_children_only(true);
          }
      }
  }

  int MeshAdapt::bucket_in_block_names(stk::mesh::Bucket& bucket, std::string& block)
  {
    stk::mesh::PartVector pv = bucket.supersets();
    for (unsigned ii=0; ii < pv.size(); ++ii)
      {
        std::string partName = pv[ii]->name();

        stk::mesh::Part& part = *pv[ii];
        bool auto_part = 0 != part.attribute<AutoPart>();
        if (stk::mesh::is_auto_declared_part(part) || auto_part)
          continue;

        if (block_names_x_map.find(partName) != block_names_x_map.end())
          {
            int mult = block_names_x_map[partName];
            block = partName;
            return mult;
          }
      }
    block = "";
    return 0;
  }

  void MeshAdapt::pre_refiner_tasks(int iBreak)
  {
    if (use_transition_elements)
      {
        const stk::mesh::BucketVector & buckets = eMeshP->get_bulk_data()->buckets( eMeshP->element_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            stk::mesh::Bucket & bucket = **k ;
            RefineFieldType_type val = static_cast<RefineFieldType_type>(0);
            if (block_selector(bucket))
              {
                std::string pname;
                int mult = bucket_in_block_names(bucket, pname);
                if (iBreak < mult)
                  {
                    val = static_cast<RefineFieldType_type>(1);
                  }
              }
            const unsigned num_elements_in_bucket = bucket.size();
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity element = bucket[iElement];
                RefineFieldType_type *fdata = stk::mesh::field_data( *eMeshP->m_refine_field , element );
                fdata[0] = static_cast<RefineFieldType_type>(val);
              }
          }
      }
    else if (block_name_inc.size())
      {
        stk::mesh::PartVector pv = eMeshP->get_fem_meta_data()->get_mesh_parts();
        typedef std::set<stk::mesh::Part *> PartSet;
        PartSet parts_all,  parts_in, parts_not_in;
        for (unsigned ii=0; ii < pv.size(); ++ii)
          {
            bool auto_part = (0 != pv[ii]->attribute<AutoPart>());
            if (auto_part)
              continue;

            if (pv[ii]->primary_entity_rank() == eMeshP->element_rank())
              {
                parts_all.insert(pv[ii]);
              }
          }
        for (StringIntMap::iterator it = block_names_x_map.begin(); it != block_names_x_map.end(); ++it)
          {
            std::string partName = it->first;
            int mult = it->second;
            stk::mesh::Part *part = eMeshP->get_fem_meta_data()->get_part(partName);
            if (debug) std::cout << "part= " << part->name() << " iBreak= " << iBreak << " mult= " << mult << std::endl;
            if (iBreak < mult)
              {
                parts_in.insert(part);
              }
          }
        std::set_difference(parts_all.begin(), parts_all.end(), parts_in.begin(), parts_in.end(), std::inserter(parts_not_in, parts_not_in.begin()));
        pv.assign(parts_not_in.begin(), parts_not_in.end());
        if (debug && !eMeshP->get_rank())
          {
            std::cout << " without use_transition_elements excluded parts= " << eMeshP->print_part_vector_string(pv,"\n") << std::endl;
          }
        if (pv.size())
          refiner->setExcludeParts(pv);
      }
  }

  // init these static vars before we use them in the function below

  double RunAdaptRunInfo::mMeshInputTime = 0.0;
  double RunAdaptRunInfo::mMeshOutputTime = 0.0;
  double RunAdaptRunInfo::mAdaptTimeOverall = 0.0;
  double RunAdaptRunInfo::mAdaptCPUTimeOverall = 0.0;

  void MeshAdapt::run_adapt_run()
  {

    // to support legacy uses, added an extension to the output_mesh name if there is none
    if ( !(has_suffix(output_mesh, ".e") || has_suffix(output_mesh, ".exo"))) {
      output_mesh = output_mesh + ".e";
    }

    // (this will already be taken care of when using the simple form of user input)
    if (next_adapted_mesh == "") {
      size_t dot_index = output_mesh.find_last_of(".");
      std::string basename(output_mesh, 0, dot_index);
      std::string suffix = output_mesh.substr(dot_index+1);
      next_adapted_mesh = basename + "_ft." + suffix;
    }

    RunAdaptRunInfo::run_adapt_run(previous_adapted_mesh, input_mesh,
                                   next_adapted_mesh, output_mesh,
                                   adapt_input, input_geometry, smooth_geometry,
                                   ioss_read_options,
                                   ioss_write_options,
                                   property_map,
                                   verify_meshes,
                                   false);

    mMeshInputTime = RunAdaptRunInfo::mMeshInputTime;
    mMeshOutputTime = RunAdaptRunInfo::mMeshOutputTime;
    mAdaptTimeOverall = RunAdaptRunInfo::mAdaptTimeOverall;
    mAdaptCPUTimeOverall = RunAdaptRunInfo::mAdaptCPUTimeOverall;
  }

  void MeshAdapt::do_dihedral_angle_check()
  {
    m_dihedral_angle_check->find_simplex_elements_with_obtuse_angles();
  }

  BlockNamesType MeshAdapt::process_block_names()
  {
    // FIXME move this next block of code to a method on UniformRefiner
    BlockNamesType block_names(percept::EntityRankEnd+1u);
    BlockNamesType block_names_rbar(percept::EntityRankEnd+1u);

#if defined(STK_BUILT_IN_SIERRA)
    if (rbar_blocks.length())
      {
        BlockNamesType rbar_names(percept::EntityRankEnd+1u);
        std::string block_name_inc_orig = block_name_inc;
        block_names_rbar = RefinerUtil::getBlockNames(block_name_inc, eMeshP->get_rank(), *eMeshP);
        if (rbar_blocks.length())
          {
            rbar_names = RefinerUtil::getBlockNames(rbar_blocks, eMeshP->get_rank(), *eMeshP);
            if (!eMeshP->get_rank())
              std::cout << "rbar_blocks= " << rbar_blocks << " rbar_names= " << rbar_names << std::endl;
          }
        for (unsigned ii=0; ii < rbar_names[eMeshP->element_rank()].size(); ii++)
          {
            std::string srb = rbar_names[eMeshP->element_rank()][ii];
            Util::replace(srb, "+", "-");
            block_names_rbar[eMeshP->element_rank()].push_back(srb);
          }
        block_name_inc = "";
        for (unsigned ii=0; ii < block_names_rbar[eMeshP->element_rank()].size(); ii++)
          {
            block_name_inc += (ii ? "," : "") + block_names_rbar[eMeshP->element_rank()][ii];
          }
        if (!eMeshP->get_rank())
          std::cout << "rbar: original block_name option = " << block_name_inc_orig << " new = " << block_name_inc << std::endl;
      }
#endif

    if (block_name_inc.length())
      {
        if (debug && !eMeshP->get_rank()) std::cout << "block_names: original (or after rbar) block_name option = " << block_name_inc << std::endl;
        block_names = RefinerUtil::getBlockNames(block_name_inc, eMeshP->get_rank(), *eMeshP);

        if (1)
          {
            int number_refines_max = 0;
            for (unsigned ii=0; ii < block_names[eMeshP->element_rank()].size(); ++ii)
              {
                std::string& bn = block_names[eMeshP->element_rank()][ii];
                //std::cout << "bn= " << bn << std::endl;
                size_t pos = bn.find(":");
                if (pos != std::string::npos)
                  {
                    std::string xp = bn.substr(pos+1);
                    size_t xpos = xp.find("x");
                    if (xpos == std::string::npos)
                      xpos = xp.find("X");
                    VERIFY_OP_ON(xpos, !=, std::string::npos, "syntax for multiplier for --blocks option is missing 'x' or 'X'");
                    xp = xp.substr(0, xpos);
                    std::string bnNew = bn.substr(0,pos);
                    if (!eMeshP->get_rank())
                      std::cout << "Note: found multiplier in block name, bn = " << bn << " mult= " << xp << " bnNew= " << bnNew << std::endl;
                    int mult = toInt(xp);
                    size_t posP = bnNew.find("+");
                    size_t posM = bnNew.find("-");
                    if (posP == std::string::npos && posM == std::string::npos)
                      {
                        VERIFY_MSG("ERROR: block name after processing not preceded by '+' or '-', bnNew= " +bnNew);
                      }
                    if (posP != std::string::npos)
                      {
                        bn = bnNew;
                        bnNew = bnNew.substr(1);
                        block_names_x_map[bnNew] = mult;
                        number_refines_max = std::max(number_refines_max, mult);
                        if (debug && !eMeshP->get_rank())
                          std::cout << "Note: bnNew= " << bnNew << std::endl;
                      }
                  }
                else
                  {
                    size_t posP = bn.find("+");
                    size_t posM = bn.find("-");
                    if (posP == std::string::npos && posM == std::string::npos)
                      {
                        VERIFY_MSG("ERROR: block name after processing not preceded by '+' or '-', bn= " +bn);
                      }
                    if (posP != std::string::npos)
                      {
                        std::string bnNew = bn;
                        bnNew = bnNew.substr(1);
                        block_names_x_map[bnNew] = number_refines;
                        if (debug && !eMeshP->get_rank())
                          std::cout << "Note: bnNew= " << bnNew << std::endl;
                      }
                  }
              }
            if (number_refines < number_refines_max)
              {
                if (!eMeshP->get_rank())
                  std::cout << "WARNING: number_refines must be set >= to the max found in --blocks option with the optional :Nx specification = " << number_refines_max
                            << "\n  NOTE: resetting number_refines to the max value." << std::endl;
              }
            if (debug && !eMeshP->get_rank())
              {
                std::cout << "block_names_x_map=\n" << block_names_x_map << std::endl;
              }
          }

        if (1)
          {
            eMeshP->commit();
            block_names = RefinerUtil::correctBlockNamesForPartPartConsistency_1(*eMeshP, block_names, input_geometry);

            eMeshP->close();
            pre_open();
            eMeshP->open(input_mesh);

            if (smooth_geometry) {
                const bool output = smooth_geometry > 1;
                eMeshP->add_coordinate_state_fields(output);
            }
#if !defined(NO_GEOM_SUPPORT)
            if (respect_spacing >= 1) {
              const bool output = respect_spacing > 1;
              eMeshP->set_respect_spacing(true);
              eMeshP->add_spacing_fields(output);
            }
#endif
            if (smooth_surfaces == 1) eMeshP->set_smooth_surfaces(true);

          }
        if (debug && !eMeshP->get_rank()) std::cout << "block_names after processing: " << block_names << std::endl;
      }
    return block_names;
  }


  void MeshAdapt::pre_open()
  {
    eMeshP->set_avoid_add_all_mesh_fields_as_input_fields(true);
  }

  int MeshAdapt::do_run_pre_commit()
  {
    if (load_balance)
      {
        throw std::runtime_error("deprecated - run decomp manually before running mesh_adapt, or, use the auto-decomp options");
      }

    if (ioss_read_options.length() || ioss_write_options.length())
      {
        if (!eMeshP->get_rank())
          {
            std::cout << "INFO: ioss_read_options= " << ioss_read_options << " ioss_write_options= " << ioss_write_options << std::endl;
          }
      }

    if (smooth_geometry)
      {
        if (smooth_surfaces == 1) eMeshP->set_smooth_surfaces(true);
        eMeshP->setProperty("smoother_niter", toString(smoother_niter));
        eMeshP->setProperty("smoother_tol", toString(smoother_tol));
        eMeshP->setProperty("smooth_use_reference_mesh", (smooth_use_reference_mesh?"1":"0"));
      }
    if (ioss_read_options.length())  eMeshP->set_ioss_read_options(ioss_read_options);
    if (ioss_write_options.length()) eMeshP->set_ioss_write_options(ioss_write_options);

    if (adapt_input.length()) {
      run_adapt_run();
      return -1;
    }

    eMeshP->set_do_print_memory(print_memory_usage);

    pre_open();
    {
      double t0 = stk::wall_time();
      if (generated_mesh)
        eMeshP->new_mesh(GMeshSpec(input_mesh));
      else
        eMeshP->open(input_mesh);
      double t1 = stk::wall_time();
      mMeshInputTime = t1 - t0;
    }
    if (smooth_surfaces == 1) eMeshP->set_smooth_surfaces(true);

    if (smooth_geometry) {
        const bool output = smooth_geometry > 1;
        eMeshP->add_coordinate_state_fields(output);
    }
#if !defined(NO_GEOM_SUPPORT)
    if (respect_spacing >= 1) {
      const bool output = respect_spacing > 1;
      eMeshP->set_respect_spacing(true);
      eMeshP->add_spacing_fields(output);
    }
#endif
    if (sync_io_regions)
      {
        if (convert.length() || enrich.length())
          {
            if (!eMeshP->get_rank())
              {
                std::cout       << "WARNING: removing property original_topology_type from input (and output mesh) since topology is changed by convert or enrich" << std::endl;
              }
            eMeshP->set_remove_io_orig_topo_type(true);
          }
      }
    eMeshP->set_sync_io_regions(sync_io_regions);
    if (!s_spatialDim) s_spatialDim = eMeshP->get_spatial_dim();

    Util::setRank(eMeshP->get_rank());

    if (doRefineMesh && ((input_geometry_type != PGEOM_ACIS) && (input_geometry_type != PGEOM_OPENNURBS)) )
    {
      m_block_names = process_block_names();
      create_refine_pattern();
    }

    int scalarDimension = 0; // a scalar

    if (proc_rank_field)
      {
        proc_rank_field_ptr = eMeshP->add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
      }

    if (fix_all_block_boundaries)
      {
        bool make_part_io_part=true;
        eMeshP->add_part("inner_skin_part", make_part_io_part);
      }

    // Mesh-based geometry Fitting - setup parts
    mesh_based_geometry_setup();

    if (input_geometry.length()) {
      eMeshP->register_and_set_smoothing_fields(); 

      setup_m2g_parts(input_geometry);
    }

    do_precheck_memory_usage();

    eMeshP->add_registered_refine_fields_as_input_fields();

    eMeshP->get_ioss_mesh_data()->add_all_mesh_fields_as_input_fields();

    if (dihedral_angle_check != 0)
      {
        m_dihedral_angle_check.reset(new DihedralAngleCheck(eMeshP.get(), (dihedral_angle_check > 0 ? dihedral_angle_check : 1)));
      }
    return 0;
  }

  int MeshAdapt::do_run_post_commit()
  {
    if (use_side_map)
      {
        eMeshP->setProperty("use_side_map", "true");
      }

    if (DO_MEMORY) {
      std::string hwm = print_memory_both(eMeshP->parallel());
      if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " initial memory after opening input mesh."  << std::endl;
    }

    if (fix_all_block_boundaries)
      {
        eMeshP->get_skin_part("inner_skin_part", true);
      }

    if (print_info)
      {
        eMeshP->print_info("PerceptMesh info:", print_info);
      }

    if (dihedral_angle_check != 0)
      {
        do_dihedral_angle_check();
      }

    if (skin_mesh)
      {
        eMeshP->skin_mesh();
      }

    if (convert_geometry_parts_OpenNURBS.size())
      {
        if (!input_geometry.size())
          throw std::runtime_error("if convert_geometry_parts_OpenNURBS is specified you must specify an input geometry --input_geometry=<file>");
        do_convert_geometry_parts_OpenNURBS(input_geometry, convert_geometry_parts_OpenNURBS);
      }

    if (input_geometry_type == MESH_BASED_GREGORY_PATCH)
      {
        stk::mesh::FieldVector from(1, eMeshP->get_coordinates_field()), to(1, eMeshP->m_unprojected_coordinates);
        PerceptMesh::copy_fields(*eMeshP, *eMeshP, from, to, eMeshP->node_rank(), EntitySelectorUCF(*eMeshP));
      }

    // Mesh-based-geometry fitting
    mesh_based_geometry_fitting();

    verify_mesh_util(true);

    do_histograms_init();

    compute_hmesh_sizes_init();

    if (print_hmesh_surface_normal)
      {
        std::string msg="before refine";
        eMeshP->print_hmesh_surface_normal(msg, std::cout);
      }

#if !defined(NO_GEOM_SUPPORT)
    if (respect_spacing)
      {
        SpacingFieldUtil sfu(*eMeshP);
        sfu.compute_spacing_field();
      }
#endif

    if (print_memory_usage)
      memory_dump(print_memory_usage, eMeshP->get_bulk_data()->parallel(), *eMeshP->get_bulk_data(), 0, "after file open");

    if (test_memory_nodes && test_memory_elements)
      {
        do_test_memory();
        return 0;
      }

    if (doRefineMesh)
      {
        t0 =  stk::wall_time();
        cpu0 = stk::cpu_time();

        create_refiner();

#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 0
        ProgressMeter pm(*refiner);
#endif

        if (input_geometry != "") {
          refiner->setGeometryFile(input_geometry);
          refiner->setSmoothGeometry(smooth_geometry);
          
          initialize_m2g_geometry(input_geometry);
        }

        refiner->setFixAllBlockBoundaries(fix_all_block_boundaries);
        refiner->setQueryPassOnly(query_only == 1);
        refiner->setDoProgressMeter(progress_meter == 1 && 0 == p_rank);
#if defined(STK_BUILT_IN_SIERRA)
        if (rbar_blocks.length())
          {
            BlockNamesType rbar_names(percept::EntityRankEnd+1u);
            if (rbar_blocks.length())
              {
                if (!eMeshP->get_rank())
                  std::cout << "For RBAR treatment: rbar_blocks= " << rbar_blocks << std::endl;

                rbar_names = RefinerUtil::getBlockNames(rbar_blocks, eMeshP->get_rank(), *eMeshP);
                if (!eMeshP->get_rank())
                  std::cout << "For RBAR treatment: rbar_names (after getBlockNames)= " << rbar_names << std::endl;
              }
            refiner->set_rbar_special_treatment(rbar_names);
          }
#endif

        for (int iBreak = 0; iBreak < number_refines; iBreak++)
          {
            if (!eMeshP->get_rank())
              {
                std::cout << "Refinement pass # " << (iBreak+1) << " start..." << std::endl;
              }

            if (save_intermediate_meshes < 0 && iBreak == 0)
              {
                eMeshP->save_as("refined_mesh.e");
              }
            if (save_intermediate_meshes > 0 && iBreak == 0)
              {
                eMeshP->save_as("refined_mesh_0.e");
              }

            pre_refiner_tasks(iBreak);


            refiner->doBreak();

            if (!eMeshP->get_rank())
              {
                std::cout << std::endl;
                int ib = iBreak;
                if (!query_only) ib = 0;
                bool printAllTopologies = false;
                refiner->getRefinementInfo().printTable(std::cout, ib , printAllTopologies);
                std::cout << std::endl;
              }
            if (print_memory_usage)
              {
                memory_dump(print_memory_usage, eMeshP->get_bulk_data()->parallel(), *eMeshP->get_bulk_data(), &refiner->getNodeRegistry(),
                            std::string("after refine pass: ")+toString(iBreak));
              }

            if (estimate_memory_usage)
              {
                std::cout << "memory_multipliers_file= " << memory_multipliers_file << std::endl;
                if (query_only)
                  {
                    if (number_refines == 1)
                      {
                        MemorySizeType tot_mem = memory_dump(false, eMeshP->get_bulk_data()->parallel(), *eMeshP->get_bulk_data(), &refiner->getNodeRegistry(),
                                                             std::string("after mesh read"));
                        if (!p_rank)
                          std::cout << "P[" << p_rank << "] tmp srk mem after mesh read= " << MegaByte(tot_mem) << std::endl;
                        bool use_new = false;
                        MemoryMultipliers::process_estimate(tot_mem, *eMeshP, refiner->getRefinementInfo(), memory_multipliers_file, input_mesh, use_new);
                      }

                    refiner->getRefinementInfo().estimateNew(iBreak);
                    MemoryMultipliers::process_estimate(0, *eMeshP, refiner->getRefinementInfo(), memory_multipliers_file, input_mesh);
                  }
                else
                  {
                    MemorySizeType tot_mem = memory_dump(false, eMeshP->get_bulk_data()->parallel(), *eMeshP->get_bulk_data(), &refiner->getNodeRegistry(),
                                                         std::string("after refine pass: ")+toString(iBreak));
                    if (!p_rank)
                      std::cout << "P[" << p_rank << "] tmp srk tot_mem= " << MegaByte(tot_mem) << std::endl;
                    MemoryMultipliers::process_estimate(tot_mem, *eMeshP, refiner->getRefinementInfo(), memory_multipliers_file, input_mesh);
                  }
              }
            if (save_intermediate_meshes < 0)
              {
                char buf[1000];
                sprintf(buf, "%04d", iBreak+1);
                eMeshP->save_as("refined_mesh.e-s"+std::string(buf));
              }
            if (save_intermediate_meshes > 0)
              {
                eMeshP->save_as("refined_mesh_"+toString(iBreak+1)+".e");
              }

            if (DO_MEMORY) {
              std::string hwm = print_memory_both(eMeshP->parallel());
              if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " memory after refining mesh for pass # " << iBreak  << std::endl;
            }


          } // iBreak


        if (delete_parents)
          {
            if (DO_MEMORY) {
              std::string hwm = print_memory_both(eMeshP->parallel());
              if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " memory before deleteParentElements" << std::endl;
            }

            refiner->deleteParentElements();

            if (DO_MEMORY) {
              std::string hwm = print_memory_both(eMeshP->parallel());
              if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " memory after deleteParentElements"  << std::endl;
            }
          }


        if (number_refines == 0 && (smooth_geometry || snap_geometry))
          {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
            refiner->setSmoothGeometry(smooth_geometry);
            refiner->snapAndSmooth((snap_geometry == 1), input_geometry, (smooth_use_reference_mesh == 1) );
#elif defined(NO_GEOM_SUPPORT)
            std::ostringstream oss;
            oss << "\nERROR: Geometry and/or smoothing is not currently supported on this platform. Try running with geometry turned off.";
            throw std::runtime_error(oss.str());
#endif
          }

        if (delete_parents)
          refiner->deleteParentElements();

        t1 =  stk::wall_time();
        cpu1 = stk::cpu_time();

        mAdaptTimeOverall = t1 - t0;
        mAdaptCPUTimeOverall = cpu1 - cpu0;

        verify_mesh_util(false);

        do_histograms_final();

        compute_hmesh_sizes_final();

        if (print_hmesh_surface_normal)
          {
            std::string msg="after refine";
            eMeshP->print_hmesh_surface_normal(msg, std::cout);
          }

        if (DEBUG_ADAPT_MAIN || 0 == p_rank) {
          std::cout << "P[" << p_rank << "]  AdaptMain:: saving mesh... " << std::endl;
        }
        if (remove_geometry_blocks) eMeshP->remove_geometry_blocks_on_output(input_geometry);

        {
          double t0 = stk::wall_time();
          eMeshP->save_as(output_mesh);
          double t1 = stk::wall_time();
          mMeshOutputTime = t1 - t0;
        }
        if (DEBUG_ADAPT_MAIN || 0 == p_rank) {
          std::cout << "P[" << p_rank << "]  AdaptMain:: mesh saved" << std::endl;
        }

        if (print_memory_usage)
          memory_dump(print_memory_usage, eMeshP->get_bulk_data()->parallel(), *eMeshP->get_bulk_data(), &refiner->getNodeRegistry(), "after final save mesh");

      } // doRefineMesh
    return 0;
  }

  int MeshAdapt::do_run_post_refine()
  {

    return 0;
  }

  int MeshAdapt::do_post_proc()
  {

#if 0
    for (int itime=0; itime < 10; itime++)
      {
        std::cout << "tmp timer[" << itime << "]= " << s_timers[itime] << " " << s_timers[itime]/s_timers[3]*100 << " %" << std::endl;
      }
#endif
    if (DO_MEMORY) {
      std::string hwm = print_memory_both(eMeshP->parallel());
      if (!eMeshP->get_rank()) std::cout << "MEM: " << hwm << " final memory after refining mesh."  << std::endl;
    }


    //stk::all_reduce( *eMeshP->get_bulk_data()->parallel(), stk::ReduceSum<1>( &failed_proc_rank ) );
    if (failed_proc_rank)
      {
        std::cout << "P[" << p_rank << "] AdaptMain::exception found on processor " << (failed_proc_rank-1) << std::endl;
        exit(1);
      }

    if (DEBUG_ADAPT_MAIN)
      {
        std::cout << "P[" << p_rank << ", " << p_size << "]  wall clock time on processor [" << p_rank << ", " << p_size << "]= " << (t1-t0) << " (sec) "
                  << " cpu time= " << (cpu1 - cpu0) << " (sec) " << std::endl;
      }

    double cpuMax = (cpu1-cpu0);
    double wallMax = (t1-t0);
    double wallMin = wallMax;
    double cpuSum = (cpu1-cpu0);
    double cpuMin = cpuMax;

    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceSum<1>( &cpuSum ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMax<1>( &cpuMax ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMin<1>( &cpuMin ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMax<1>( &wallMax ) );
    stk::all_reduce( eMeshP->get_bulk_data()->parallel(), stk::ReduceMin<1>( &wallMin ) );

    if (0 == p_rank)
      {
        std::cout << "P[" << p_rank << ", " << p_size << "]  min wall clock time = " << wallMin << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  max wall clock time = " << wallMax << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  min cpu  clock time = " << cpuMin << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  max cpu  clock time = " << cpuMax << " (sec)" << std::endl;
        std::cout << "P[" << p_rank << ", " << p_size << "]  sum cpu  clock time = " << cpuSum << " (sec)" << std::endl;
      }

    if (print_timers)
      {
        std::ostringstream str;
        refiner->rootTimer().stop();
        stk::diag::printTimersTable(str, refiner->rootTimer(),
                                    stk::diag::METRICS_ALL, false,
                                    eMeshP->get_bulk_data()->parallel() );
        if (0 == p_rank)
          {
            std::cout << str.str() << std::endl;
          }
      }

    return 0;
  }

  // version - return true if build in Sierra
  bool MeshAdapt::get_version(std::string* v)
  {
    if (v) *v = version_prefix+version;
#if defined(STK_BUILT_IN_SIERRA)
    return true;
#else
    return false;
#endif
  }

  void MeshAdapt::log_usage( bool status )
  {
#if defined(STK_BUILT_IN_SIERRA)
    const bool disable_audit = !sierra::Env::get_param("noaudit").empty() || std::getenv("SIERRA_USAGE_METRICS_OFF") != NULL;

    size_t hwm_max = 0, hwm_min = 0, hwm_avg = 0;

    stk::get_memory_high_water_mark_across_processors(eMeshP->parallel(), hwm_max, hwm_min, hwm_avg);

    stk::all_reduce( eMeshP->parallel(), stk::ReduceSum<1>( &mAdaptCPUTimeOverall ) );

    if (eMeshP->get_parallel_rank() == 0) {
      // Audit log
      bool runtest = sierra::Env::get_param("runtest").empty() ? false : true;
      if ((!disable_audit) && (!runtest)) {
        const double bytes_in_MB = 1024*1024;
        auditdata data;
        AuditLogDefaults(&data, "mesh_adapt", sierra::ProductRegistry::version(), eMeshP->get_parallel_size());
        strcpy(data.starttime,sierra::format_time(sierra::Env::start_time(), "%Y%m%d%H%M%S").c_str());
        strcpy(data.purpose,"meshing");

        data.num_proc       = eMeshP->get_parallel_size();

        data.time_elapsed   = mAdaptTimeOverall;
        data.cputimesum     = mAdaptCPUTimeOverall;
        data.meshinputtime  = mMeshInputTime;
        data.meshoutptutime = mMeshOutputTime;
        data.hwm_min        = hwm_min / bytes_in_MB;
        data.hwm_max        = hwm_max / bytes_in_MB;
        data.hwm_avg        = hwm_avg / bytes_in_MB;
        strcpy(data.status, status ? "success" : "fail" );

#if 0
        printf( "product %s\n", data.product );
        printf( "version %s\n", data.version );
        printf( "purpose %s\n", data.purpose );
        printf( "hostname %s\n", data.hostname );
        printf( "username %s\n", data.username );
        printf( "num_proc %d\n", data.num_proc );
        printf( "starttime %s\n", data.starttime );
        printf( "time_elapsed %lf\n", data.time_elapsed );
        printf( "cputimesum %lf\n", data.cputimesum );
        printf( "meshinputtime %lf\n", data.meshinputtime );
        printf( "meshoutputtime %lf\n", data.meshoutptutime );
        printf( "hwm_max %lf\n", data.hwm_max );
        printf( "hwm_min %lf\n", data.hwm_min );
        printf( "hwm_avg %lf\n", data.hwm_avg );
        printf( "status %s\n", data.status );
        //printf( "audit_filename %s\n", audit_filename.c_str() );
#endif
        OutputAuditLog(&data);
      }
    }
#endif    
  }

  int MeshAdapt::adapt_main_full_options(int argc, char **argv)
  {
#if defined( STK_HAS_MPI )
    stk::ParallelMachine comm(stk::parallel_machine_init(&argc, &argv));
#endif
#if defined(WITH_KOKKOS)
    Kokkos::initialize(argc, argv);
#endif
    EXCEPTWATCH;

    Teuchos::CommandLineProcessor clp;
    p_rank = stk::parallel_machine_rank(comm);
    p_size = stk::parallel_machine_size(comm);

    setup_options(clp, comm, argc, argv);

    return adapt_main_full_options_normal(argc, argv);
  }

  int MeshAdapt::adapt_main_full_options_normal(int argc, char **argv)
  {
    t0   = 0.0;
    t1   = 0.0;
    cpu0 = 0.0;
    cpu1 = 0.0;

#if defined( STK_PERCEPT_HAS_GEOMETRY )
    if (dump_geometry_file && input_geometry.find(".3dm") != std::string::npos)
      {
        GeometryKernelOpenNURBS gko;
        gko.debug_dump_file(input_geometry);
      }
#endif

    std::string input_mesh_save = input_mesh;
    std::string output_mesh_save = output_mesh;

    eMeshP.reset(new percept::PerceptMesh);
    if (output_active_elements_only)
      eMeshP->output_active_children_only(true);

    if (1)
      {
        if (property_map.length())
          {
            eMeshP->parse_property_map_string(property_map);
          }
        const char * env_val = std::getenv("Percept_property_map");
        if (env_val)
          {
            std::string vv(env_val);
            if (eMeshP->get_rank() == 0) std::cout << "found Percept_property_map = " << vv << std::endl;
            eMeshP->parse_property_map_string(vv);
          }
        if (eMeshP->getProperty("MeshAdapt.debug") == "true")
          debug = 1;
      }

    if (progress_meter && eMeshP->get_rank() == 0)
      {
        std::cout << "Stage: Open mesh..." << " cpu: " << eMeshP->cpu_time() << " [sec]" << std::endl;
      }

    int res1 = do_run_pre_commit();

    // trap for RunAdaptRun, etc.
    if (res1 < 0)
      {
        log_usage();

        eMeshP.reset();

#if defined( STK_HAS_MPI )
        stk::parallel_machine_finalize();
#endif
#if defined(WITH_KOKKOS)
        Kokkos::finalize();
#endif
        return 0;
      }

    if (progress_meter && eMeshP->get_rank() == 0)
      {
        std::cout << "Stage: Commit mesh..." << " cpu: " << eMeshP->cpu_time() << " [sec]" << std::endl;
      }

    eMeshP->commit();
    if (progress_meter && eMeshP->get_rank() == 0)
      {
        std::cout << "Stage: Commit mesh...done" << " cpu: " << eMeshP->cpu_time() << " [sec]" << std::endl;
      }

    do_run_post_commit();

    do_post_proc();

    log_usage();

    eMeshP.reset();

#if defined( STK_HAS_MPI )
    stk::parallel_machine_finalize();
#endif
#if defined(WITH_KOKKOS)
    Kokkos::finalize();
#endif
    return result;
  }


  /// removes beams and shells in favor of node sets
  void MeshAdapt::do_convert_geometry_parts_OpenNURBS(std::string geometry_file, std::string newExodusFile)
  {
#if defined(STK_PERCEPT_HAS_GEOMETRY)

    eMeshP->reopen();

    GeometryKernelOpenNURBS kernel;
    MeshGeometry mesh_geometry(*eMeshP, &kernel);
    GeometryFactory factory(&kernel, &mesh_geometry);
    factory.read_file(geometry_file, &*eMeshP);
    const std::vector<GeometryEvaluator*>& evals = mesh_geometry.getGeomEvaluators();
    stk::mesh::PartVector newParts, oldParts;
    for (unsigned ii=0; ii < evals.size(); ++ii)
      {
        stk::mesh::Part *part = evals[ii]->mPart;
        oldParts.push_back(part);
        std::string partName = part->name();
        std::string newPartName = partName + "_n";
        stk::mesh::Part& newPart = eMeshP->get_fem_meta_data()->declare_part(newPartName, stk::topology::NODE_RANK);
        stk::io::put_io_part_attribute(newPart);
        stk::io::remove_io_part_attribute(*part);
        newParts.push_back(&newPart);
      }
    eMeshP->commit();

    eMeshP->get_bulk_data()->modification_begin();
    for (unsigned ii=0; ii < newParts.size(); ++ii)
      {
        stk::mesh::Part& newPart = *newParts[ii];
        stk::mesh::Part& oldPart = *oldParts[ii];
        stk::mesh::PartVector addParts(1, &newPart);

        if (0) std::cout << "oldPart= " << oldPart.name() << " newPart= " << newPart.name() << std::endl;

        stk::mesh::Selector sp(oldPart);
        std::vector<stk::mesh::Entity> vec;
        //stk::mesh::Selector on_locally_owned_part =  ( eMeshP->get_fem_meta_data()->locally_owned_part() );
        stk::mesh::get_selected_entities(sp, eMeshP->get_bulk_data()->buckets(eMeshP->element_rank()), vec);
        for (size_t jj=0; jj < vec.size(); ++jj)
          {
            stk::mesh::Entity element = vec[jj];
            const MyPairIterRelation elem_nodes(*eMeshP, element, eMeshP->node_rank() );
            for (unsigned kk=0; kk < elem_nodes.size(); ++kk)
              {
                stk::mesh::Entity node = elem_nodes[kk].entity();
                if (eMeshP->owned(node))
                  {
                    eMeshP->get_bulk_data()->change_entity_parts(node, addParts);
                  }
              }
          }
        for (size_t jj=0; jj < vec.size(); ++jj)
          {
            stk::mesh::Entity element = vec[jj];
            MyPairIterRelation elem_sides(*eMeshP, element, eMeshP->side_rank());
            for (unsigned kk=0; kk < elem_sides.size(); ++kk)
              {
                if ( ! eMeshP->get_bulk_data()->destroy_relation(element, elem_sides[kk].entity(), elem_sides[kk].relation_ordinal()))
                  {
                    throw std::logic_error("convert_geometry_parts_OpenNURBS failed to delete relation");
                  }
              }

            if ( ! eMeshP->get_bulk_data()->destroy_entity( element ) )
              {
                throw std::logic_error("convert_geometry_parts_OpenNURBS failed to delete element");
              }
          }
        if (0)
          {
            std::vector<stk::mesh::Entity> vecNodes;
            stk::mesh::get_selected_entities(stk::mesh::Selector(newPart), eMeshP->get_bulk_data()->buckets(eMeshP->node_rank()), vecNodes);
            std::cout << "part= " << newPart.name() << " vecNodes.size= " << vecNodes.size() << std::endl;
          }
      }
    if (1)
      {
        stk::mesh::PartVector excludeParts;
        SidePartMap side_part_map;
        std::string geomFile = geometry_file;
        bool avoidFixSideSetChecks = false;
        FixSideSets fss(0, *eMeshP, excludeParts, side_part_map, geomFile, avoidFixSideSetChecks);
        eMeshP->get_bulk_data()->modification_begin();
        fss.fix_side_sets_2(false, 0, 0, "MeshAdapt");
        stk::mesh::fixup_ghosted_to_shared_nodes(*eMeshP->get_bulk_data());
        eMeshP->get_bulk_data()->modification_end();
        eMeshP->save_as("junk.e");
      }

    stk::mesh::fixup_ghosted_to_shared_nodes(*eMeshP->get_bulk_data());
    eMeshP->get_bulk_data()->modification_end();
    //eMeshP->print_info("new",2);

    if (0)
      {
        for (unsigned ii=0; ii < newParts.size(); ++ii)
          {
            stk::mesh::Part& newPart = *newParts[ii];
            std::vector<stk::mesh::Entity> vecNodes;
            stk::mesh::get_selected_entities(stk::mesh::Selector(newPart), eMeshP->get_bulk_data()->buckets(eMeshP->node_rank()), vecNodes);
            std::cout << "part= " << newPart.name() << " vecNodes.size= " << vecNodes.size() << std::endl;
          }
      }


    eMeshP->save_as(newExodusFile);
#endif
  }

void MeshAdapt::setup_m2g_parts(std::string input_geometry)
{
  if ( (input_geometry_type != PGEOM_ACIS) && (input_geometry_type != PGEOM_OPENNURBS) )	  return;

#ifdef HAVE_CUBIT        
#ifdef HAVE_ACIS
  //Madison Brewer: what I don't like about this is that we're not really utilizing the geometry interface layer/kernel to its full extent.
  //It essentially just gets called for snapping and disappears after that.
  //BIG QUESTION: Do we want to try and refactor percept to truly interact with its geometry through these kernels? How much refactoring would this incur?

  m_PGA = new PGeomACIS;
  m_PGeomPntr = m_PGA;
#else
  m_PGeomPntr = new PGeom;
#endif

  if (input_geometry_type == PGEOM_ACIS) {
    m_PGeomPntr->initialize(ACIS_GEOMETRY_ENGINE);
    m_PGeomPntr->import_acis_file(input_geometry.c_str());
  }
  else if (input_geometry_type == PGEOM_OPENNURBS) {
    m_PGeomPntr->initialize(OPENNURBS_GEOMETRY_ENGINE);
    m_PGeomPntr->import_open_nurbs_file(input_geometry.c_str());
  }

  stk::mesh::MetaData * md = eMeshP->get_fem_meta_data();

  std::vector<int> surfIDs;
  m_PGeomPntr->get_surfaces(surfIDs);
  std::vector<std::string> quadNames(surfIDs.size());
  std::vector<std::string> triNames(surfIDs.size());

  //make parts that we'll store new mesh entities on
  for(unsigned i = 0; i < surfIDs.size(); i++) { 
    std::string name = "geom_surface_quad_";
    name = name + std::to_string(surfIDs[i]);
    stk::mesh::Part& part = md->declare_part_with_topology(name,
                                                           stk::topology::SHELL_QUAD_4);
    if (dump_geometry_file) stk::io::put_io_part_attribute(part);
    quadNames[i] = name;
  } 

  for(unsigned i = 0; i<surfIDs.size();i++){
    std::string name = "geom_surface_tri_";
    name = name + std::to_string(surfIDs[i]);
    stk::mesh::Part& part = md->declare_part_with_topology(name,
                                                           stk::topology::SHELL_TRI_3);
    if (dump_geometry_file) stk::io::put_io_part_attribute(part);
    triNames[i] = name;
  }

  std::vector<int> curveIDs;
  m_PGeomPntr->get_curves(curveIDs);
  std::vector<std::string> curveNames(curveIDs.size());
  for (unsigned i = 0; i < curveIDs.size(); i++) { 
    std::string name = "geom_curve_";
    name = name + std::to_string(curveIDs[i]);
    stk::mesh::Part& part = md->declare_part_with_topology(name,
                                                           stk::topology::BEAM_2);
    if (dump_geometry_file) stk::io::put_io_part_attribute(part);
    curveNames[i] = name;
  }

  //setup refinement: bcarnes: why is this here?
  m_block_names = process_block_names();
  create_refine_pattern();
#endif
}
  
void MeshAdapt::initialize_m2g_geometry(std::string input_geometry)
{
  if( (input_geometry_type != PGEOM_ACIS) && (input_geometry_type != PGEOM_OPENNURBS) ) return;
#ifdef HAVE_CUBIT  

  std::string m2gFile = input_geometry.substr(0,input_geometry.length()-3) + "m2g";

  int THIS_PROC_NUM = stk::parallel_machine_rank( MPI_COMM_WORLD);

  stk::mesh::MetaData* md = eMeshP->get_fem_meta_data();
  stk::mesh::BulkData* bd = eMeshP->get_bulk_data();

  std::vector<int> curveIDs;
  m_PGeomPntr->get_curves(curveIDs);

  std::vector<int> surfIDs;
  m_PGeomPntr->get_surfaces(surfIDs);

  eMeshP->initializeIdServer();

  PGeomAssoc<stk::mesh::BulkData, stk::mesh::Entity, stk::mesh::Entity,
    stk::mesh::Entity, stk::mesh::Entity> geom_assoc(m_PGeomPntr);
  geom_assoc.set_node_callback(get_node_from_id);
  geom_assoc.set_edge_callback(get_beam_from_ids);
  geom_assoc.set_face_callback(get_shell_from_ids);
  geom_assoc.set_elem_callback(get_hex_from_id);
  geom_assoc.set_validate_nodes_callback(validate_node_ownership);
  geom_assoc.set_validate_element_callback(validate_element_ownership);
  geom_assoc.set_mesh(bd);
  geom_assoc.set_fill_curve_and_surface_maps_during_import(false);
  const bool geometry_exists = true;
  geom_assoc.import_m2g_file(m2gFile.c_str(), geometry_exists);

  bd->modification_begin();

  for (unsigned i = 0; i < curveIDs.size(); i++) { //create beams and put them into corresponding curve parts

    std::vector<stk::mesh::Part *> add_parts_beams(1, static_cast<stk::mesh::Part*>(0));

    add_parts_beams[0] = md->get_part("geom_curve_" + std::to_string(curveIDs[i]));

    std::vector<std::vector<int>> edge_node_ids;
    geom_assoc.get_curve_edge_nodes(curveIDs[i], edge_node_ids);

    for (unsigned ii = 0; ii < edge_node_ids.size(); ii++) {

      bool toDeclare = true;

      std::vector<stk::mesh::EntityId> beamNodeIDs;
      std::vector<int> procsSharedTo;
      std::vector<stk::mesh::EntityKey> keysToCheck;
      int lowestRank = std::numeric_limits<int>::max();

      for (unsigned j = 0; j < edge_node_ids[ii].size(); j++)
        beamNodeIDs.push_back((stk::mesh::EntityId) edge_node_ids[ii][j]);

      for (unsigned j = 0; j < edge_node_ids[ii].size(); j++) {

        stk::mesh::Entity cur_node = bd->get_entity(stk::topology::NODE_RANK, beamNodeIDs[j]);

        stk::mesh::EntityKey key = bd->entity_key(cur_node);
        keysToCheck.push_back(key);
      }

      bd->shared_procs_intersection(keysToCheck, procsSharedTo);
      procsSharedTo.push_back(THIS_PROC_NUM); //find all processes that own or have this node shared to it
      for (size_t iii = 0; iii < procsSharedTo.size(); iii++) {
        if (procsSharedTo[iii] < lowestRank)
          lowestRank = procsSharedTo[iii]; //lowest ranking process is responsible for creation of this entity
        //		QUESTION: does this create a significant load imbalance for creation?
      }

      if (lowestRank != THIS_PROC_NUM)
        toDeclare = false;

      if (toDeclare) {

        stk::mesh::EntityId id2 = eMeshP->getNextId(stk::topology::ELEMENT_RANK);
        stk::mesh::declare_element(*eMeshP->get_bulk_data(),
                                   add_parts_beams, id2, beamNodeIDs);
      }
    }
  }

  std::vector<std::pair<std::vector<int>, int>> uncreatedEnts;
  for (unsigned i = 0; i < surfIDs.size(); i++) {
    std::vector<stk::mesh::Part *> add_parts_shells(1,static_cast<stk::mesh::Part*>(0));

    std::vector<std::vector<int>> face_node_ids;
    geom_assoc.get_surface_edge_nodes(surfIDs[i], face_node_ids);

    for (unsigned ii = 0; ii < face_node_ids.size(); ii++) {

      std::vector<stk::mesh::EntityId> shellNodeIDs;
      for (unsigned j = 0; j < face_node_ids[ii].size(); j++)
        shellNodeIDs.push_back((stk::mesh::EntityId) face_node_ids[ii][j]);

      if (shellNodeIDs.size() == 3)
        add_parts_shells[0] = md->get_part("geom_surface_tri_"
                                           + std::to_string(surfIDs[i])); //store these parts in a partvector for FASTER access
      else {
        add_parts_shells[0] = md->get_part("geom_surface_quad_"
                                           + std::to_string(surfIDs[i]));
      }

      bool toDeclare = true;

      int lowestRank = std::numeric_limits<int>::max();
      std::vector<stk::mesh::EntityKey> keysToCheck;
      std::vector<int> procsSharedTo;

      stk::mesh::Entity cur_node;
      for (unsigned j = 0; j < face_node_ids[ii].size(); j++) {

        cur_node = bd->get_entity(stk::topology::NODE_RANK, shellNodeIDs[j]);

        stk::mesh::EntityKey key = bd->entity_key(cur_node);
        keysToCheck.push_back(key);
      }

      bd->shared_procs_intersection(keysToCheck, procsSharedTo);
      procsSharedTo.push_back(THIS_PROC_NUM);//find all processes that either own or have these nodes shared to them
      for (size_t iii = 0; iii < procsSharedTo.size(); iii++) {
        if (procsSharedTo[iii] < lowestRank)
          lowestRank = procsSharedTo[iii]; //lowest ranking process is responsible for creation
        //		QUESTION: does this create a significant load imbalance for creation?
      }

      if (lowestRank != THIS_PROC_NUM) {
        toDeclare = false;
      }

      if (toDeclare) {
        stk::mesh::EntityId id2 = eMeshP->getNextId(stk::topology::ELEMENT_RANK);

        stk::mesh::declare_element(*bd, add_parts_shells,
                                   id2, shellNodeIDs);
      }
    }
  }
  bd->modification_end();

  geom_assoc.fill_curve_and_surface_maps();
  delete m_PGeomPntr; //bad things happen if you don't explicitly reallocate the memory here
#endif
}

  int MeshAdapt::main(int argc, char **argv) {

    int res=0;
    res = adapt_main(argc, argv);
    return res;
  }

}

