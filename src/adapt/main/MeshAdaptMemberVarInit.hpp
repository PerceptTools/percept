// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


// Note: this file is included in MeshAdapt.hpp and .cpp to provide a path to using c++11 member var init

// options start...
// ====================================================================================================================================
CPP11MVITYPE(std::string) options_description_desc CPP11MVI( "adapt options");

// NOTE: Options --directory --output-log --runtest are handled/defined in RunEnvironment
CPP11MVITYPE(std::string) input_mesh CPP11MVI("");
CPP11MVITYPE(std::string) input_geometry CPP11MVI("");
#if IN_MESH_ADAPT
CPP11MVITYPE(GeometryType) input_geometry_type CPP11MVIEQ(= GEOM_NONE);
#endif
CPP11MVITYPE(std::string) fit_geometry_file CPP11MVI("");
CPP11MVITYPE(double) fit_angle CPP11MVIEQ (= 115.0);
CPP11MVITYPE(std::string) fit_3d_file CPP11MVI( "");
CPP11MVITYPE(std::string) output_mesh CPP11MVI("");
CPP11MVITYPE(std::string) block_name_inc CPP11MVI( "");
CPP11MVITYPE(std::string) block_name_exc CPP11MVI( "");
CPP11MVITYPE(int) use_transition_elements CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) previous_adapted_mesh CPP11MVI("");
CPP11MVITYPE(std::string) next_adapted_mesh CPP11MVI("");

// for Salinas
#if defined(STK_BUILT_IN_SIERRA)
CPP11MVITYPE(std::string) rbar_blocks CPP11MVI( "");
#endif
// for Salinas and other codes
//std::string ignore_blocks = "";
// just use block_name_inc to exclude....

#if defined(STK_BUILT_IN_SIERRA)
CPP11MVITYPE(std::string) version_prefix CPP11MVI( "Sierra_");
#else
CPP11MVITYPE(std::string) version_prefix CPP11MVI( "NonSierra_");
#endif

CPP11MVITYPE(std::string) version CPP11MVI( "1.0");
CPP11MVITYPE(int) print_version CPP11MVIEQ(= 0);

//CPP11MVITYPE(int) convert_geometry_parts_OpenNURBS CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) convert_geometry_parts_OpenNURBS CPP11MVI("");

CPP11MVITYPE(std::string) convert CPP11MVI("");
CPP11MVITYPE(std::string) refine CPP11MVI("");
//std::string) refine="";
CPP11MVITYPE(std::string) enrich CPP11MVI("");
CPP11MVITYPE(bool) doRefineMesh CPP11MVIEQ(= true);
CPP11MVITYPE(int) load_balance CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) convert_Hex8_Tet4_24 CPP11MVI( "Hex8_Tet4_24");
CPP11MVITYPE(int) print_info CPP11MVIEQ(=0);
CPP11MVITYPE(int) dihedral_angle_check CPP11MVIEQ(=0);
CPP11MVITYPE(int) skin_mesh CPP11MVIEQ(=0);
CPP11MVITYPE(int) serialized_io_group_size CPP11MVIEQ(= 0);
CPP11MVITYPE(int) remove_original_elements CPP11MVIEQ(= 1);
CPP11MVITYPE(int) output_active_elements_only CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) verify_meshes CPP11MVI("0");
CPP11MVITYPE(int) number_refines CPP11MVIEQ(= 1);
CPP11MVITYPE(int) proc_rank_field CPP11MVIEQ(= 0);
CPP11MVITYPE(int) query_only CPP11MVIEQ(= 0);
CPP11MVITYPE(int) progress_meter CPP11MVIEQ(= 0);
CPP11MVITYPE(int) print_timers CPP11MVIEQ(= 0);
CPP11MVITYPE(int) smooth_geometry CPP11MVIEQ(= 0);
CPP11MVITYPE(int) smoother_niter CPP11MVIEQ(= 1000);
CPP11MVITYPE(double) smoother_tol CPP11MVIEQ(= 1.e-4);
CPP11MVITYPE(int) save_intermediate_meshes CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) adapt_input CPP11MVI( "");
CPP11MVITYPE(int) smooth_use_reference_mesh CPP11MVIEQ(= 1);
CPP11MVITYPE(int) fix_all_block_boundaries CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) ioss_write_options CPP11MVI( "");
CPP11MVITYPE(std::string) ioss_read_options CPP11MVI( "");
CPP11MVITYPE(int) generated_mesh CPP11MVIEQ(= 0);
CPP11MVITYPE(int) snap_geometry CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) internal_test CPP11MVI( "");
#if !defined(NO_GEOM_SUPPORT)
CPP11MVITYPE(int) respect_spacing CPP11MVIEQ(= 1);
#endif
CPP11MVITYPE(int) smooth_surfaces CPP11MVIEQ(= 0);
//double min_spacing_factor CPP11MVIEQ(= 0.25); // range [0,0.5]
CPP11MVITYPE(int) remove_geometry_blocks CPP11MVIEQ(= 0);
CPP11MVITYPE(int) dump_geometry_file CPP11MVIEQ(= 0);
CPP11MVITYPE(int) sync_io_regions CPP11MVIEQ(= 1);
CPP11MVITYPE(int) delete_parents CPP11MVIEQ(= 1);
CPP11MVITYPE(int) print_memory_usage CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) property_map CPP11MVIEQ(= "");
// a list of comma-separated names like Entity, Relation, Field, etc.
CPP11MVITYPE(std::string) memory_multipliers_file CPP11MVI("");
CPP11MVITYPE(int) estimate_memory_usage CPP11MVIEQ(=0);
CPP11MVITYPE(int) precheck_memory_usage CPP11MVIEQ(=0);
CPP11MVITYPE(int) streaming_size CPP11MVIEQ(=0);
CPP11MVITYPE(int) streaming_rank CPP11MVIEQ(=0);
CPP11MVITYPE(int) streaming_pass_start CPP11MVIEQ(= -2);  // FIXME - change to not start from -1 below
CPP11MVITYPE(int) streaming_pass_end CPP11MVIEQ(= -2);
//std::string streaming_instruction CPP11MVIEQ(="");
CPP11MVITYPE(int) streaming_W CPP11MVIEQ(= 0);
CPP11MVITYPE(int) streaming_iW CPP11MVIEQ(= 0);
CPP11MVITYPE(std::string) compute_hmesh CPP11MVI( "");
CPP11MVITYPE(int) print_hmesh_surface_normal CPP11MVIEQ(= 0);
CPP11MVITYPE(int) save_internal_fields CPP11MVIEQ(= 0);

CPP11MVITYPE(int) use_side_map CPP11MVIEQ(= 0); // experimental

CPP11MVITYPE(double) hmesh_factor CPP11MVIEQ(= 0.0);
#define CPP11_COMMA ,
CPP11MV1_ARRAY_INIT(double, hmesh_min_max_ave_factor, 3, tmpvals, {0 CPP11_COMMA 0 CPP11_COMMA 0} );
CPP11MVITYPE(std::string) histograms_root CPP11MVI("cout");
//std::string histogram_options = " CPP11MVI(mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume] }");
CPP11MVITYPE(std::string) histogram_options CPP11MVI( "");

//  Hex8_Tet4_24 (default), Quad4_Quad4_4, Qu
CPP11MVITYPE(std::string) block_name_desc CPP11MVI(
    "block name(s) to convert: there are several options\n"
      "  (1) empty string or option not specified: convert all blocks in the input mesh file\n"
      "  (2) file:my_filename.my_ext (e.g. file:filelist.dat) which will read input block names\n"
      "            from the given file\n"
      "  (3) [+]block_name_1,[+]block_name_2, etc ,block_name_n to include only these blocks, plus sign is optional\n"
      "  (4) a single input block name (e.g. block_3) to be converted\n"
      "  (5) -block_3,-block_5 to exclude blocks from those included (all blocks or include-only blocks), minus sign is mandatory\n"
      "  (6) block_1..block_10 include the range of blocks #1 to #10\n"
      "  (7) any combination of [+] and - options and range (..) option can be specified\n"
      "  (8) Note: you can add the optional specification :Nx or :NX of the number of times to refine\n"
      "        a particular block to any block name specification, e.g. --blocks=1..3:2x,5,6:3x would refine\n"
      "        blocks 1,2 and 3 twice, block 5 'number_refines' times, block 6 three times.\n"
      "Note: wherever you specify block_# this can be replaced with just the #, e.g. \"1,2,4,5\""
      );

#if IN_MESH_ADAPT
CPP11MVITYPE(std::string) convert_options  CPP11MVI(UniformRefinerPatternBase::s_convert_options);
CPP11MVITYPE(std::string) refine_options  CPP11MVI(UniformRefinerPatternBase::s_refine_options);
CPP11MVITYPE(std::string) enrich_options  CPP11MVI(UniformRefinerPatternBase::s_enrich_options);
#endif

CPP11MVITYPE(int) test_memory_elements CPP11MVIEQ(= 0);
CPP11MVITYPE(int) test_memory_nodes CPP11MVIEQ(= 0);

//convert_options = "DEFAULT or one of "+convert_options;
//refine_options = "DEFAULT or one of "+refine_options;
//enrich_options = "DEFAULT or one of "+enrich_options;

// : if not specified, use input mesh name appended with _{converted,refined_#refines,enriched}");

CPP11MVITYPE(std::string) block_name_desc_inc  CPP11MVI( "which blocks to include, specified as: "+block_name_desc);
CPP11MVITYPE(std::string) block_name_desc_exc  CPP11MVI( "which blocks to exclude, specified as: "+block_name_desc);

CPP11MVITYPE(int) help = 0;
// ====================================================================================================================================
// options end
#ifndef NDEBUG
CPP11MVITYPE(int) debug CPP11MVIEQ(= 1);
#else
CPP11MVITYPE(int) debug CPP11MVIEQ(= 0);
#endif

CPP11MVITYPE(int) s_spatialDim CPP11MVIEQ(=0);

CPP11MVITYPE(double) t0   CPP11MVIEQ(= 0.0);
CPP11MVITYPE(double) t1   CPP11MVIEQ(= 0.0);
CPP11MVITYPE(double) cpu0 CPP11MVIEQ(= 0.0);
CPP11MVITYPE(double) cpu1 CPP11MVIEQ(= 0.0);

CPP11MVITYPE(int) i_pass CPP11MVIEQ(= 0);
CPP11MVITYPE(int) m_M CPP11MVIEQ(= 1);
CPP11MVITYPE(int) m_W CPP11MVIEQ(= 1);
CPP11MVITYPE(int) m_iW CPP11MVIEQ(= 0);
CPP11MVITYPE(int) m_M_0 CPP11MVIEQ(= 0);
CPP11MVITYPE(int) m_M_1 CPP11MVIEQ(= 0);
CPP11MVITYPE(int) m_iM CPP11MVIEQ(= 0);

CPP11MVITYPE(std::string) histogram_basic_options  CPP11MVI( "{file_root: "+histograms_root + ", mesh: [edge_length, quality_edge, quality_vol_edge_ratio, volume] }" );
CPP11MVITYPE(int) result CPP11MVIEQ(= 0);
CPP11MVITYPE(unsigned) failed_proc_rank CPP11MVIEQ(= 0u);

CPP11MVITYPE(unsigned) p_rank CPP11MVIEQ(= 0);
CPP11MVITYPE(unsigned) p_size CPP11MVIEQ(= 0);

CPP11MVITYPE(int) first_extra_option_index CPP11MVIEQ(= -1);

CPP11MVITYPE(double) mMeshInputTime       CPP11MVIEQ(= 0.0);
CPP11MVITYPE(double) mMeshOutputTime      CPP11MVIEQ(= 0.0);
CPP11MVITYPE(double) mAdaptTimeOverall    CPP11MVIEQ(= 0.0);
CPP11MVITYPE(double) mAdaptCPUTimeOverall CPP11MVIEQ(= 0.0);
