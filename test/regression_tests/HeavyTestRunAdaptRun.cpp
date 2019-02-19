// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/Percept.hpp>

#include <gtest/gtest.h>
#include <percept/PerceptMesh.hpp>
#include <percept/util/Loops.hpp>

#include <unit_tests/UnitTestSupport.hpp>
#include <percept/fixtures/QuadFixture.hpp>

#include <adapt/main/RunAdaptRun.hpp>
#include <adapt/ElementRefinePredicate.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/RefinerUtil.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <regression_tests/RegressionTestLocalRefiner.hpp>


namespace percept {
  using namespace regression_tests;

  namespace heavy_tests {

#include <regression_tests/RegressionTestFileLoc.hpp>

    static std::string post_fix(int i) {
      return "np" + toString(i);
    }

    //naluFlowPastSquare.e-before-adapt.e-s0002.8.0
    TEST(heavy_rar, run_adapt_run_nalu)
    {
      std::string type = "tri";
      //std::string type = "quad";

      EXCEPTWATCH;

      stk::ParallelMachine pm = MPI_COMM_WORLD ;
      bool debug = false;

      const unsigned p_rank = stk::parallel_machine_rank( pm );
      (void)p_rank;
      const unsigned p_size = stk::parallel_machine_size( pm );
      if (p_size == 1 || p_size == 3 || p_size == 8)
        {
          percept::PerceptMesh eMesh(2u);
          eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);

          eMesh.open("naluFlowPastSquare.e-before-adapt.e-s0002");
          //eMesh.open("before-ref.1.e");

          eMesh.register_and_set_refine_fields();
          eMesh.add_registered_refine_fields_as_input_fields();
          ErrorFieldType * error_field = &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "error_indicator");
          stk::mesh::put_field_on_mesh( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1, nullptr);
          stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);
          eMesh.add_input_field(error_field);
          UniformRefinerPatternBase *localBreakPattern = 0;
          if (type == "quad")
            localBreakPattern = new Local_Quad4_Quad4_N_Transition(eMesh);
          else
            localBreakPattern = new Local_Tri3_Tri3_N_HangingNode(eMesh);

          eMesh.commit();
          eMesh.output_active_children_only(false);

          stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

          ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);
          TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, *localBreakPattern, 0);

          breaker.setRemoveOldElements(false);
          breaker.setAlwaysInitializeNodeRegistry(false);

          RefinerUtil::rebuild_family_tree(eMesh, debug);
          AdaptedMeshVerifier::check_parent_element_field(eMesh, "after rebuild", debug);

          breaker.initializeRefine();
          erp.setMarkNone(true);
          breaker.initializeDB(true);
          erp.setMarkNone(false);

          //SetOfEntities emptySet(*eMesh.get_bulk_data());
          //breaker.replaceNodeRegistryOwnership(emptySet, eMesh.element_rank());

          AdaptedMeshVerifier adaptedMeshVerifier(debug);
          if (!adaptedMeshVerifier.isValid(eMesh, true))
            throw std::runtime_error("Invalid initial mesh");

          std::string input = "{family_tree_extension: ft, error_indicator_field: error_indicator, "
            "marker: {type: physical, physical_error_criterion: 1.0}, "
            "max_number_elements_fraction: 1.2, max_refinement_level: 2, extra_output: yes"
            // following are optional and have defaults
            ", refine_field: refine_field, refine_level_field: refine_level,"
            "transition_element_field: transition_element, parent_element_field: parent_element,"
            "new_nodes_field: new_nodes"
            "}";

          std::string output = "rar.yaml";

          percept::RunAdaptRunInfo rar(eMesh, input, output, debug);
          rar.create();
          rar.create_marker();

          rar.m_marker->mark();

          if (1)
            {
              RefinerUtil::save_node_registry(eMesh, breaker.getNodeRegistry(), "HeavyTestLocalRefiner");
              PerceptMesh& eMeshNR = eMesh;
              NodeRegistry nr1(eMeshNR);
              RefinerUtil::rebuild_node_registry(eMeshNR, nr1, true, &eMesh, &breaker.getNodeRegistry(), true);
            }

          //eMesh.save_as("before-ref.1.e");

          breaker.refine();
          std::cout << "tet_transition number elements = "
                    << eMesh.get_number_elements() << std::endl;

          eMesh.save_as("naluFlowPastSquare_ft.e");
          eMesh.output_active_children_only(true);
          eMesh.save_as("naluFlowPastSquare.e");

          delete localBreakPattern;
        }
    }

    template<class LocalBreakPattern>
    static void do_het_local_corner_refine_sidesets(PerceptMesh& eMesh, int num_time_steps, std::vector<double>& vbbox, const std::string& file,
                                                    bool save_intermediate=false, bool delete_parents=false, int num_ref_passes = 6,  int num_unref_passes = 4,
                                                    int num_loops = 1)
    {
      EXCEPTWATCH;
      stk::ParallelMachine pm = MPI_COMM_WORLD ;

      const unsigned p_size = stk::parallel_machine_size( pm );
      if (p_size==1 || p_size == 3 || p_size >= 8)
        {
          LocalBreakPattern break_het_to_het_N(eMesh);
          int scalarDimension = 0; // a scalar
          stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
          RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

          // for plotting, use doubles, for internal use, use int
          RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
          stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
          stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

          eMesh.register_and_set_refine_fields();
          //eMesh.set_remove_io_orig_topo_type(true);
          //eMesh.set_sync_io_regions(true);

          eMesh.commit();
          //eMesh.print_info("input", 2);

          std::cout << "het_local initial number elements= " << eMesh.get_number_elements() << std::endl;

          eMesh.save_as( output_files_loc+"het_local_"+percept::heavy_tests::post_fix(p_size)+".e.0");

          stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

          ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
          TransitionElementAdapter<ElementRefinePredicate> *breaker_p;

          std::vector<std::string> do_wedge_bl;
          {
            const stk::mesh::PartVector & mesh_parts = eMesh.get_fem_meta_data()->get_mesh_parts();
            for (unsigned ipart=0; ipart < mesh_parts.size(); ipart++) {
              stk::mesh::Part* part = mesh_parts[ipart];

              if (part->primary_entity_rank() != stk::topology::ELEMENT_RANK
                  || stk::mesh::is_auto_declared_part(*part)) continue;

              stk::topology stk_topo = stk::mesh::get_topology(eMesh.get_fem_meta_data()->get_cell_topology(*part));

              if (stk_topo == stk::topology::WEDGE_6) {
                do_wedge_bl.push_back(part->name());
              }
            }
          }

          stk::mesh::Selector wedge_selector;
            {
              stk::mesh::PartVector pv;
              for (unsigned ii = 0; ii < do_wedge_bl.size(); ++ii)
                {
                  stk::mesh::Part* part = eMesh.get_fem_meta_data()->get_part(do_wedge_bl[ii]);
                  VERIFY_OP_ON(part, !=, 0, "can't find wedge part "+do_wedge_bl[ii]);
                  VERIFY_OP_ON(part->topology(), ==, stk::topology::WEDGE_6, "block has wrong topology: "+do_wedge_bl[ii]);
                  if (eMesh.get_rank() == 0) std::cout << "TEA_SpecialWedgeRefinement:: found part = " << part->name() << std::endl;
                  pv.push_back(part);
                }
              wedge_selector = stk::mesh::selectUnion(pv);
              //breaker_p = new TEA_SpecialWedgeRefinement<ElementRefinePredicate>(erp, eMesh, break_het_to_het_N, proc_rank_field, &wedge_selector);
              breaker_p = new TEA_SpecialWedgeRefinement<ElementRefinePredicate>(erp, eMesh, break_het_to_het_N, proc_rank_field, &wedge_selector, true, false);
            }
          TransitionElementAdapter<ElementRefinePredicate>& breaker = *breaker_p;
          if (eMesh.get_rank() == 0)
            std::cout << "type= " << eMesh.demangle(typeid(breaker).name()) << std::endl;

          // special for local adaptivity
          breaker.setRemoveOldElements(false);
          breaker.setAlwaysInitializeNodeRegistry(false);

          // FIXME
          //breaker.setAddChildrenToParts(false);

          SetElementFieldBoundingBox set_ref_field(eMesh, vbbox, refine_field);

          int iplot=0;
          if (1)
            {
              char buf[1000];
              sprintf(buf, "%04d", iplot);
              if (iplot == 0)
                eMesh.save_as(file+"_anim.e", double(iplot));
              else
                eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
              ++iplot;
            }
          for (int loop=0; loop < num_loops; ++loop)
            {
              for (int ipass=0; ipass < num_ref_passes; ipass++)
                {

                  int nel = eMesh.get_number_elements();
                  if (eMesh.get_rank() == 0)
                    std::cout << "P[" << eMesh.get_rank() << "] before refine ipass= " << ipass << " local number elements= "
                              << nel << std::endl;

                  elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                  if (1)
                    {
                      char buf[1000];
                      sprintf(buf, "%04d", ipass);
                      eMesh.save_as(file+"_set_field.e-s"+std::string(buf), double(iplot));
                      ++iplot;
                    }

                  breaker.refine();

                  //eMesh.print_info("refined", 2);

                  nel = eMesh.get_number_elements();
                  if (eMesh.get_rank() == 0)
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " local number elements= "
                              << nel << std::endl;
                  //eMesh.dump_elements_compact();
                  if (1)
                    {
                      char buf[1000];
                      sprintf(buf, "%04d", iplot);
                      eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                      ++iplot;
                    }

                }

              for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                {
                  elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                  breaker.unrefine();
                  int nel = eMesh.get_number_elements();
                  if (eMesh.get_rank() == 0)
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " local number elements= " << nel << std::endl;
                  if (1)
                    {
                      char buf[1000];
                      sprintf(buf, "%04d", iplot);
                      eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                      ++iplot;
                    }
                }
            }

          SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1);

          eMesh.save_as(output_files_loc+file+"_tmp_square_sidesets_local_unref_"+percept::heavy_tests::post_fix(p_size)+".e");

          eMesh.save_as(output_files_loc+file+"_sidesets_final_local_"+percept::heavy_tests::post_fix(p_size)+".e.s-"+toString(num_time_steps) );

          for (int iunref=0; iunref < (num_unref_passes? 10 : 0); iunref++)
            {
              elementOpLoop(*eMesh.get_bulk_data(), set_ref_field_val_unref_all, refine_field);
              if (eMesh.get_rank() == 0)
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
              breaker.unrefine();

              if (1)
                {
                  char buf[1000];
                  sprintf(buf, "%04d", iplot);
                  eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                  ++iplot;
                }

              int nel = eMesh.get_number_elements();
              if (eMesh.get_rank() == 0)
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " local number elements= " << nel << std::endl;
            }


          if (delete_parents)
            breaker.deleteParentElements();
          int nel = eMesh.get_number_elements();
          if (eMesh.get_rank() == 0)
            std::cout << "local final number elements= " << nel << std::endl;
          eMesh.save_as(output_files_loc+file+"_sidesets_final_unrefed_local_"+percept::heavy_tests::post_fix(p_size)+".e."+toString(num_time_steps) );

          delete breaker_p;

          // end_demo
        }
    }

  } // namespace unit_tests
} // namespace percept


