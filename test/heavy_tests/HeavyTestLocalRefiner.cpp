// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include <math.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <gtest/gtest.h>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <gtest/gtest.h>
#include <boost/lexical_cast.hpp>
#include <stk_io/IossBridge.hpp>

#include <percept/Percept.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>

#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <percept/RunEnvironment.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

#include <adapt/IEdgeBasedAdapterPredicate.hpp>
#include <adapt/ElementRefinePredicate.hpp>
#include <adapt/PredicateBasedElementAdapter.hpp>
#include <adapt/PredicateBasedEdgeAdapter.hpp>

#include <percept/function/ElementOp.hpp>

#include <regression_tests/RegressionTestLocalRefiner.hpp>

  namespace percept
  {
    namespace heavy_tests
    {
      using namespace regression_tests;

#if 1
      static const std::string path_sep = "._.";
      static const std::string input_files_loc="./input_files"+path_sep;
      static const std::string output_files_loc="./output_files"+path_sep;
#else
      static const std::string input_files_loc="./input_files/";
      static const std::string output_files_loc="./output_files/";
#endif

      //=============================================================================
      //=============================================================================
      //=============================================================================


      static void do_moving_shock_test_large_test(int num_time_steps, bool save_intermediate=false, bool delete_parents=false)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //         shock_width = 1./50.; //1./25.0;
        //         shock_diff_criterion = 0.1;
        shock_width = 1./500.; //1./25.0;
        double shock_diff_criterion = 0.1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 2 || p_size == 8 || p_size == 16)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"OUO_large_testh1tet.g");

            Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.register_and_set_refine_fields();
            eMesh.commit();

            eMesh.output_active_children_only(true);

            int num_elems = eMesh.get_number_elements();
            if (eMesh.get_rank()==0) std::cout << "moving_shock initial number elements= " << num_elems << std::endl;

            std::string fileName = output_files_loc+"OUO_large_scale_moving_shock.e";
            eMesh.save_as(fileName, 0.0);

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            PlaneShock shock;
            shock.plane_point_init[0] = 0.0;
            shock.plane_point_init[1] = -0.07;
            shock.plane_point_init[2] = -0.35;
            shock.plane_normal[0] = 0;
            shock.plane_normal[1] = .25/std::sqrt(.25*.25+1.);
            shock.plane_normal[2] = 1./std::sqrt(.25*.25+1.);

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, &univ_selector, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate>
              breaker(srp,
                      eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            double threshold = 1.0;  // force rebal every time
            breaker.setDoRebalance(true, threshold);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.2*0.08*(10./((double)num_time_steps));
            double shock_displacement = 0.0;
            int num_ref_passes = 2;
            int num_unref_passes = 3;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                if (eMesh.get_rank()==0) std::cout << " istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    if (eMesh.get_rank()==0) std::cout << " ipass= " << ipass <<  std::endl;
                    breaker.doBreak();
                    num_elems = eMesh.get_number_elements();
                    if (eMesh.get_rank()==0) std::cout << " done... ipass= " << ipass << " moving_shock number elements= " << num_elems << std::endl;
                  }

                breaker.getNodeRegistry().init_entity_repo();
                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
                    breaker.buildUnrefineList(elements_to_unref);

                    if (eMesh.get_rank()==0) std::cout << " iunref_pass= " << iunref_pass <<  std::endl;
                    breaker.unrefineTheseElements(elements_to_unref);
                    num_elems = eMesh.get_number_elements();
                    if (eMesh.get_rank()==0) std::cout << " done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << num_elems << std::endl;
                  }

                if (delete_parents && istep == num_time_steps-1)
                  breaker.deleteParentElements();

                std::ostringstream ostr;
                ostr << std::setfill('0') << std::setw(5) << (istep+1);
                fileName = std::string(output_files_loc+"OUO_large_scale_moving_shock.e-s") + ostr.str();

                eMesh.save_as(fileName, (double)(istep+1));

                shock_displacement += delta_shock_displacement;
              }

            for (int iunref=0; iunref < 10; iunref++)
              breaker.unrefineAll();

            breaker.deleteParentElements();

            num_elems = eMesh.get_number_elements();
            if (eMesh.get_rank()==0) std::cout << "moving_shock final number elements= " << num_elems << std::endl;
            eMesh.save_as(output_files_loc+"OUO_large_scale_final_unrefed_moving_shock_"+"np"+toString(p_size)+".e."+toString(num_time_steps) );

            stk::diag::printTimersTable(sierra::Env::outputP0(), breaker.rootTimer(),
                                        stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME,
                                        false, eMesh.get_bulk_data()->parallel());

            // end_demo
          }
      }

      TEST(heavy_localRefiner, break_tet_to_tet_N_5_EdgeBased_moving_shock_large_test)
      {
        const bool do_full_demo = false;
        if (do_full_demo)
          {
            int num_time_steps = 10;  // 10 for stress testing
            for (int istep=1; istep <= num_time_steps; istep++)
              do_moving_shock_test_large_test(istep, false, true);
          }
        else
          // normal regression testing
          {
            int num_time_steps = 5;  // 30 for making animation
            do_moving_shock_test_large_test(num_time_steps, true);
            //int num_time_steps = 3;  // 10 for stress testing
            //do_moving_shock_test_large_test(num_time_steps);
          }
      }

    }
  }

