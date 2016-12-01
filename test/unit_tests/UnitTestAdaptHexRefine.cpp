// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
#include <adapt/PredicateBasedElementAdapter.hpp>
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

namespace percept {
  namespace unit_tests {

    /**  to do:
     * 1. unit test of DiscretizeHex::discretize
     *    1a. all mark patterns
     *    2a. compute volume; expect certain numbers of elements of each type, plot
     * 2. check of random marks for 3x3x3 cube
     * 3. larger, parallel
     */

#if 0
      TEST(regr_localRefiner, hex_transition)
      {
        bool do_test = false;
        //stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          bool do_bootstrap_mesh = true;
          if (do_bootstrap_mesh)
          {
            const unsigned n = 3;
            const unsigned nx = n , ny = n, nz = n;

            percept::PerceptMesh eMesh(3u);
            std::string gmesh_spec = toString(nx)+"x"+toString(ny)+"x"+toString(nz)+"|bbox:-1,-1,-1,1,1,1";
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();

            eMesh.save_as(input_files_loc+"hex_cube_hex8_0.e");
            //eMesh.print_info("test1",2);
            //eMesh.dump_vtk("sqhex3.vtk",false);
          }
          //if (1) return;

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"hex_cube_hex8_0.e");

            eMesh.output_active_children_only(false);

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            do_hex_local_corner_refine_sidesets<Local_Hex8_Hex8_N_Transition>(eMesh, 1, false, true, 2, 0);

            if (0)
              {
                percept::RebalanceMesh rb(eMesh, 0, true);
                double imb_before, imb_after;
                rb.rebalance(imb_before, imb_after);
              }

          }
        }
      }
#endif

  } // namespace unit_tests
} // namespace percept


