// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <gtest/gtest.h>

#include <percept/function/FieldFunction.hpp>
#include <percept/function/FunctionOperator.hpp>
#include <percept/function/ConstantFunction.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/RunEnvironment.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <adapt/Colorer.hpp>
#include <percept/ExceptionWatch.hpp>


  namespace percept
  {
    namespace unit_tests
    {

#include "RegressionTestFileLoc.hpp"

      static stk::diag::Writer &
      dw()
      {
        //static stk::diag::Writer s_diagWriter(dwout().rdbuf(), 0);
        int dw_enabled = 1;
        static stk::diag::Writer s_diagWriter(std::cout.rdbuf(), dw_enabled);

        s_diagWriter.setPrintMask(percept::LOG_NORM | percept::LOG_ALWAYS);

        return s_diagWriter;
      }

#define EXTRA_PRINT 0

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      TEST(mesh_colorer, test1)
      {
        EXCEPTWATCH;

        dw().m(percept::LOG_MESH_COLORER) << "mesh_colorer::test1 " << stk::diag::dendl;

        percept::PerceptMesh eMesh(3u);
        if (eMesh.get_parallel_size() <= 3)
          {
            const size_t numxyz=3;
            const size_t num_x = numxyz;
            const size_t num_y = numxyz;
            const size_t num_z = numxyz;
            std::string config_mesh = 
              Ioss::Utils::to_string(num_x) + "x" + 
              Ioss::Utils::to_string(num_y) + "x" +
              Ioss::Utils::to_string(num_z) + "|bbox:0,0,0,1,1,1";
	
            eMesh.new_mesh(percept::GMeshSpec(config_mesh));
            int vectorDimension = 0;
            stk::mesh::FieldBase *element_color_field = eMesh.add_field("element_colors", stk::topology::ELEMENT_RANK, vectorDimension);
            eMesh.commit();

            std::vector<stk::mesh::EntityRank> mer;  mer.push_back(stk::topology::ELEMENT_RANK);
            Colorer meshColorer(mer);
            unsigned elementType = 0u;
            meshColorer.color(eMesh, &elementType, 0, element_color_field);
            eMesh.save_as(output_files_loc+"cube_colored.e");
            //std::cout << "Mesh coloring info: " << meshColorer.getElementColors() << std::endl;
          }
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================

      TEST(mesh_colorer, test_quad)
      {
        EXCEPTWATCH;

        dw().m(percept::LOG_MESH_COLORER) << "mesh_colorer::test_quad " << stk::diag::dendl;

        percept::PerceptMesh eMesh(2u);
        if (eMesh.get_parallel_size() == 1 || eMesh.get_parallel_size() == 3)
          {
            eMesh.open(input_files_loc+"break_test._.quad._.square._.square_quad4.e");
            int vectorDimension = 0;
            stk::mesh::FieldBase *element_color_field = eMesh.add_field("element_colors", stk::topology::ELEMENT_RANK, vectorDimension);
            eMesh.commit();

            std::vector<stk::mesh::EntityRank> mer;  mer.push_back(eMesh.face_rank());
            Colorer meshColorer(mer);
            unsigned elementType = 0u;
            meshColorer.color(eMesh, &elementType, 0, element_color_field);
            //std::cout << "Mesh coloring info: " << meshColorer.getElementColors() << std::endl;
            eMesh.print_info();
            eMesh.dump();
            eMesh.save_as(output_files_loc+"square_quad4_colored.e");
          }
      }

    }
  }

