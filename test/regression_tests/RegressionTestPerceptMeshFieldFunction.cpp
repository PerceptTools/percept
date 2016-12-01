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
#include <string>
#include <typeinfo>

#include <math.h>

#include <stk_mesh/base/FieldParallel.hpp>

#include <gtest/gtest.h>

#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept
  {
    namespace regression_tests
    {

#define EXTRA_PRINT 0

      /// configuration: you can choose where to put the generated Exodus files (see variables input_files_loc, output_files_loc)
      /// The following defines where to put the input and output files created by this set of functions

#if 1
      const std::string path_sep = "._.";
      const std::string input_files_loc="./input_files"+path_sep;
      const std::string output_files_loc="./output_files"+path_sep;
#else
      const std::string input_files_loc="./input_files/";
      const std::string output_files_loc="./output_files/";
#endif

      static double pressure_value = 123.4;

#if 1
      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      TEST(perceptMesh, open_new_close_PerceptMesh)
      {
        EXCEPTWATCH;

        // start_demo_open_new_close_PerceptMesh
        PerceptMesh eMesh(3u);
        eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
        int scalarDimension = 0; // a scalar
        int vectorDimension = 3;

	stk::mesh::FieldBase* pressure_field = eMesh.add_field("pressure", stk::topology::NODE_RANK, scalarDimension);
        eMesh.add_field("velocity", stk::topology::NODE_RANK, vectorDimension);
        eMesh.add_field("element_volume", stk::topology::ELEMENT_RANK, scalarDimension);

        eMesh.commit();

        // create a field function from the new pressure field
        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);

        // set the value of the pressure field to a constant everywhere
        ConstantFunction initPressureValues(pressure_value, "initPVal");
        ff_pressure.interpolateFrom(initPressureValues);

        if (eMesh.get_rank()== 0) eMesh.print_fields("Pressure");
        //exit(1);

        // here we could evaluate this field function
        double x=0.123, y=0.234, z=0.345, time=0.0;
        std::cout << "P[" << eMesh.get_rank() << "] "
                  << "before write ff_pressure = " << eval(x,y,z,time, ff_pressure) << std::endl;

        //eval_print(x, y, z, time, ff_pressure);

        double pval = eval(x, y, z, time, ff_pressure);
        EXPECT_DOUBLE_EQ(pval, pressure_value);

        eMesh.save_as(output_files_loc+"cube_with_pressure.e");
        eMesh.close();

        // end_demo

        // start_demo_open_new_close_PerceptMesh_1
        // open the file we previously saved with the new fields
        eMesh.open_read_only(output_files_loc+"cube_with_pressure.e");

        // get the pressure field
        pressure_field = eMesh.get_field(stk::topology::NODE_RANK, "pressure");

        // FIXME
        std::vector< const stk::mesh::FieldBase * > sync_fields( 1 , pressure_field );
        stk::mesh::communicate_field_data( eMesh.get_bulk_data()->aura_ghosting() , sync_fields );
        // FIXME

        //if (1 || eMesh.get_rank()== 0) eMesh.print_fields("Pressure");

        FieldFunction ff_pressure_1("ff_pressure", pressure_field, eMesh, 3, 1);
        ff_pressure_1.add_alias("P");
        StringFunction sf_pressure("P");
        std::cout << "P[" << eMesh.get_rank() << "] "
                  << "after read ff_pressure = " << eval(x,y,z,time, ff_pressure_1) << std::endl;

        // a point-source at the origin
#define EXACT_SOL log(sqrt( x*x + y*y + z*z) + 1.e-10)

        StringFunction sf_exact_solution(EXPAND_AND_QUOTE(EXACT_SOL), Name("sf_exact_solution"), 3, 1);
        StringFunction sf_error = sf_exact_solution - sf_pressure;

        std::cout << "P[" << eMesh.get_rank() << "] "
                  << "sf_pressure = " << eval(x,y,z,time, sf_pressure) << std::endl;
        //!eval_print(x,y,z,time, sf_error);
        std::cout << "P[" << eMesh.get_rank() << "] "
                  << "sf_error = " << eval(x,y,z,time, sf_error) << std::endl;
        double val_cpp = EXACT_SOL - pressure_value;
        double val_sf  = eval(x,y,z,time, sf_error);
        EXPECT_DOUBLE_EQ(val_sf, val_cpp);
        // end_demo

      }

#endif

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      TEST(perceptMesh, open_new_close_PerceptMesh_2)
      {
        EXCEPTWATCH;

        double x=0.123, y=0.234, z=0.345, time=0.0;

        // start_demo_open_new_close_PerceptMesh_2
        PerceptMesh eMesh(3u);
        // open the file we previously saved with the new fields
        eMesh.open_read_only(input_files_loc+"cube_with_pressure.e");

        eMesh.print_info("Info after reading mesh");

        //eMesh.print_fields();
        stk::mesh::FieldBase *f_coords = eMesh.get_field(stk::topology::NODE_RANK, "coordinates");

        // create a field function from the existing coordinates field
        FieldFunction ff_coords("ff_coords", f_coords, eMesh, 3, 3);

        // here we could evaluate this field function
        eval_vec3_print(x, y, z, time, ff_coords);

        // get the pressure field
        stk::mesh::FieldBase* pressure_field = eMesh.get_field(stk::topology::NODE_RANK, "pressure");

        // FIXME
        std::vector< const stk::mesh::FieldBase * > sync_fields( 1 , pressure_field );
        stk::mesh::communicate_field_data( eMesh.get_bulk_data()->aura_ghosting() , sync_fields );
        // FIXME

        //double * pdata = eMesh.node_field_data(pressure_field, 1);

        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);
        ff_pressure.add_alias("P");
        StringFunction sf_pressure("P");

        // a point-source at the origin
#define EXACT_SOL log(sqrt( x*x + y*y + z*z) + 1.e-10)

        StringFunction sf_exact_solution(EXPAND_AND_QUOTE(EXACT_SOL), Name("sf_exact_solution"), 3, 1);
        StringFunction sf_error = sf_exact_solution - sf_pressure;

        eval_print(x,y,z,time, sf_error);
        double val_cpp = EXACT_SOL - pressure_value;
        double val_sf  = eval(x,y,z,time, sf_error);
        EXPECT_DOUBLE_EQ(val_sf, val_cpp);
        // end_demo

      }

#define EXPECT_CATCH(expression, name)                                  \
      {                                                                 \
        bool didCatch_ ## name = false;                                 \
        try {                                                           \
          expression ;                                                  \
        }                                                               \
        catch ( const std::exception & X ) {                            \
          std::cout << # name << "  expected to catch this exception: " << X.what() << std::endl; \
          didCatch_ ## name = true;                                     \
        }                                                               \
        catch( ... ) {                                                  \
          std::cout << " Caught unknown exception"                      \
                    << std::endl ;                                      \
          std::cout.flush();                                            \
          didCatch_ ## name = false;                                    \
        }                                                               \
        EXPECT_TRUE(didCatch_ ## name);                                 \
      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      TEST(perceptMesh, open_new_close_PerceptMesh_3)
      {
        EXCEPTWATCH;

        // start_demo_open_new_close_PerceptMesh_3
        PerceptMesh eMesh(3u);
        eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
        int scalarDimension = 0; // a scalar
        int vectorDimension = 3;

        stk::mesh::FieldBase* pressure_field = eMesh.add_field("pressure", stk::topology::NODE_RANK, scalarDimension);
        eMesh.add_field("velocity", stk::topology::NODE_RANK, vectorDimension);
        eMesh.add_field("element_volume", stk::topology::ELEMENT_RANK, scalarDimension);

        eMesh.commit();

        EXPECT_CATCH( eMesh.commit() , commit_again);

        // create a field function from the new pressure field
        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);

        // set the value of the pressure field to a constant everywhere
        ConstantFunction initPressureValues(pressure_value, "initPVal");
        ff_pressure.interpolateFrom(initPressureValues);

        // save
        eMesh.save_as(output_files_loc+"cube_with_pressure_3.e");
        eMesh.close();

        EXPECT_CATCH( eMesh.print_info("bad", 1) , mesh_closed_try_print);

        // end_demo

      }


      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      TEST(perceptMesh, open_new_reopen_PerceptMesh)
      {
        EXCEPTWATCH;

        // start_demo_open_new_reopen_PerceptMesh
        PerceptMesh eMesh(3u);
        eMesh.new_mesh(GMeshSpec("3x3x3|bbox:0,0,0,1,1,1"));  // create a 3x3x3 hex mesh in the unit cube
        int scalarDimension = 0; // a scalar
        int vectorDimension = 3;

        eMesh.add_field("pressure", stk::topology::NODE_RANK, scalarDimension);
        eMesh.add_field("velocity", stk::topology::NODE_RANK, vectorDimension);
        eMesh.add_field("element_volume", stk::topology::ELEMENT_RANK, scalarDimension);

        eMesh.commit();

        /// reopen the mesh to allow for more fields to be added - note that this involves a db write/read operation
        eMesh.reopen(output_files_loc+"optional_temp_filename.e");
        stk::mesh::FieldBase* momentum_field = eMesh.add_field("momentum", stk::topology::NODE_RANK, vectorDimension);
        eMesh.commit();

        // create a field function from the new pressure field
        stk::mesh::FieldBase *pressure_field = eMesh.get_field(stk::topology::NODE_RANK, "pressure");
        FieldFunction ff_pressure("ff_pressure", pressure_field, eMesh, 3, 1);

        // set the value of the pressure field to a constant everywhere
        ConstantFunction initPressureValues(pressure_value, "initPVal");
        ff_pressure.interpolateFrom(initPressureValues);

        // set the momentum field
        std::vector<double> momentum_values(3);
        momentum_values[0] = 2034.5;
        momentum_values[1] = 2134.5;
        momentum_values[2] = 2234.5;
        ConstantFunctionVec initMomentumValues(momentum_values, "initMomVal");

        // create a field function from the new momentum field
        FieldFunction ff_momentum("ff_momentum", momentum_field, eMesh, 3, 3);
        ff_momentum.interpolateFrom(initMomentumValues);

        // save
        eMesh.save_as(output_files_loc+"cube_with_pressure_and_momentum.e");
        eMesh.close();


        // end_demo

      }

      //======================================================================================================================
      //======================================================================================================================
      //======================================================================================================================
      TEST(perceptMesh, inner_skin)
      {
        EXCEPTWATCH;

        //const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD);
        const unsigned p_size = stk::parallel_machine_size( MPI_COMM_WORLD);

        if (p_size == 1 || p_size == 3)
          {
            // start_demo_open_new_reopen_PerceptMesh
            PerceptMesh eMesh(3u);
            eMesh.open("input_files._.hex-7block_R2.g");
            bool make_part_io_part = true;
            eMesh.add_part("inner_skin_part", make_part_io_part);

            eMesh.commit();

            eMesh.initializeIdServer();
            double coords[3] = { 100, 200, 300};
            std::cout << " get_number_edges= " << eMesh.get_number_edges()
                      << " get_number_elements_locally_owned= " << eMesh.get_number_elements_locally_owned()
                      << std::endl;

            eMesh.get_skin_part("inner_skin_part");

            // save
            eMesh.save_as(output_files_loc+"hex-7block-skinned.g");

            if (p_size == 1)
              {
                stk::mesh::Entity node = eMesh.get_node(coords[0], coords[1], coords[2]);
                stk::mesh::Entity node1 = eMesh.get_bulk_data()->get_entity(eMesh.node_rank(), 1ull);
                EXPECT_TRUE(eMesh.is_valid(node));
                node = eMesh.get_closest_node(coords[0], coords[1], coords[2], 0, 0);
                EXPECT_TRUE(eMesh.is_valid(node));


                stk::mesh::EntityId id = eMesh.getNextId(eMesh.node_rank());
                eMesh.get_bulk_data()->modification_begin();
                eMesh.createOrGetNode(id, coords);
                stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
                eMesh.get_bulk_data()->modification_end();


                eMesh.dump_elements();
                eMesh.dump_elements_compact();
                eMesh.dump_vtk(node1, "junk.vtk");
                eMesh.dump_vtk("junk1.vtk", true);
                eMesh.close();

                eMesh.dump(output_files_loc+"hex-7block-skinned.g");

                PerceptMesh eMesh1(3u);
                eMesh1.open_read_only("input_files._.hex-7block_R2.g");
                eMesh1.print_all();
                eMesh1.delete_side_sets();
                eMesh1.renumber_sides_for_exodus();
              }

            // end_demo
          }
      }


    }//    namespace unit_tests
  }//  namespace percept

