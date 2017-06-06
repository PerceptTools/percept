// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#include <percept/Percept.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )

#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#if ENABLE_SMOOTHER3
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherNewton.hpp>
#endif
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherAlgebraic.hpp>
#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_total_element_metric.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_update_coordinates.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>

#include <iostream>
#include <cstdlib>
#include <array>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

#define DO_TESTS 1
#if DO_TESTS

#include <percept/PerceptUtils.hpp>
#include <percept/Util.hpp>

namespace std {
}

namespace percept
{

  namespace heavy_tests
  {


#define EXTRA_PRINT 0
      static int s_par_size_max = 2;

#if 1
      static const std::string path_sep = "._.";
      static const std::string input_files_loc="./input_files"+path_sep;
      static const std::string output_files_loc="./output_files"+path_sep;
#else
      static const std::string input_files_loc="./input_files/";
      static const std::string output_files_loc="./output_files/";
#endif

#define EXTRA_PRINT 0

      static std::string procs_string[9] = {"np0", "np1", "np2", "np3", "np4", "np5", "np6", "np7", "np8"};

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run PMMParallelShapeImprover on a domain with one side perturbed
      percept::PerceptMesh *setup_quad_4_tests()
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
        percept::PerceptMesh * eMeshP = 0;

        if (1 || p_size <= par_size_max)
          {
            const unsigned nele = 12;
            //const unsigned nx = nele , ny = nele , nz = p_size*nele ;
            const unsigned nx = nele , ny = nele;

            bool sidesets=true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets);
            fixture.set_bounding_box(0,1,0,1);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            eMeshP = new percept::PerceptMesh(&fixture.meta_data, &fixture.bulk_data);
            percept::PerceptMesh& eMesh = *eMeshP;

            //eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_4_smooth.0.e");
            eMesh.set_save_internal_fields(false);

            eMesh.reopen();
            eMesh.add_coordinate_state_fields();
            eMesh.add_spacing_fields();
            eMesh.register_and_set_refine_fields();
            eMesh.register_and_set_smoothing_fields();
            stk::mesh::FieldBase *proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, 0);
            eMesh.commit();
            eMesh.set_proc_rank_field(proc_rank_field);

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( stk::topology::NODE_RANK );

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_elements_in_bucket = bucket.size();

                  for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity node = bucket[iEntity];

                      double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                      double iy = data[1]; // /double(nele);
                      iy = iy*iy;
                      data[1] = iy; // *double(nele);
                    }
                }
              }
            SpacingFieldUtil sfu(eMesh);
            sfu.compute_spacing_field();
            eMesh.save_as(input_files_loc+"quad_4_smooth.1.e");

            // save state of original mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1"), eMesh.get_coordinates_field());

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                if (boundarySelector_1(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity node = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                        double ix = data[0]; // /double(nele);
                        //double bump_size=2.8; // 0.8
                        double bump_size=2.8; // 0.8
                        data[1] += (ix)*(1.0-ix)*bump_size; //*double(nele);
                        //std::cout << "tmp srk surface 1 node = " << data[0] << " " << data[1] << std::endl;
                      }
                  }
              }

            // save state of projected mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"quad_4_smooth.0_perturbed.e");
          }
        return eMeshP;
      }

#define SMOOTHER  percept::ReferenceMeshSmootherConjugateGradientImpl<STKMesh>
      //#define SMOOTHER  percept::ReferenceMeshSmootherConjugateGradient

      TEST(heavy_perceptMeshSmoother, quad_4)
      {
        PerceptMesh& eMesh = *setup_quad_4_tests();
        eMesh.setProperty("smoother_metric_type", "shape");
        eMesh.setProperty("ReferenceMeshSmoother_do_anim","1");
        stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
        stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
        stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
        stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
        stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

        int  msq_debug             = 0; // 1,2,3 for more debug info
        bool always_smooth         = true;

        SMOOTHER pmmpsi(&eMesh, &boundarySelector, 0, 1001, 1.e-4, 1);
        pmmpsi.run( always_smooth, msq_debug);

        eMesh.save_as(output_files_loc+"quad_4_si_smooth.1.e");
        delete &eMesh;
      }

#if ENABLE_SMOOTHER3
      TEST(heavy_perceptMeshSmoother, quad_4_Newton)
      {
        PerceptMesh& eMesh = *setup_quad_4_tests();

        stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
        stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
        stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
        stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
        stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

        int  msq_debug             = 0; // 1,2,3 for more debug info
        bool always_smooth         = true;
        int do_anim = 0;

        percept::ReferenceMeshSmootherNewton pmmpsi(&eMesh, &boundarySelector, 0, 1001, 1.e-4, 1);
        pmmpsi.m_do_animation = do_anim;
        pmmpsi.run( always_smooth, msq_debug);

        eMesh.save_as(output_files_loc+"quad_4_si_smooth.1.e");
        delete &eMesh;
      }
#endif

      TEST(heavy_perceptMeshSmoother, quad_4_Algebraic)
      {
        PerceptMesh& eMesh = *setup_quad_4_tests();

        stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
        stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
        stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
        stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
        stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

        int  msq_debug             = 0; // 1,2,3 for more debug info
        bool always_smooth         = true;
        int do_anim = 1;

        double drop_off_coeffs[3] = {1,1,1};  // FIXME
        int nlayers_drop_off = 10;
        int niter = 1;
        percept::ReferenceMeshSmootherAlgebraic pmmpsi(&eMesh, &boundarySelector, 0, niter, 1.e-4, drop_off_coeffs, nlayers_drop_off);
        pmmpsi.m_do_animation = do_anim;
        pmmpsi.run( always_smooth, msq_debug);

        eMesh.save_as(output_files_loc+"quad_4_si_smooth.1.e");
        delete &eMesh;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run new parallel smoother on an interior-node perturbed domain
      TEST(DISABLED_heavy_perceptMeshSmoother, quad_4_1)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = 1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
        if (p_size <= par_size_max)
          {
            //const unsigned p_rank = stk::parallel_machine_rank( pm );
            //const unsigned p_size = stk::parallel_machine_size( pm );

            const unsigned n = 2;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets=true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets);

            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);

            eMesh.reopen();
            //eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields();
            stk::mesh::FieldBase *proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, 0);
            //eMesh.addParallelInfoFields(true,true);
            eMesh.commit();
            eMesh.set_proc_rank_field(proc_rank_field);

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            //eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            //eMesh.print_info("quad fixture",  2);
            eMesh.save_as(input_files_loc+"quad_4_1_smooth.0.e");

            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1"), eMesh.get_coordinates_field());

            unsigned center_node_id = 5;
            stk::mesh::Entity node = eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, center_node_id);
            double *data = 0;
            if (eMesh.is_valid(node))
              {
                data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                //std::cout << "tmp srk  center node= " << data[0] << " " << data[1] << std::endl;
                if (1)
                  {
                    data[0] += .2;
                    data[1] += .3;
                  }
              }

            eMesh.save_as(input_files_loc+"quad_4_1_smooth.0_perturbed.e");
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), eMesh.get_coordinates_field());

            {
              SMOOTHER pmmpsi(&eMesh, &boundarySelector, 0, 1001, 1.e-4, 1);
              pmmpsi.run(true, 0);
            }

            if (eMesh.is_valid(node))
              {
                EXPECT_DOUBLE_EQ_APPROX(data[0], 1.0);
                EXPECT_DOUBLE_EQ_APPROX(data[1], 1.0);
              }

            eMesh.save_as(output_files_loc+"quad_4_1_smooth.1.e");
          }
      }


      //=============================================================================
      //=============================================================================
      //=============================================================================

      // run PMMParallelShapeImprover on a domain with one side perturbed
      TEST(heavy_perceptMeshSmoother, tri_4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 2)
        if (1 || p_size <= par_size_max)
          {
            const unsigned nele = 12;
            //const unsigned nx = nele , ny = nele , nz = p_size*nele ;
            const unsigned nx = nele , ny = nele;

            bool sidesets=true;
            percept::QuadFixture<double, Triangle<3> > fixture( pm , nx , ny, sidesets);
            fixture.set_bounding_box(0,1,0,1);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            //eMesh.print_info("tri fixture",  2);
            eMesh.save_as(input_files_loc+"tri_4_smooth.0.e");

            eMesh.reopen();
            eMesh.add_coordinate_state_fields();
            eMesh.add_spacing_fields();
            stk::mesh::FieldBase *proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, 0);
            eMesh.commit();
            eMesh.set_proc_rank_field(proc_rank_field);

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( stk::topology::NODE_RANK );

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_elements_in_bucket = bucket.size();

                  for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity node = bucket[iEntity];

                      double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                      double iy = data[1]; // /double(nele);
                      iy = iy*iy;
                      data[1] = iy; // *double(nele);
                    }
                }
              }
            SpacingFieldUtil sfu(eMesh);
            sfu.compute_spacing_field();
            eMesh.save_as(input_files_loc+"tri_4_smooth.1.e");

            // save state of original mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1"), eMesh.get_coordinates_field());

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                if (boundarySelector_1(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity node = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , node );
                        double ix = data[0]; // /double(nele);
                        //double bump_size=2.8; // 0.8
                        double bump_size=2.8; // 0.8
                        data[1] += (ix)*(1.0-ix)*bump_size; //*double(nele);
                        //std::cout << "tmp srk surface 1 node = " << data[0] << " " << data[1] << std::endl;
                      }
                  }
              }

            // save state of projected mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"tri_4_smooth.0_perturbed.e");

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            {
              SMOOTHER pmmpsi(&eMesh, &boundarySelector, 0, 1001, 1.e-4, 1);
              pmmpsi.run( always_smooth, msq_debug);
            }

            eMesh.save_as(output_files_loc+"tri_4_si_smooth.1.e");


          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom, new parallel smoother

      void do_test_hex4(unsigned n=12, const std::string& filename = "hex_4_si_smooth.1.e")
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        unsigned par_size_max = s_par_size_max;

        const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= par_size_max)
          {
            std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for hex_1 case, n = " << n << std::endl;
            //unsigned nn = n+1;
            //std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n*p_size)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");

            bool test_edge_len = false;
            if (test_edge_len)
              {
                n=1;
              }

            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields();
            eMesh.commit();
            if (test_edge_len)
              {
                stk::mesh::Entity element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);

                double edge_length_ave = eMesh.edge_length_ave(element);
                std::cout << "edge_length_ave= " << edge_length_ave << std::endl;
                exit(1);
              }

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( stk::topology::NODE_RANK );

            // cluster the mesh towards the bump
            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (boundarySelector_5(**k))
                {
                  stk::mesh::Bucket & bucket = **k ;

                  const unsigned num_elements_in_bucket = bucket.size();

                  for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                    {
                      stk::mesh::Entity entity = bucket[iEntity];

                      double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                      data[2] = data[2]*data[2];
                    }
                }
              }
            eMesh.save_as(input_files_loc+"hex_4_smooth.0.e");

            // save state of original mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1"), eMesh.get_coordinates_field());

            // create the bump
            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                if (boundarySelector_5(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        double ix = data[0];
                        double iy = data[1];
                        data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*4.;
                      }
                  }
              }
            // save state of projected mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"hex_4_smooth.0_perturbed.e");
            eMesh.save_as("hex_uns_init.e");
            //eMesh.dump_vtk("hex_uns_init.vtk");

            std::cout << "tmp srk doing Shape smoothing for hex_4 case..." << std::endl;

            //bool do_jacobi = true;

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 1001;

            if (1)
              {
                SMOOTHER pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-5, 1);
                pmmpsi.m_do_animation = 1;
                pmmpsi.run( always_smooth, msq_debug);
              }

            eMesh.save_as(output_files_loc+filename);

          }
      }

      TEST(heavy_perceptMeshSmoother, hex_4)
      {
        do_test_hex4();
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom, new parallel smoother

      void do_hex_4_unstructured(unsigned nele = 3)
      {
        do_test_hex4(nele, "compare_hex_4_si_smooth.1.e");
      }

		void do_hex_4_structured(unsigned nele = 3)
		{
			stk::ParallelMachine pm = MPI_COMM_WORLD;
			MPI_Barrier( MPI_COMM_WORLD );

			const unsigned p_size = stk::parallel_machine_size( pm );
			if (p_size > 1) return;

			PerceptMesh eMesh(3u);

			const unsigned nxyz = nele+1;
			std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};

			std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
			eMesh.set_block_structured_grid(bsg);

			eMesh.add_coordinate_state_fields();

			typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
			VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
			typename StructuredGrid::MTField *coord_field_N = bsg->m_fields["coordinates_N"].get();
			VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
			typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
			VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

			// cluster nodes towards the bottom face
			using Node = typename StructuredGrid::MTNode;
			std::vector<Node> nodes;
			bsg->get_nodes(nodes);
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				data[2] = data[2]*data[2];
				set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
			}

			bsg->dump_vtk("hex_str_smooth_0");

			// save state of original mesh
			// field, dst, src:
			eMesh.copy_field(coord_field_NM1, coord_field);

			// now, make a bump on the bottom face, which pushes the face into the interior
			//   and thus creates a tangled mesh
			const double scale=2.0*4.*(5.0/(double)nele);
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				if (std::abs(data[2]) < 1e-5) // if (boundarySelector_5(node))
				{
					double ix = data[0];
					double iy = data[1];
					data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*scale;
					set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				}
			}

			// save state of projected mesh
			// field, dst, src:
			eMesh.copy_field(coord_field_N, coord_field);

			bsg->dump_vtk("hex_str_init");

			if (nele==5) {
			// test standalone total metric 2
				long double mtot=0;
				bool valid = false;
				size_t num_invalid=0;
				size_t n_invalid=0;

				SmootherMetricUntangleImpl<StructuredGrid> untangle_metric(&eMesh);

				GenericAlgorithm_total_element_metric<StructuredGrid> ga2(&untangle_metric, &eMesh, valid, &num_invalid, mtot, n_invalid);
				ga2.run();

				EXPECT_NEAR(mtot, 0.0179978, 1.e-6);
				EXPECT_EQ (n_invalid, (size_t)25);

				if (0==eMesh.get_parallel_rank()) {
					std::cout << "GenericAlgorithm_total_element_metric(SmootherMetricUntangleImpl): mtot = " << mtot
					<< " n_invalid = " << n_invalid << std::endl;
				}
			}

			int msq_debug = 0;// 1,2,3 for more debug info
			bool always_smooth = true;
			int innerIter = 1001;

			eMesh.setProperty("in_filename","messyMesh");
			stk::diag::Timer root_timer(eMesh.getProperty("in_filename"), rootTimerStructured());
			stk::diag::TimeBlock root_block(root_timer);

			{
				CubeBoundarySelector boundarySelector (&eMesh);

				ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-5, 1);
				pmmpsi.m_do_animation = 0;          	  	  //1;
				pmmpsi.run( always_smooth, msq_debug);
			}

			bsg->dump_vtk("hex_str_smooth_1");

			// generate unstructured smoothed case for comparison
			std::cout << "About to do_test_hex4"<<std::endl;
			do_test_hex4(nele, "compare_hex_4_si_smooth.1.e");

			PerceptMesh eMesh1;
			eMesh1.open_read_only(output_files_loc+"compare_hex_4_si_smooth.1.e");
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);

				stk::mesh::Entity u_node = eMesh1.get_node(data[0], data[1], data[2]);
				double * u_data = stk::mesh::field_data( *eMesh1.get_coordinates_field() , u_node );
				double d = Math::distance_3d(data, u_data);

				EXPECT_NEAR(d, 0.0, 1.e-12);
			}

		}

      TEST(heavy_perceptMeshSmoother, hex_4_sgrid)
      {
    	const std::vector<unsigned> nele_array = {5}; //{5,10,20,40};
    	for (auto nele : nele_array) {
    		do_hex_4_structured(nele);
    	}
      }

      TEST(heavy_perceptMeshSmoother, test_edge_len_sgrid_vs_ugrid)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size != 1)
          return;

        double edge_length_ave_uns = 0;
        double edge_length_ave_str = 0;

        double perturb[][3] = {{0.01,0.02,0.03},{-0.01,-0.03,0.015},{0.001,-0.002,-0.03},{-0.003,0.002,0.018},
                               {0.011,0.012,0.013},{-0.011,-0.013,0.0115},{0.0101,-0.0102,-0.013},{-0.0103,0.0102,0.0118}};


        // unstructured
        {
          unsigned n=1;
          std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
          PerceptMesh eMesh(3);
          eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
          eMesh.addParallelInfoFields(true,true);
          eMesh.add_coordinate_state_fields();
          eMesh.commit();

          stk::mesh::Entity element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);

          MyPairIterRelation elem_nodes(eMesh, element, eMesh.node_rank());
          for (int i=0; i < 8; ++i)
            {
              double *coord = stk::mesh::field_data(*eMesh.get_coordinates_field(), elem_nodes[i].entity());
              std::cout << "i= " << i << " coord= " << Math::print_3d(coord) << std::endl;
              for (int j=0; j < 3; ++j) coord[j] += perturb[i][j];
            }

          edge_length_ave_uns = eMesh.edge_length_ave(element);
          eMesh.dump_vtk("uns_cube.vtk");
          std::cout << "edge_length_ave_uns= " << edge_length_ave_uns << std::endl;
        }

        // structured
        {
          unsigned nele = 1;

          PerceptMesh eMesh(3u);

          unsigned nxyz = nele+1;
          std::array<unsigned, 3> sizes{{nxyz,nxyz,nxyz}};
          std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
          eMesh.set_block_structured_grid(bsg);
          eMesh.add_coordinate_state_fields();

          typename StructuredGrid::MTField *coord_field     = bsg->m_fields["coordinates"].get();
          VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
          typename StructuredGrid::MTField *coord_field_N   = bsg->m_fields["coordinates_N"].get();
          VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
          typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
          VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

          if (1)
            {
              StructuredCellIndex element{{0,0,0,0}};
              //void get_element_nodes(const StructuredCellIndex& element, std::vector<StructuredCellIndex>& nodes);

              typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
              std::vector<typename StructuredGrid::MTNode> nodes;
              bsg->get_element_nodes(element, nodes);
              for (int i = 0; i < 8; ++i)
                {
                  double nc[3];
                  int ii = BlockStructuredGrid::permute_to_unstructured[i];
                  //ii = i;
                  percept::get_field<StructuredGrid>(nc, 3, &eMesh, coord_field, nodes[i]);
                  for (int j = 0; j < 3; ++j)
                    {
                      nc[j] += perturb[ii][j];
                    }
                  percept::set_field<StructuredGrid>(nc, 3, &eMesh, coord_field, nodes[i]);
                }
              bsg->dump_vtk("bsg_cube");
              edge_length_ave_str = eMesh.edge_length_ave(element);
              std::cout << "edge_length_ave_str= " << edge_length_ave_str << std::endl;
            }
        }
        EXPECT_NEAR(edge_length_ave_str, edge_length_ave_uns, 1.e-8);
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom, new parallel smoother, convert to tet

      // heavy_perceptMeshSmoother.tet_4
      TEST(heavy_perceptMeshSmoother, tet_4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size <= 2)
          {
            unsigned n = 12;
            std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for tet_1 case, n = " << n << std::endl;
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields();
            int scalarDimension=0;
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            Hex8_Tet4_6_12 convert_pattern(eMesh);
            eMesh.commit();

            UniformRefiner converter(eMesh, convert_pattern, proc_rank_field);
            converter.doBreak();
            eMesh.save_as(input_files_loc+"tet_4_smooth.0.0.e");

            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;

            eMesh.populateParallelInfoFields(true,true,&boundarySelector);

            const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( stk::topology::NODE_RANK );

            // cluster the mesh towards the bump
            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (boundarySelector_5(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        data[2] = data[2]*data[2];
                      }
                  }
              }
            eMesh.save_as(input_files_loc+"tet_4_smooth.0.e");

            // save state of original mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1"), eMesh.get_coordinates_field());

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                if (boundarySelector_5(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        double ix = data[0];
                        double iy = data[1];
                        data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*4.;
                      }
                  }
              }
            // save state of projected mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"tet_4_smooth.0_perturbed.e");

            std::cout << "tmp srk doing Shape smoothing for tet_4 case..." << std::endl;

            //bool do_jacobi = true;

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 1001;

            {
              SMOOTHER pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-4, 1);
              pmmpsi.run( always_smooth, msq_debug);
            }

            eMesh.save_as(output_files_loc+"tet_4_si_smooth.1.e");

          }

      }
      
///////////////////////////////////////Madison's doodles

		double randDbl(double fMin, double fMax)
		{
			double f = (double)rand()/RAND_MAX;
			return fMin +f*(fMax-fMin);
		}

		void do_hex_4_structuredDOODLE(unsigned nele, double deformCoeff_1, double deformCoeff_2, int numIters)
		{
			stk::ParallelMachine pm = MPI_COMM_WORLD;
			MPI_Barrier( MPI_COMM_WORLD );

			const unsigned p_size = stk::parallel_machine_size( pm );
			if (p_size > 1) return;

			PerceptMesh eMesh(3u);

			const unsigned nxyz = nele+1;
			std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
			std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
			eMesh.set_block_structured_grid(bsg);

			eMesh.add_coordinate_state_fields();

			typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
			VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
			typename StructuredGrid::MTField *coord_field_N = bsg->m_fields["coordinates_N"].get();
			VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
			typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
			VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

			// cluster nodes towards the bottom face
			using Node = typename StructuredGrid::MTNode;
			std::vector<Node> nodes;
			bsg->get_nodes(nodes);
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				data[2] = data[2]*data[2];
				set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
			}

			bsg->dump_vtk("hex_str_smooth_0");

			// save state of original mesh
			// field, dst, src:
			eMesh.copy_field(coord_field_NM1, coord_field);

			// now, make a bump on the bottom face, which pushes the face into the interior
			//   and thus creates a tangled mesh
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				if (std::abs(data[2]) < 1e-5) // if (boundarySelector_5(node))
				{
					double ix = data[0];
					double iy = data[1];
					data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*deformCoeff_1*deformCoeff_2; //change scale to see if I can get more elements to run
					set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				}
			}

			// save state of projected mesh
			// field, dst, src:
			eMesh.copy_field(coord_field_N, coord_field);

			bsg->dump_vtk("hex_str_init");

			int innerIter = 1001;
			// selects all 6 faces of the cube
			struct CubeBoundarySelector : public StructuredGrid::MTSelector {
				PerceptMesh *m_eMesh;
				CubeBoundarySelector(PerceptMesh *eMesh) : m_eMesh(eMesh) {}

				bool operator()(StructuredCellIndex& index)
				{
					//commenting this out didn't help with poor scaling
					unsigned iblock = index[3];
					std::shared_ptr<StructuredBlock> sgrid = m_eMesh->get_block_structured_grid()->m_sblocks[iblock];
					unsigned sizes[3] = {sgrid->m_sizes.node_size[0], sgrid->m_sizes.node_size[1], sgrid->m_sizes.node_size[2]};

					if (index[0] == 0 || index[0] == sizes[0]-1 ||
							index[1] == 0 || index[1] == sizes[1]-1 ||
							index[2] == 0 || index[2] == sizes[2]-1 )
					{
						return true;
					}
					return false;

				}
			};

			eMesh.setProperty("in_filename","messyMesh");
			stk::diag::Timer root_timer(eMesh.getProperty("in_filename"), rootTimerStructured());
			stk::diag::TimeBlock root_block(root_timer);
			{
				CubeBoundarySelector boundarySelector (&eMesh); //to enforce consistency between blocks since there's apparently redundant nodes there, would I use this to keep them fixed?

				ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-5, 1);
				pmmpsi.m_do_animation = 0;//1;

			    long double alpha=0.1;
				GenericAlgorithm_update_coordinates<StructuredGrid> ga1(&pmmpsi,&eMesh,alpha);

				for(int i = 0; i<numIters;i++)//
				{ //Running this loop over the run fuction seemed to produce backward scaling results. I think this is due to excessive kokkos calls
					ga1.run();//focus just on running metric
				}
			}

			printTimersTableStructured();
		}

		TEST(/*DISABLED_*/heavy_perceptMeshSmoother, hex_4_sgridDOODLE)
		{
			unsigned nelebele = 100;
			double dc_1 = 2.0;
			double dc_2 = 4.0;
			int numIters=20;
			do_hex_4_structuredDOODLE(nelebele, dc_1, dc_2,numIters);
		}

#ifdef KOKKOS_HAVE_OPENMP
		void initMesh(unsigned nele, double deformCoeff_1, double deformCoeff_2, int numIters)
		{
			PerceptMesh eMesh(3u);
			const unsigned nxyz = nele+1;
			std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
//			std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(&eMesh, sizes);
			std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
			eMesh.set_block_structured_grid(bsg);

			eMesh.add_coordinate_state_fields();

			typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
			VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
			typename StructuredGrid::MTField *coord_field_N = bsg->m_fields["coordinates_N"].get();
			VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
			typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
			VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

//    	  if(!(bsg->m_fields.find("cg_s")== bsg->m_fields.end()))
//    		 std::cout <<"you do have a cg_s field" << std::endl;

			// cluster nodes towards the bottom face
			using Node = typename StructuredGrid::MTNode;
			std::vector<Node> nodes;
			bsg->get_nodes(nodes);
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				data[2] = data[2]*data[2];
				set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
			}

			// save state of original mesh
			// field, dst, src:
			eMesh.copy_field(coord_field_NM1, coord_field);

			// now, make a bump on the bottom face, which pushes the face into the interior
			//   and thus creates a tangled mesh
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				if (std::abs(data[2]) < 1e-5) // if (boundarySelector_5(node))
				{
					double ix = data[0];
					double iy = data[1];
					data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*deformCoeff_1*deformCoeff_2; //change scale to see if I can get more elements to run
					set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				}
			}

			// save state of projected mesh
			// field, dst, src:
			eMesh.copy_field(coord_field_N, coord_field);

			// now, make a bump on the bottom face, which pushes the face into the interior
			//   and thus creates a tangled mesh
			for (auto node : nodes) {
				double data[3];
				get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				if (std::abs(data[2]) < 1e-5) // if (boundarySelector_5(node))
				{
					double ix = data[0];
					double iy = data[1];
					data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*deformCoeff_1*deformCoeff_2; //change scale to see if I can get more elements to run
					set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
				}
			}

			eMesh.copy_field(coord_field_N, coord_field);

			//        int  msq_debug             = 0; // 1,2,3 for more debug info
			//        bool always_smooth         = true;
			int innerIter = 1001;
			// selects all 6 faces of the cube
			struct CubeBoundarySelector : public StructuredGrid::MTSelector {
				PerceptMesh *m_eMesh;
				CubeBoundarySelector(PerceptMesh *eMesh) : m_eMesh(eMesh) {}

				bool operator()(StructuredCellIndex& index)
				{
					//commenting this out didn't help with poor scaling
					unsigned iblock = index[3];
					std::shared_ptr<StructuredBlock> sgrid = m_eMesh->get_block_structured_grid()->m_sblocks[iblock];
					unsigned sizes[3] = {sgrid->m_sizes.node_size[0], sgrid->m_sizes.node_size[1], sgrid->m_sizes.node_size[2]};

					if (index[0] == 0 || index[0] == sizes[0]-1 ||
							index[1] == 0 || index[1] == sizes[1]-1 ||
							index[2] == 0 || index[2] == sizes[2]-1 )
					{
						return true;
					}
					return false;

				}
			};

			eMesh.setProperty("in_filename","messyMeshSimple");
			stk::diag::Timer root_timer(eMesh.getProperty("in_filename"), rootTimerStructured());
			stk::diag::TimeBlock root_block(root_timer);

			CubeBoundarySelector boundarySelector (&eMesh); //to enforce consistency between blocks since there's apparently redundant nodes there, would I use this to keep them fixed?

			ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-5, 1);
			pmmpsi.m_do_animation = 0;//1;

		    long double alpha=0.1;
			GenericAlgorithm_update_coordinates<StructuredGrid> ga1(&pmmpsi,&eMesh,alpha);

			simplified_gatm_1_BSG sgatm(&pmmpsi,&eMesh,alpha);

			for(int i=0;i< numIters;i++)
			{
				sgatm.run();
			}
			printTimersTableStructured();

		} //initMesh

		TEST(simplified, gatm1)
		{
			unsigned nele = 100;
			double dc_1 = 2.0;
			double dc_2 = 4.0;
			int numIters=5000;

			initMesh(nele,dc_1,dc_2,numIters);
		}
#endif

      void build_perturbed_cube(PerceptMesh& eMesh, const unsigned nele)
      {
        const unsigned nxyz = nele+1;
        std::array<unsigned, 3> sizes{{nxyz,nxyz,nxyz}};
        std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
        eMesh.set_block_structured_grid(bsg);
        eMesh.add_coordinate_state_fields();

        typename StructuredGrid::MTField *coord_field     = bsg->m_fields["coordinates"].get();
        VERIFY_OP_ON(coord_field,     !=, 0, "bad coord_field");
        typename StructuredGrid::MTField *coord_field_N   = bsg->m_fields["coordinates_N"].get();
        VERIFY_OP_ON(coord_field_N,   !=, 0, "bad coordinates_N");
        typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
        VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");
        
        // save state of original mesh
        // field, dst, src:
        eMesh.copy_field(coord_field_NM1, coord_field);
        
        using Node = typename StructuredGrid::MTNode;
        std::vector<Node> nodes;
        bsg->get_nodes(nodes);
        
        // randomly perturb interior nodes
        const double epsilon = 1.0/(double)nele;
        const double tol = 1e-6;
        std::srand(1);
        for (auto node : nodes) {
          double data[3];
          get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
          const double x = data[0], y = data[1], z = data[2];
          if (x > tol && x< 1-tol && 
              y > tol && y< 1-tol && 
              z > tol && z< 1-tol) {
            for (int d=0; d<3; d++) {
              data[d] = data[d] + epsilon*( (double)std::rand()/(double)RAND_MAX - 0.5);
            }
            set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
          }
        }
        
        //bsg->dump_vtk("hex_str_random");
        
        // save state of projected mesh
        // field, dst, src:
        eMesh.copy_field(coord_field_N, coord_field);
      }

      TEST(heavy_perceptMeshSmoother, total_metric_2)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;
        
        const std::vector<unsigned> nele_array = {100,200,300,400};

        const std::vector<size_t> n_invalid_array = {391365,     3189265, 10820052,   25719893};
        const std::vector<long double> mtot_array = {5.52566e-08,7.08e-09,2.10961e-09,8.9339e-10};

        for (unsigned i=0; i<1; i++) {

          const unsigned nele = nele_array[i];

          std::cout << "total_metric_2: nele = " << nele << std::endl;

          const size_t n_invalid_expected = n_invalid_array[i];
          const long double mtot_expected = mtot_array[i];

          PerceptMesh eMesh(3u);
        
          build_perturbed_cube(eMesh, nele);

/*
          std::vector<typename StructuredGrid::MTElement> elements;
          eMesh.get_block_structured_grid()->get_elements(elements);

          Kokkos::View<int*, DataLayout , MemSpace > element_invalid_flags;
          Kokkos::Experimental::resize(element_invalid_flags, elements.size());
*/

          int num_iterations = 10;
          for (int iter=0; iter<num_iterations; iter++)
          {
            // test standalone total metric
            //bool valid = false;
            //size_t num_invalid=0;
        	long double mtot=0;
            size_t n_invalid=0;
            
            //SmootherMetricUntangleImpl<StructuredGrid> untangle_metric(&eMesh);
            //SGridSmootherMetricUntangle untangle_metric(&eMesh);

            //GenericAlgorithm_total_element_metric<StructuredGrid> ga2(&untangle_metric, &eMesh, valid, &num_invalid, mtot, n_invalid);
            //ga2.run();

            SGridGenericAlgorithm_total_element_metric<SmootherMetricUntangleImpl<StructuredGrid> > ga2(&eMesh, mtot, n_invalid);
            ga2.run();

/*
            {
            	stk::diag::Timer reduce1(std::string("GATM ") + std::to_string(nele), rootTimerStructured());
                stk::diag::TimeBlock reduce1_block(reduce1);

                Kokkos::parallel_reduce(elements.size(), KOKKOS_LAMBDA(unsigned index, long double& mtot_loc)
                {
                  typename StructuredGrid::MTElement element = elements[index];

                  bool local_valid = false;
                  long double mm = const_cast<SmootherMetricUntangleImpl<StructuredGrid>*>(&untangle_metric)->metric(element, local_valid);
                  element_invalid_flags(index) = !local_valid;

                  mtot_loc += mm;
                }
                , mtot);
            }
            
            Kokkos::parallel_reduce(elements.size(), KOKKOS_LAMBDA(unsigned index, size_t& ninvalid_loc)
            {
            	ninvalid_loc += element_invalid_flags(index);
            }, n_invalid);
*/

            EXPECT_NEAR(mtot,      mtot_expected, 1.e-7);
            EXPECT_EQ  (n_invalid, n_invalid_expected);
            
            if (0 && 0==eMesh.get_parallel_rank() && iter==num_iterations-1) {
              std::cout << "GenericAlgorithm_total_element_metric(SmootherMetricUntangleImpl): mtot = " << mtot
                        << " n_invalid = " << n_invalid << std::endl;
            }
          }      
        }
      }

      struct Centroid {

        Centroid(typename StructuredGrid::MTField::Array4D coords_in, const unsigned nele_in, double &avgx_in) :
          coords(coords_in), nele(nele_in), avgx(avgx_in) {}
        
        typename StructuredGrid::MTField::Array4D coords;
        const unsigned nele;
        double &avgx;

        KOKKOS_INLINE_FUNCTION
        void 
        operator()(const unsigned& index, double& avgx_loc) const
        {
          // index = indx[0]*nele*nele + indx[1]*nele + indx[2]
          const unsigned indx[3] = {index / (nele*nele),
                                    (index / nele) % nele,
                                    index % nele};

          for (unsigned i=0; i<2; i++) {
            for (unsigned j=0; j<2; j++) {
              for (unsigned k=0; k<2; k++) {
                avgx_loc += coords(indx[0]+i,indx[1]+j,indx[2]+k,0);
              }
            }
          }
        }
        
        void run()
        {
          stk::diag::Timer root_timer("Centroid " + std::to_string(nele), rootTimerStructured());
          stk::diag::TimeBlock root_block(root_timer);

          const unsigned total_num_elements = nele*nele*nele;
          Kokkos::parallel_reduce(total_num_elements, *this, avgx);
          avgx /= (double)total_num_elements;
        }
      };

      TEST(heavy_perceptMeshSmoother, centroid)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;

        const std::vector<unsigned> nele_array = {100}; //{100,200,300,300,400};

        for (auto nele : nele_array) {

          PerceptMesh eMesh(3u);
        
          build_perturbed_cube(eMesh, nele);

          int num_iterations = 10;
          for (int iter=0; iter<num_iterations; iter++) {

            double avgx=0.0;
            
            const unsigned iblock = 0;
            
            typename StructuredGrid::MTField *coord_field = 
              eMesh.get_block_structured_grid()->m_fields["coordinates"].get();
            
            Centroid centroid(*coord_field->m_block_fields[iblock], nele, avgx);
            centroid.run();
            
            if ( iter+1==num_iterations) {
              std::cout << "avgx = " << avgx << " " << std::endl;
              
              EXPECT_NEAR(avgx, 4.00001, 1.e-5);
            }
          }
        }
      }

  }//heavy_tests
}//percept


#endif
#endif

