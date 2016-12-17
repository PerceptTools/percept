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
      TEST(heavy_perceptMeshSmoother, quad_4_1)
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
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields();
            eMesh.commit();

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

            std::cout << "tmp srk doing Shape smoothing for hex_4 case..." << std::endl;

            //bool do_jacobi = true;

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 1001;

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

      TEST(heavy_perceptMeshSmoother, hex_4_structured)
      {
        unsigned nele = 3;
        //if (true) return;
        bool doCompare = true;
        if (doCompare)
          do_test_hex4(nele, "compare_hex_4_si_smooth.1.e");

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );
        //unsigned par_size_max = s_par_size_max;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        //if (p_size == 1 || p_size == 3)
        if (p_size == 1)
          {
            PerceptMesh eMesh(3u);

            unsigned nxyz = nele+1;
            std::array<unsigned, 3> sizes{{nxyz,nxyz,nxyz}};
            std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(&eMesh, sizes);
            eMesh.add_coordinate_state_fields();

            bsg->dump_vtk("bsg_cube");

            //StructuredCellIndex element{{0,0,0,0}};

            typename StructuredGrid::MTField *coord_field     = bsg->m_fields["coordinates"].get();
            VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
            typename StructuredGrid::MTField *coord_field_N   = bsg->m_fields["coordinates_N"].get();
            VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
            typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
            VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

#if 0
            stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
            stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
            stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
            stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
            stk::mesh::Selector boundarySelector_5(*eMesh.get_non_const_part("surface_5") );
            stk::mesh::Selector boundarySelector_6(*eMesh.get_non_const_part("surface_6") );
            stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4 | boundarySelector_5 | boundarySelector_6;
#endif
            using Node = typename StructuredGrid::MTNode;
            std::vector<Node> nodes;
            bsg->get_nodes(nodes);
            for (auto node : nodes)
              {
                double data[3];
                get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
                data[2] = data[2]*data[2];
                set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
              }
            bsg->dump_vtk("hex_str_smooth_0");

            // save state of original mesh
            // field, dst, src:
            eMesh.copy_field(coord_field_NM1, coord_field);

            if (0)
              {
                std::vector<typename StructuredGrid::MTElement> elements;
                bsg->get_elements(elements);
                SmootherMetricUntangleImpl<StructuredGrid> smu(&eMesh);
                std::cout << "elem0= " << elements[0] << std::endl;
                for (auto elem : elements)
                  {
                    bool valid = true;
                    smu.metric(elem, valid);
                    //std::cout << "elem= " << elem << " valid= " << valid << std::endl;
                  }
              }

            // using Element = typename StructuredGrid::MTElement;
            // std::vector<Element> elements;
            // bsg->get_elements(elements);

            for (auto node : nodes)
              {
                //if (boundarySelector_5(node))
                {
                  double data[3];
                  get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
                  if (std::abs(data[2]) < 1e-5)
                    {
                      double ix = data[0];
                      double iy = data[1];
                      data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*2.0*4.;
                      set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
                    }
                }
              }

            // save state of projected mesh
            // field, dst, src:
            eMesh.copy_field(coord_field_N, coord_field);

            bsg->dump_vtk("hex_str_smooth_perturbed");

            std::cout << "tmp srk doing Shape smoothing for hex_4 case..." << std::endl;

            //bool do_jacobi = true;

            int  msq_debug             = 0; // 1,2,3 for more debug info
            bool always_smooth         = true;
            int innerIter = 1001;

            struct CubeBoundarySelector : public StructuredGrid::MTSelector {
              PerceptMesh *m_eMesh;
              CubeBoundarySelector(PerceptMesh *eMesh) : m_eMesh(eMesh){}

              bool operator()(StructuredCellIndex& index)
              {
                unsigned iblock = index[3];
                std::shared_ptr<StructuredBlock> sgrid = m_eMesh->get_block_structured_grid()->m_sblocks[iblock];
                //const unsigned A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];
                unsigned sizes[3] = {sgrid->m_sizes.node_size[0], sgrid->m_sizes.node_size[1], sgrid->m_sizes.node_size[2]};
                //unsigned Asizes[3] = {m_sizes.node_size[A0], m_sizes.node_size[A1], m_sizes.node_size[A2]};

                if (index[0] == 0 || index[0] == sizes[0]-1 ||
                    index[1] == 0 || index[1] == sizes[1]-1 ||
                    index[2] == 0 || index[2] == sizes[2]-1 )
                  {
                    return true;
                  }
                return false;
              }
            };

            if (1)
              {
                CubeBoundarySelector boundarySelector (&eMesh);
                ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-5, 1);
                pmmpsi.m_do_animation = 1;
                pmmpsi.run( always_smooth, msq_debug);
              }

            bsg->dump_vtk("hex_str_smooth_1");

            if (doCompare)
              {
                PerceptMesh eMesh1;
                eMesh1.open_read_only(output_files_loc+"compare_hex_4_si_smooth.1.e");
                // std::vector<stk::mesh::Entity> vec;
                // stk::mesh::get_selected_entities(eMesh1.get_fem_meta_data()->universal_part() , eMesh1.get_bulk_data()->buckets(eMesh.node_rank()), vec);
                      /// find and return pointer to node closest to given point - in parallel, check return for null (if null, closest node is on another proc)
                for (auto node : nodes)
                  {
                    double data[3];
                    get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);

                    stk::mesh::Entity u_node = eMesh1.get_node(data[0], data[1], data[2]);
                    double * u_data = stk::mesh::field_data( *eMesh1.get_coordinates_field() , u_node );
                    double d = Math::distance_3d(data, u_data);
                    //std::cout << "d= " << d << std::endl;
                    EXPECT_NEAR(d, 0.0, 1.e-5);
                  }
              }
          }
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




    }

  }

#endif
#endif
