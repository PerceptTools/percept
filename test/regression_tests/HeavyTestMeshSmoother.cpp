// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#include <percept/Percept.hpp>
#include "kokkos_misuse_bugs_isolation.hpp"
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
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherAlgebraic.hpp>
#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_total_element_metric.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_update_coordinates.hpp>
#include <percept/mesh/mod/smoother/gradient_functors.hpp>
#include <percept/mesh/mod/smoother/get_alpha_0_refmesh.hpp>
#include <percept/mesh/mod/smoother/get_edge_len_avg.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>
#include <percept/structured/StructuredGridRefiner.hpp>

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

#include <percept/PerceptUtils.hpp>
#include <percept/Util.hpp>

#include <percept/math/DenseMatrix.hpp>

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



      void build_perturbed_cube(PerceptMesh& eMesh, const unsigned nele, bool add_all_fields=false)
        {
            std::cout << "about to build unperturbed cube\n";
            const unsigned nxyz = nele+1;
            std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
            std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
            eMesh.set_block_structured_grid(bsg);

            if(!add_all_fields) {
                bsg->register_field("coordinates_N", eMesh.get_spatial_dim());
                bsg->register_field("coordinates_NM1", eMesh.get_spatial_dim());
                bsg->register_field("cg_g",eMesh.get_spatial_dim());
            }
            else
            eMesh.add_coordinate_state_fields();

            typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
            VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
            typename StructuredGrid::MTField *coord_field_N = bsg->m_fields["coordinates_N"].get();
            VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
            typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
            VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

            // save state of original mesh
            // field, dst, src:
            eMesh.copy_field(coord_field_NM1, coord_field);

            Array4D coord_field_0 = *(coord_field->m_block_fields[0]);//hard coded to single block case
            Array4D::HostMirror host_coord_field = Kokkos::create_mirror_view(coord_field_0);
            Kokkos::deep_copy(host_coord_field,coord_field_0);

            // randomly perturb interior nodes
            const double epsilon = 1.0/(double)nele;
            const double tol = 1e-6;
            std::srand(1);
            double data[3];

            std::vector<std::array<double, 3> > rand_coords;
            std::vector<double> rand_coeffs;

            std::cout << "about to perturbed cube\n";
            for (unsigned k=0; k<sizes[0]; k++) {
                for (unsigned j=0; j<sizes[1]; j++) {
                    for (unsigned i=0; i<sizes[2]; i++) {
                        for (int d=0; d<3; d++) {
                            data[d]=host_coord_field(i,j,k,d);
                        }
                        const double x = data[0], y = data[1], z = data[2];
                        if (x > tol && x< 1-tol &&
                                y > tol && y< 1-tol &&
                                z > tol && z< 1-tol) {
                            for (int d=0; d<3; d++) {
                                data[d] = data[d] +epsilon*( (double)std::rand()/(double)RAND_MAX - 0.5);

//                                std::array<double,3> maniped_data = { {data[0],data[1],data[2]}};
//                                double rand_coeff = (double)std::rand()/(double)RAND_MAX;
//                                rand_coords.push_back(maniped_data);
//                                rand_coeffs.push_back(rand_coeff);

                            }
                        }
                        for (int d=0; d<3; d++) {
                            host_coord_field(i,j,k,d)=data[d];
                        }
                    }
                }
            }
            std::cout << "coordinates perturbed\n";

            Kokkos::deep_copy(coord_field_0, host_coord_field);
            std::cout <<"deep copy completed ...";
            //bsg->dump_vtk("hex_str_random");

            // save state of projected mesh
            // field, dst, src:
            eMesh.copy_field(coord_field_N, coord_field);
            std::cout << "fields copied\n";

//            std::cout << "Size of rand_coeffs = " << rand_coeffs.size() << " size of rand_coords = " << rand_coords.size() << std::endl;
//            for(unsigned iRand=0;iRand<rand_coords.size();iRand++)
//            {
//                std::cout << "(x, y, z,) = (" << rand_coords[iRand][0] << ", " << rand_coords[iRand][1] << ", " << rand_coords[iRand][2] << ") "
//                << " and coeff = " << rand_coeffs[iRand] << std::endl;
//            }
        }//build_perturbed_cube

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

            eMesh.reopen();
            eMesh.add_coordinate_state_fields(false);
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

        SMOOTHER pmmpsi(&eMesh, &boundarySelector, NULL, 0, 1001, 1.e-4, 1);
        pmmpsi.run();

        eMesh.save_as(output_files_loc+"quad_4_si_smooth.1.e");
        delete &eMesh;
      }

      TEST(heavy_perceptMeshSmoother, quad_4_Algebraic)
      {
        PerceptMesh& eMesh = *setup_quad_4_tests();

        stk::mesh::Selector boundarySelector_1(*eMesh.get_non_const_part("surface_1") );
        stk::mesh::Selector boundarySelector_2(*eMesh.get_non_const_part("surface_2") );
        stk::mesh::Selector boundarySelector_3(*eMesh.get_non_const_part("surface_3") );
        stk::mesh::Selector boundarySelector_4(*eMesh.get_non_const_part("surface_4") );
        stk::mesh::Selector boundarySelector = boundarySelector_1 | boundarySelector_2 | boundarySelector_3 | boundarySelector_4;

        int do_anim = 1;

        double drop_off_coeffs[3] = {1,1,1};  // FIXME
        int nlayers_drop_off = 10;
        int niter = 1;
        percept::ReferenceMeshSmootherAlgebraic pmmpsi(&eMesh, &boundarySelector,NULL, 0, niter, 1.e-4, drop_off_coeffs, nlayers_drop_off);
        pmmpsi.m_do_animation = do_anim;
        pmmpsi.run();

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
            eMesh.add_coordinate_state_fields(true);
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
              SMOOTHER pmmpsi(&eMesh, &boundarySelector,NULL, 0, 1001, 1.e-4, 1);
              pmmpsi.run();
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
            eMesh.add_coordinate_state_fields(true);
            eMesh.add_spacing_fields(true);
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

            {
              SMOOTHER pmmpsi(&eMesh, &boundarySelector,NULL, 0, 1001, 1.e-4, 1);
              pmmpsi.run();
            }

            eMesh.save_as(output_files_loc+"tri_4_si_smooth.1.e");


          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      // A cube with an indented bump on the bottom, new parallel smoother

      void do_test_hex4(unsigned n=6, const std::string& filename = "hex_4_si_smooth.1.e",const double smooth_tol=1e-5)
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
            eMesh.add_coordinate_state_fields(true);
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
            const double scale=2.0*4.*(5.0/(double)(n));
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
//                if (boundarySelector_5(**k))
//                  {
                    stk::mesh::Bucket & bucket = **k ;

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iEntity = 0; iEntity < num_elements_in_bucket; iEntity++)
                      {
                        stk::mesh::Entity entity = bucket[iEntity];

                        double * data = stk::mesh::field_data( *eMesh.get_coordinates_field() , entity );
                        double ix = data[0];
                        double iy = data[1];
                        if (std::abs(data[2]) < 1e-5) {
                            data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*scale;
                        }
                      }
//                  }
              }
            // save state of projected mesh
            // field, dst, src:
            eMesh.copy_field(eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), eMesh.get_coordinates_field());

            eMesh.save_as(input_files_loc+"hex_4_smooth.0_perturbed.e");
            eMesh.save_as("hex_uns_init.e");

            int innerIter = 1001;

            if (1)
              {
                SMOOTHER pmmpsi(&eMesh, &boundarySelector,NULL, 0, innerIter, 1.e-5, 1);
                pmmpsi.m_do_animation = 0;
                pmmpsi.run();
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

      void make_different_fixtures()
      {
          stk::ParallelMachine pm = MPI_COMM_WORLD;
          MPI_Barrier( MPI_COMM_WORLD );

          const unsigned p_size = stk::parallel_machine_size( pm );
          if (p_size > 1) return;

          PerceptMesh eMesh(3u);

          std::array<unsigned, 3> sizes_b0 { {5+1,5+1,1+1}};//nele+1
          std::array<double, 3> widths_b0 { {1,1,.2}};
          std::array<double, 3> offsets_b0 { {0.0,0.0,0.0}};

          std::shared_ptr<BlockStructuredGrid> bsg_b0 = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes_b0, widths_b0, offsets_b0);
          bsg_b0->dump_vtk("test_b0");

          std::array<unsigned, 3> sizes_b1 { {5+1,1+1,4+1}};//nele+1
          std::array<double, 3> widths_b1 { {1,.2,.8}};
          std::array<double, 3> offsets_b1 { {0.0,0.0,0.2}};

          std::shared_ptr<BlockStructuredGrid> bsg_b1 = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes_b1, widths_b1, offsets_b1);
          bsg_b1->dump_vtk("test_b1");

          std::array<unsigned, 3> sizes_b2 { {5+1,4+1,4+1}};//nele+1
          std::array<double, 3> widths_b2 { {1,.8,.8}};
          std::array<double, 3> offsets_b2 { {0.0,0.2,0.2}};

          std::shared_ptr<BlockStructuredGrid> bsg_b2 = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes_b2, widths_b2, offsets_b2);
          bsg_b2->dump_vtk("test_b2");
      }

      TEST(DISABLED_diff,fixt)
      {
          make_different_fixtures();
      }

      void do_hex_4_structured(unsigned nele = 3, const double smooth_tol=1e-5)
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

//NEW
            const double scale=2.0*4.*(5.0/(double)nele);

            Array4D coord_field_0 = *(coord_field->m_block_fields[0]);//hard coded to single block case
            Array4D::HostMirror host_coord_field = Kokkos::create_mirror_view(coord_field_0);
            Kokkos::deep_copy(host_coord_field,coord_field_0);

            std::cout << "about to perturbed cube\n";
            for (unsigned k=0; k<sizes[0]; k++) {
                for (unsigned j=0; j<sizes[1]; j++) {
                    for (unsigned i=0; i<sizes[2]; i++) {
                        double data[3];

                        for (int d=0; d<3; d++) {
                            data[d]=host_coord_field(i,j,k,d);
                        }
                        data[2]=data[2]*data[2];
                        for (int d=0; d<3; d++) {
                            host_coord_field(i,j,k,d)=data[d];
                        }
                    }
                }
            }

            Kokkos::deep_copy(coord_field_0, host_coord_field);
            eMesh.copy_field(coord_field_NM1, coord_field);

            for (unsigned k=0; k<sizes[0]; k++) {
                for (unsigned j=0; j<sizes[1]; j++) {
                    for (unsigned i=0; i<sizes[2]; i++) {
                        double data[3];

                        for (int d=0; d<3; d++) {
                            data[d]=host_coord_field(i,j,k,d);
                        }
                        double ix = data[0];
                        double iy = data[1];

                        if (std::abs(data[2]) < 1e-5) {
                            data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*scale;
                        }
                        for (int d=0; d<3; d++) {
                            host_coord_field(i,j,k,d)=data[d];
                        }
                    }
                }
            }

            std::cout << "coordinates perturbed\n";

            Kokkos::deep_copy(coord_field_0, host_coord_field);
            std::cout <<"deep copy completed ...";
//NEW
            bsg->dump_vtk("hex_str_smooth_0_"+std::to_string(nele)+"x"+std::to_string(nele)+"x"+std::to_string(nele));

            eMesh.copy_field(coord_field_N, coord_field);

            bsg->dump_vtk("hex_str_init_"+std::to_string(nele)+"x"+std::to_string(nele)+"x"+std::to_string(nele));

            if (0/*nele==5*/) {//fails gpu because it is an old metric
                // test standalone total metric 2
                Double mtot=0;
                bool valid = false;
                size_t num_invalid=0;
                size_t n_invalid=0;

                SmootherMetricUntangleImpl<StructuredGrid> untangle_metric(&eMesh);

                GenericAlgorithm_total_element_metric<StructuredGrid> ga2(&untangle_metric, &eMesh, valid, &num_invalid, mtot, n_invalid);
                ga2.run(0);//hard coded to signel block case

                mtot=ga2.mtot;
                n_invalid=ga2.n_invalid;

                EXPECT_NEAR(mtot, 0.0179978, 1.e-6);
                EXPECT_EQ (n_invalid, (size_t)25);

                if (0==eMesh.get_parallel_rank()) {
                    std::cout << "GenericAlgorithm_total_element_metric(SmootherMetricUntangleImpl): mtot = " << mtot
                    << " n_invalid = " << n_invalid << std::endl;
                }
            }

            int innerIter = 1001;

            eMesh.setProperty("in_filename","messyMesh_structured__"+std::to_string(nele)+"x"+std::to_string(nele)+"x"+std::to_string(nele));
            stk::diag::Timer root_timer(eMesh.getProperty("in_filename"), rootTimerStructured());
            stk::diag::TimeBlock root_block(root_timer);

            {
                SGridBoundarySelector boundarySelector (&eMesh);

                ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);
                pmmpsi.m_do_animation = 0;                    //1;
                pmmpsi.run();
                std::cout << "\n\n~~finished smoother for structured grid~~\n\n";
            }

            bsg->dump_vtk("hex_str_smooth_1_"+std::to_string(nele)+"x"+std::to_string(nele)+"x"+std::to_string(nele));

            // generate unstructured smoothed case for comparison
            std::cout << "About to do_test_hex4"<<std::endl;
            do_test_hex4(nele, "compare_hex_4_si_smooth.1.e",smooth_tol);

            PerceptMesh eMesh1;
            eMesh1.open_read_only(output_files_loc+"compare_hex_4_si_smooth.1.e");

            unsigned noBad=0;
            double maxDist=0;
            coord_field_0 = *(coord_field->m_block_fields[0]);
            Kokkos::deep_copy(host_coord_field,coord_field_0);

            const SGridSizes block_sizes = bsg->m_sblocks[0]->m_sizes;
            const std::array<unsigned, 3> loop_orderings =
                    bsg->m_sblocks[0]->m_loop_ordering;


            for (/*auto node : nodes*/unsigned iNode=0;iNode<nodes.size();iNode++) {
                double data[3];
//                get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);

                unsigned indx[3];
                const int L0 = loop_orderings[0], L1 =
                        loop_orderings[1], L2 = loop_orderings[2];
                const unsigned int size2[3] = { 1 + block_sizes.node_max[L0]
                        - block_sizes.node_min[L0], 1
                        + block_sizes.node_max[L1]
                        - block_sizes.node_min[L1], 1
                        + block_sizes.node_max[L2]
                        - block_sizes.node_min[L2] };
                indx[L2] = block_sizes.node_min[L2]
                        + (iNode / (size2[L0] * size2[L1]));
                indx[L1] = block_sizes.node_min[L1]
                        + ((iNode / size2[L0]) % size2[L1]);
                indx[L0] = block_sizes.node_min[L0] + (iNode % size2[L0]);

                for(unsigned iDim=0;iDim<3;iDim++)
                { data[iDim]=host_coord_field(indx[0],indx[1],indx[2],iDim);}

                stk::mesh::Entity u_node = eMesh1.get_node(data[0], data[1], data[2]);
                double * u_data = stk::mesh::field_data( *eMesh1.get_coordinates_field() , u_node );
                double d = Math::distance_3d(data, u_data);

              //  EXPECT_NEAR(d, 0.0, 1.e-12);// tolerance of 1e-12 seems to be too high if Double = double instead of Double = long double; test fails if it's 1e-12 and Double = double
                if(std::abs(d)>maxDist)
                    maxDist=std::abs(d);
                if(d>1e-11/*1e-12*/ || d<-1e-11/*-1e-12*/)//madbrew:    truncation error in parallel adds causes this to fail the diff under a tolerance of 1e-12
                    noBad++;
            }

        std::ofstream dump;
        dump.open("post_smooth_results.txt", ios::app);
        dump << std::endl;
        dump << std::endl;
        dump << std::endl;
        dump << "Tolerance Used on smoother" << smooth_tol << std::endl;
        dump << "Size of mesh" << nele*nele*nele <<std::endl;
        dump << "Number coords out of tolerance: " << noBad << " out of "<< nodes.size()<< std::endl;
        dump << "Max dissonance between nodes is  " << maxDist << std::endl;
        dump << std::endl;
        dump << std::endl;
        dump << std::endl;
        dump.close();

        EXPECT_EQ  (noBad, (unsigned)0);
        }//do_hex_4_structured

      TEST(heavy_perceptMeshSmoother, hex_4_sgrid)
      {
        const std::vector<unsigned> nele_array = {5};//{5,10,20,40};
        const std::vector<double> smoother_tol_arr = {1e-5};//{1e-5,1e-6,1e-7};
        for(auto tol : smoother_tol_arr){
            for (auto nele : nele_array) {
                do_hex_4_structured(nele,tol);
            }
        }
      }

      TEST(heavy_perceptMeshSmoother, test_edge_len_sgrid_vs_ugrid)
      {
#ifdef KOKKOS_ENABLE_CUDA
        return;
#endif
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
          eMesh.add_coordinate_state_fields(true);
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
            unsigned n = 6;
            std::cout << "P["<<p_rank<<"] " << "tmp srk doing Laplace smoothing for tet_1 case, n = " << n << std::endl;
            std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
            PerceptMesh eMesh(3);
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.addParallelInfoFields(true,true);
            eMesh.add_coordinate_state_fields(true);
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

            int innerIter = 1001;

            {
              SMOOTHER pmmpsi(&eMesh, &boundarySelector,NULL, 0, innerIter, 1.e-4, 1);
              pmmpsi.run();
            }

            eMesh.save_as(output_files_loc+"tet_4_si_smooth.1.e");

          }

      }


        void do_hex_4_structured_isolate_update_nodes(unsigned nele /*, double deformCoeff_1, double deformCoeff_2,*/, int numIters)
        {
            stk::ParallelMachine pm = MPI_COMM_WORLD;
            MPI_Barrier( MPI_COMM_WORLD );

            const unsigned p_size = stk::parallel_machine_size( pm );
            if (p_size > 1) return;

            PerceptMesh eMesh(3u);
            build_perturbed_cube(eMesh,nele,true);

            printf("setting uo boundary selector\n");
            int innerIter = 1001;
            // selects all 6 faces of the cube
            struct SGridBoundarySelector : public StructuredGrid::MTSelector {
                PerceptMesh *m_eMesh;
                SGridBoundarySelector(PerceptMesh *eMesh_in) : m_eMesh(eMesh_in) {}

                bool operator()(StructuredCellIndex& index)
                {

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
                SGridBoundarySelector boundarySelector (&eMesh);
                printf("setting up update_coords\n");
                ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, NULL,&boundarySelector, 0, innerIter, 1.e-5, 1);
                pmmpsi.m_do_animation = 0;//1;

                /*long double*/Double alpha=0.1;
                GenericAlgorithm_update_coordinates<StructuredGrid> ga1(&pmmpsi,&eMesh,alpha);

                printf("running update\n");
                for(int i = 0; i<numIters;i++)//
                {
                    ga1.run();
                }
            }

            Kokkos::View<Double*,DataLayout,MemSpace> v("v",70000);
            Kokkos::View<Double*,DataLayout,MemSpace>::HostMirror vm;
            vm = Kokkos::create_mirror_view(v);

            for(unsigned i=0;i<v.size();i++)
            {
                vm(i)=90.8;
                if(i==4578)
                    vm(i)=10000.0;
            }

            Kokkos::deep_copy(v,vm);
            max_scanner<Double> ms(v);
            std::cout << ms.find_max() << std::endl;


        }

        TEST(heavy_perceptMeshSmoother, hex_4_sgrid_isolate_update_nodes)
        {
            unsigned nelebele = 10;
            int numIters=200;
            do_hex_4_structured_isolate_update_nodes(nelebele/*, dc_1, dc_2,*/,numIters);
        }

#ifdef KOKKOS_ENABLE_OPENMP
        void test_simplified_update_nodes(unsigned nele, double deformCoeff_1, double deformCoeff_2, int numIters)
        {
            PerceptMesh eMesh(3u);
            const unsigned nxyz = nele+1;
            std::array<unsigned, 3> sizes { {nxyz,nxyz,nxyz}};
//          std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(&eMesh, sizes);
            std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
            eMesh.set_block_structured_grid(bsg);

            eMesh.add_coordinate_state_fields();

            typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
            VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");
            typename StructuredGrid::MTField *coord_field_N = bsg->m_fields["coordinates_N"].get();
            VERIFY_OP_ON(coord_field_N, !=, 0, "bad coordinates_N");
            typename StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
            VERIFY_OP_ON(coord_field_NM1, !=, 0, "bad coordinates_NM1");

//        if(!(bsg->m_fields.find("cg_s")== bsg->m_fields.end()))
//           std::cout <<"you do have a cg_s field" << std::endl;

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
                    data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*deformCoeff_1*deformCoeff_2;
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
                    data[2] = (ix)*(1.0-ix)*(iy)*(1.0-iy)*deformCoeff_1*deformCoeff_2;
                    set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
                }
            }

            eMesh.copy_field(coord_field_N, coord_field);

            //        int  msq_debug             = 0; // 1,2,3 for more debug info
            //        bool always_smooth         = true;
            int innerIter = 1001;
            // selects all 6 faces of the cube
            struct SGridBoundarySelector : public StructuredGrid::MTSelector {
                PerceptMesh *m_eMesh;
                SGridBoundarySelector(PerceptMesh *eMesh_) : m_eMesh(eMesh_) {}

                bool operator()(StructuredCellIndex& index)
                {
                    //commenting this out didn't help with poor scaling
                    unsigned iblock = index[3];
                    std::shared_ptr<StructuredBlock> sgrid = m_eMesh->get_block_structured_grid()->m_sblocks[iblock];
                    unsigned sizes_[3] = {sgrid->m_sizes.node_size[0], sgrid->m_sizes.node_size[1], sgrid->m_sizes.node_size[2]};

                    if (index[0] == 0 || index[0] == sizes_[0]-1 ||
                            index[1] == 0 || index[1] == sizes_[1]-1 ||
                            index[2] == 0 || index[2] == sizes_[2]-1 )
                    {
                        return true;
                    }
                    return false;
                }
            };

            eMesh.setProperty("in_filename","messyMeshSimple");
            stk::diag::Timer root_timer(eMesh.getProperty("in_filename"), rootTimerStructured());
            stk::diag::TimeBlock root_block(root_timer);

            SGridBoundarySelector boundarySelector (&eMesh);

            ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, NULL, &boundarySelector, 0, innerIter, 1.e-5, 1);
            pmmpsi.m_do_animation = 0;//1;

            /*long double*/Double alpha=0.1;
            GenericAlgorithm_update_coordinates<StructuredGrid> ga1(&pmmpsi,&eMesh,alpha);

            simplified_gatm_1_BSG sgatm(&pmmpsi,&eMesh,alpha);

            for(int i=0;i< numIters;i++)
            {
                sgatm.run();
            }
            printTimersTableStructured();

        } //test_simplified_update_nodes

        TEST(DISABLED_heavy_perceptMeshSmoother, simplfied_update_nodes)
        {
            unsigned nele = 100;
            double dc_1 = 2.0;
            double dc_2 = 4.0;
            int numIters=5000;

            test_simplified_update_nodes(nele,dc_1,dc_2,numIters);
        }
#endif
        struct zero_a4d
        {
            StructuredGrid::MTField::Array4D a4d;

            zero_a4d(StructuredGrid::MTField::Array4D a4d_in){a4d=a4d_in;}

            void zero_out_fields()
            {
                Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0,(unsigned)a4d.size()),*this);
            }

            KOKKOS_INLINE_FUNCTION
            void operator()(const unsigned& index) const
            {
                a4d.data()[index]=0;
            }

        };

      TEST(heavy_perceptMeshSmoother, total_metric_and_gradient)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;
        
        const std::vector<unsigned> nele_array = {100,200,300,400};

        const std::vector<size_t> n_invalid_array = {391365,     3189265, 10820052,   25719893};
        const std::vector<Double> mtot_array      = {5.52566e-08,7.08e-09,2.10961e-09,8.9339e-10};

        const std::vector<double> grad_untangle     = {1.82363504678e-07, 1.63634415837e-08, 3.97221589862e-09, 1.45422487689e-09};
        const std::vector<double> grad_smooth       = {389575833.965,     1124552264.63,     89552984978.3,     93567250196.3};

        const bool do_untang=true;

        const unsigned lower=0;
        const unsigned upper=1;

        for (unsigned iMeshSize=lower; iMeshSize<upper; iMeshSize++) {

          const unsigned nele = nele_array[iMeshSize];

          std::cout << "total_metric_2: nele = " << nele << std::endl;

          const size_t n_invalid_expected = n_invalid_array[iMeshSize];
          const Double mtot_expected = mtot_array[iMeshSize];

          PerceptMesh eMesh(3u);

#if !defined(KOKKOS_ENABLE_CUDA)
          bool add_all_fields = true;
#else
          bool add_all_fields = false;
#endif
          build_perturbed_cube(eMesh, nele, add_all_fields);
          StructuredGrid::MTField::Array4D grad_fld_data = *eMesh.get_block_structured_grid()->m_fields["cg_g"].get()->m_block_fields[0];
          StructuredGrid::MTField * grad_field           =  eMesh.get_block_structured_grid()->m_fields["cg_g"].get();

          sGrid_GenericAlgorithm_get_gradient_1 gg_1_sgrid (&eMesh, NULL);

          int grad_num_iters = 1;
          for(int iter=0;iter<grad_num_iters;iter++){

            // grad for smoothing
            gg_1_sgrid.m_metric.m_untangling=false;
            gg_1_sgrid.run();
            double g = std::sqrt(eMesh.nodal_field_dot(grad_field, grad_field));
            //std::cout << std::setprecision(12) << "grad smoothing = " << g << std::endl;
            EXPECT_NEAR(g,grad_smooth[iMeshSize],1.0);

            // grad for untangling
            gg_1_sgrid.m_metric.m_untangling=true;
            zero_a4d za_1(grad_fld_data);
            za_1.zero_out_fields();
            gg_1_sgrid.run();
            g = std::sqrt(eMesh.nodal_field_dot(grad_field, grad_field));
            //std::cout << std::setprecision(12) << "grad untangling = " << g << std::endl;
            EXPECT_NEAR(g,grad_untangle[iMeshSize],1e-6);

            za_1.zero_out_fields();
          }

          int metric_num_iterations = 1;
          for (int iter=0; iter<metric_num_iterations; iter++)
          {
            Double mtot=0;
            size_t n_invalid=0;

            SGridGenericAlgorithm_total_element_metric< HexMeshSmootherMetric >  ga_tem(&eMesh, mtot, n_invalid);
            ga_tem.m_metric.m_untangling = do_untang;
            ga_tem.m_metric.m_use_ref_mesh = true;

            ga_tem.run(0);

            mtot = ga_tem.mtot;
            n_invalid=ga_tem.n_invalid;

            EXPECT_NEAR(mtot,      mtot_expected, 1.e-7);
            EXPECT_EQ  (n_invalid, n_invalid_expected);

          }//iterations metric
        }//mesh sizes
      }

      TEST(heavy_perceptMeshSmoother, sgrid_perturbed_smooth)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;

        const std::vector<unsigned> nele_array = {25,50,100,200};

        const unsigned lower=0;
        const unsigned upper=1; //nele_array.size();

        for (unsigned iMeshSize=lower; iMeshSize<upper; iMeshSize++) {

          const unsigned nele = nele_array[iMeshSize];

          PerceptMesh eMesh(3u);

          build_perturbed_cube(eMesh, nele, true);

          int innerIter = 1001;
          const double smooth_tol = 1.e-3;

          eMesh.setProperty("in_filename","perturbed_"+std::to_string(nele)+"x"+std::to_string(nele)+"x"+std::to_string(nele));
          stk::diag::Timer root_timer(eMesh.getProperty("in_filename"), rootTimerStructured());
          stk::diag::TimeBlock root_block(root_timer);

          {
              SGridBoundarySelector boundarySelector (&eMesh);

              ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> pmmpsi(&eMesh, NULL,&boundarySelector, 0, innerIter, smooth_tol, 1);
              pmmpsi.m_do_animation = 0;
              pmmpsi.run();
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
              
              EXPECT_NEAR(avgx, 4.00001, 5.e-4);
            }
          }
        }
      }

      TEST(DISABLED_isolate, bugs)
      {
//          unsigned p_size = 10000000;
//          release_vs_debug_openmp_fail_reduce r(p_size);
//          unsigned amt = r.reduce();
//          std::cout << amt << std::endl;
//          EXPECT_NEAR(amt, 3*p_size, 0);//passes this
//
//          release_vs_debug_openmp_fail_for f(p_size);
//          bool isGood = f.fill();
//          std::cout << isGood << std::endl;
//          EXPECT_NEAR(isGood,1,0);



//          printf("size of ulli on host is %d",sizeof(unsigned long long int));

          const uint64_t STEP= 10000000;
          const uint64_t LIM = 180000000;//200000000; //GPU runs OUT OF MEMORY AFTER 18E7
          const uint64_t ITERS=10;

          for(uint64_t ihm=STEP;ihm<=LIM;ihm+=STEP)
          {
              /*det_arr_tester_w_data*/function_user datwd(ihm);


              stk::diag::Timer arr_test(/*"arr_test_w_data"*/"function_user"+std::to_string(ihm)+"_"+std::to_string(ITERS), rootTimerStructured());
              stk::diag::TimeBlock view_arr_tb(arr_test);

              for(uint64_t i=0;i<ITERS;i++)
                  datwd.run();

              std::cout << ihm << std::endl;
          }

          for(uint64_t ihm=STEP;ihm<=LIM;ihm+=STEP)
          {
              det_arr_tester dat(ihm);

              stk::diag::Timer arr_test("arr_test"+std::to_string(ihm)+"_"+std::to_string(ITERS), rootTimerStructured());
              stk::diag::TimeBlock view_arr_tb(arr_test);

              for(uint64_t i=0;i<ITERS;i++)
                  dat.run();

              std::cout << ihm << std::endl;
          }

      }
  }//heavy_tests
}//percept

#endif
