// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/SingleTetFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

  namespace percept
  {
    namespace unit_tests
    {

      //=============================================================================
      //=============================================================================
      //=============================================================================

#if !defined(NO_GEOM_SUPPORT)
      /// test stretch_eigens

      TEST(perceptMesh, jacobian_stretch_eigens)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        MPI_Barrier( MPI_COMM_WORLD );

        const unsigned p_size = stk::parallel_machine_size( pm );
        const unsigned p_rank = stk::parallel_machine_rank( pm );
        if (p_size <= 2)
          {
            // create a 12x12 quad mesh with sidesets
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets_on);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("quad fixture",  0);
            //eMesh.save_as("./output_files/quad_fixture.e");

            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

            const stk::mesh::BucketVector & buckets = bulkData.buckets( eMesh.element_rank() );

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (in_surface_selector(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        double eigens[3] = { 1.0, 1.0, 1.0 };

                        JacobianUtil jac;
                        jac.stretch_eigens(eMesh, element, eigens, coordField, cell_topo_data);
                        const double tol=1.e-6;
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[0], 1.0, tol);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[1], 1.0, tol);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[2], 1.0, tol);
                        if (0) std::cout << "P[" << p_rank << ":" << p_size << "] element= " << eMesh.identifier(element)
                                         << " eigens = " << eigens[0] << " " << eigens[1] << " " << eigens[2] << std::endl;
                      }
                  }
              }

          }

        if (p_size <= 2)
          {
            // create a 12x12 quad mesh with sidesets
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, sidesets_on);
            fixture.meta_data.commit();
            fixture.set_bounding_box(0,1,0,2);
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("quad fixture",  0);
            eMesh.save_as("jac-test-1.e");

            double min_max_ave[3]={0,0,0};
            double hmesh = eMesh.hmesh_stretch_eigens(min_max_ave);
            const double tol=1.e-6;
            EXPECT_DOUBLE_EQ_APPROX_TOL(hmesh, 2.0/12.0, tol);
            EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[0], 2.0/12.0, tol);
            EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[1], 2.0/12.0, tol);
            EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[2], 2.0/12.0, tol);

            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

            const stk::mesh::BucketVector & buckets = bulkData.buckets( eMesh.element_rank() );

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (in_surface_selector(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        double eigens[3];

                        JacobianUtil jac;
                        jac.stretch_eigens(eMesh, element, eigens, coordField, cell_topo_data);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[0], 1.0, tol);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[1], 2.0/12.0, tol);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[2], 1.0/12.0, tol);
                        if (0) std::cout << "P[" << p_rank << ":" << p_size << "] element= " << eMesh.identifier(element)
                                         << " eigens = " << eigens[0] << " " << eigens[1] << " " << eigens[2] << std::endl;
                      }
                  }
              }

          }

        /// triangles
        if (1 && p_size <= 2)
          {
            // create a 12x12 quad mesh with sidesets
            const unsigned n = 12;
            //const unsigned nx = n , ny = n , nz = p_size*n ;
            const unsigned nx = n , ny = n;

            bool sidesets_on = true;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, sidesets_on);
            fixture.meta_data.commit();
            fixture.set_bounding_box(0,1,0,2);
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.print_info("tri fixture",  0);
            eMesh.save_as("jac-test-tri.e");

            double min_max_ave[3]={0,0,0};
            double hmesh = eMesh.hmesh_stretch_eigens(min_max_ave);

            const double tol=1.e-6;
            const double sqp = 1./6.*std::sqrt(1./6.*(5.+std::sqrt(13.)));
            const double sqm = 1./6.*std::sqrt(1./6.*(5.-std::sqrt(13.)));
            EXPECT_DOUBLE_EQ_APPROX_TOL(hmesh, sqp, tol);
            EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[0], sqp, tol);
            EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[1], sqp, tol);
            EXPECT_DOUBLE_EQ_APPROX_TOL(min_max_ave[2], sqp, tol);

            stk::mesh::BulkData& bulkData = *eMesh.get_bulk_data();
            CoordinatesFieldType* coordField = eMesh.get_coordinates_field();

            const stk::mesh::BucketVector & buckets = bulkData.buckets( eMesh.element_rank() );

            for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
              {
                //if (in_surface_selector(**k))
                  {
                    stk::mesh::Bucket & bucket = **k ;

                    // in case the cell topology is needed
                    const CellTopologyData * const cell_topo_data = eMesh.get_cell_topology(bucket);
                    shards::CellTopology cell_topo(cell_topo_data);

                    const unsigned num_elements_in_bucket = bucket.size();

                    for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                      {
                        stk::mesh::Entity element = bucket[iElement];
                        double eigens[3];

                        JacobianUtil jac;
                        jac.stretch_eigens(eMesh, element, eigens, coordField, cell_topo_data);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[0], 1.0, tol);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[1], sqp, tol);
                        EXPECT_DOUBLE_EQ_APPROX_TOL(eigens[2], sqm, tol);
                        if (0) std::cout << "P[" << p_rank << ":" << p_size << "] element= " << eMesh.identifier(element)
                                         << " eigens = " << eigens[0] << " " << eigens[1] << " " << eigens[2] << std::endl;
                      }
                  }
              }
          }
      }

#endif


    }
  }
