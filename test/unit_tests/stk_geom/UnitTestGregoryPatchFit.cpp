// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <gtest/gtest.h>

#include <percept/mesh/geometry/stk_geom/3D/FitGregoryPatches.hpp>
#include <percept/mesh/geometry/stk_geom/3D/EvaluateGregoryPatch.hpp>

#include <percept/PerceptMesh.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/TriQuadSurfaceMesh3D.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/SingleTetFixture.hpp>
#include <percept/fixtures/TetWedgeFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/PyramidFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

namespace percept
{
  namespace unit_tests
  {
    //=============================================================================
    //=============================================================================
    //=============================================================================

    TEST(unit_stk_geom_3d, test_1)
    {
      const bool debug_print = false;
      //SplineFit::s_debug_print = debug_print;
      EXCEPTWATCH;
      stk::ParallelMachine pm = MPI_COMM_WORLD ;
      MPI_Barrier( MPI_COMM_WORLD );

      const unsigned p_rank = stk::parallel_machine_rank( pm );
      const unsigned p_size = stk::parallel_machine_size( pm );

      if (p_size <= 1)
        {
          unsigned n = 2;
          std::cout << "P["<<p_rank<<"] " << "tmp srk test_1 n = " << n << std::endl;
          std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
          PerceptMesh eMesh(3);
          eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));

          FitGregoryPatches::SurfaceSets surfaceSets;
          FitGregoryPatches::AngleMap angleMap;
          FitGregoryPatches fitter(eMesh, surfaceSets, angleMap, 15.0);
          fitter.register_or_set_fields();

          eMesh.commit();
          eMesh.setProperty("FitGregoryPatches::noEdgeFitting", "true");
          std::vector<stk::mesh::PartVector> currentParts;
          fitter.getCurrentParts(currentParts);
          if (debug_print)
            std::cout << "P["<<p_rank<<"] " << " currentParts.size= " << currentParts.size() << std::endl;
          for (unsigned iSurfaceSet = 0; iSurfaceSet < currentParts.size(); ++iSurfaceSet)
            {
              // normals
              fitter.getNormals(&currentParts[iSurfaceSet]);
            }

          eMesh.save_as("stk_geom_3d_test1.e");
        }
    }

    TEST(unit_stk_geom_3d, test_2)
    {
      const bool debug_print = false;
      //SplineFit::s_debug_print = debug_print;
      EXCEPTWATCH;
      stk::ParallelMachine pm = MPI_COMM_WORLD ;
      MPI_Barrier( MPI_COMM_WORLD );

      const unsigned p_rank = stk::parallel_machine_rank( pm );
      const unsigned p_size = stk::parallel_machine_size( pm );

      if (p_size <= 1)
        {
          unsigned n = 2;
          if (debug_print) std::cout << "P["<<p_rank<<"] " << "tmp srk test_1 n = " << n << std::endl;
          std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
          PerceptMesh eMesh(3);
          eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));

          FitGregoryPatches fitter(eMesh);
          fitter.register_or_set_fields();

          eMesh.commit();
          eMesh.setProperty("FitGregoryPatches::noEdgeFitting", "true");

          fitter.computeControlPoints();
          int faces[] = {11,14,15,21};

          double expected_data[][20][3]={

            {{0,0,0},{0.111111,-0.0555556,-0.0555556},{0.333333,0,0},{0.5,0,0},{-0.0555556,-0.0555556,0.111111},
             {0.037037,-0.157407,0.00925927},{0.00925927,-0.157407,0.037037},{0.305556,-0.0833333,0.0833333},
             {0.333333,-0.0833333,0.0833333},{0.5,-0.0833333,0.0833333},{0,0,0.333333},{0.0833333,-0.0833333,0.333333},
             {0.0833333,-0.0833333,0.305556},{0.333333,2.1027e-19,0.333333},{0.333333,2.1027e-19,0.333333},
             {0.5,0,0.333333},{0,0,0.5},{0.0833333,-0.0833333,0.5},{0.333333,0,0.5},{0.5,0,0.5}},

            {{0,0,0},{-0.0555556,-0.0555556,0.111111},{0,0,0.333333},{0,0,0.5},{-0.0555556,0.111111,-0.0555556},
             {-0.157407,0.00925927,0.037037},{-0.157407,0.037037,0.00925927},{-0.0833333,0.0833333,0.305556},
             {-0.0833333,0.0833333,0.333333},{-0.0833333,0.0833333,0.5},{0,0.333333,0},{-0.0833333,0.333333,0.0833333},
             {-0.0833333,0.305556,0.0833333},{6.30808e-19,0.333333,0.333333},{6.30808e-19,0.333333,0.333333},
             {0,0.333333,0.5},{0,0.5,0},{-0.0833333,0.5,0.0833333},{0,0.5,0.333333},{0,0.5,0.5}},

            {{0,0,0},{-0.0555556,0.111111,-0.0555556},{0,0.333333,0},{0,0.5,0},{0.111111,-0.0555556,-0.0555556},
             {0.00925927,0.037037,-0.157407},{0.037037,0.00925927,-0.157407},{0.0833333,0.305556,-0.0833333},
             {0.0833333,0.333333,-0.0833333},{0.0833333,0.5,-0.0833333},{0.333333,0,0},{0.333333,0.0833333,-0.0833333},
             {0.305556,0.0833333,-0.0833333},{0.333333,0.333333,2.1027e-19},{0.333333,0.333333,2.1027e-19},
             {0.333333,0.5,0},{0.5,0,0},{0.5,0.0833333,-0.0833333},{0.5,0.333333,0},{0.5,0.5,0}},

            {{0.5,0,0},{0.666667,0,0},{0.888889,-0.0555556,-0.0555556},{1,0,0},{0.5,-0.0833333,0.0833333},
             {0.694444,-0.0833333,0.0833333},{0.666667,-0.0833333,0.0833333},{0.962963,-0.157407,0.00925927},
             {0.990741,-0.157407,0.037037},{1.05556,-0.0555556,0.111111},{0.5,0,0.333333},{0.666667,2.1027e-19,0.333333},
             {0.666667,2.1027e-19,0.333333},{0.916667,-0.0833333,0.333333},{0.916667,-0.0833333,0.305556},
             {1,0,0.333333},{0.5,0,0.5},{0.666667,0,0.5},{0.916667,-0.0833333,0.5},{1,0,0.5}}

          };

          int nfaces = sizeof(faces)/sizeof(int);
          for (unsigned stage=0; stage < 2; ++stage)
            {
              if (stage == 1)
                {
                  fitter.m_reverseAll = !fitter.m_reverseAll;
                  fitter.computeControlPoints();
                }
              std::cout << "face data for fitter.m_reverseAll= " << fitter.m_reverseAll << std::endl;
              for (int jf=0; jf < nfaces; ++jf)
                {
                  stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.side_rank(),faces[jf]);
                  if (eMesh.is_valid(face_i))
                    {
                      //eMesh.print_entity(std::cout, face_i);
                      // to generate the data above, use these settings:
                      //bool convert=false, printHeader=false;
                      bool convert=true, printHeader=true;
                      std::cout << fitter.printForMathematica(face_i, convert, printHeader, 80) << std::endl;
                    }
                }
              for (int jf=0; jf < nfaces; ++jf)
                {
                  stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.side_rank(),faces[jf]);
                  if (eMesh.is_valid(face_i))
                    {
                      double *Cp = stk::mesh::field_data( *eMesh.m_gregory_control_points_field , face_i);
                      for (unsigned k=0; k < 20; ++k)
                        {
                          for (unsigned jc=0; jc < 3; ++jc)
                            {
                              EXPECT_NEAR(expected_data[jf][k][jc], Cp[k+jc*FitGregoryPatches::MaxControlPoints()], 1.e-5);
                            }
                        }
                    }
                }
            }
          eMesh.save_as("stk_geom_3d_test2.e");
          PerceptMesh::create_refined_mesh(eMesh, "ref_test_2.e", 10);
        }
    }

    TEST(unit_stk_geom_3d, test_3)
    {
      stk::ParallelMachine pm = MPI_COMM_WORLD ;

      unsigned nNodes = 5, nTriFaces = 1, nQuadFaces = 1;
      TriQuadSurfaceMesh3D::Point coords[] = {
        {{0,0.0,0.0}},
        {{1.0,0,0.25}},
        {{0,1,0.25}},
        {{1.5,0.5,0}},
        {{1,1.5,0.0}}
      };
      // note: 1-based
      TriQuadSurfaceMesh3D::TriIds tris[] { TriQuadSurfaceMesh3D::TriIds {2,3,1} };
      TriQuadSurfaceMesh3D::QuadIds quads[] { TriQuadSurfaceMesh3D::QuadIds {2,4,5,3} };
      //const unsigned p_rank = stk::parallel_machine_rank( pm );
      const unsigned p_size = stk::parallel_machine_size( pm );
      if (p_size <= 1)
        {
          TriQuadSurfaceMesh3D mesh(pm, false);
          bool isCommitted = false;
          percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
          FitGregoryPatches fitter(eMesh);
          fitter.register_or_set_fields();
          eMesh.commit();
          eMesh.setProperty("FitGregoryPatches::noEdgeFitting", "true");

          mesh.populate(nNodes, nTriFaces, nQuadFaces, coords, tris, quads);

          eMesh.save_as("tri_quad_0.e");

          stk::mesh::Entity face_1 = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);
          std::cout << "is_valid = " << eMesh.is_valid(face_1) <<  std::endl;
          std::cout << "num_faces= " << eMesh.get_bulk_data()->num_faces(face_1) << std::endl;

          fitter.computeControlPoints();
          int faces[] = {1,2};

          double expected_data[][20][3]={

            {{1,0,0.25},{0.749921,0.249724,0.257008},{0.5,0.5,0.25},{0.250079,0.750276,0.242992},
             {0,1,0.25},{0.749325,-0.00234722,0.247069},{0.502475,0.237955,0.250456},{0.499048,0.246686,0.271597},
             {0.252488,0.488721,0.239483},{0.249259,0.497423,0.252909},{-0.000595196,0.747929,0.240061},
             {0.49955,-0.00156481,0.164713},{0.250026,0.250092,0.122664},{0.249974,0.249908,0.127336},
             {-0.000396797,0.498619,0.160041},{0.25,0,0.0625},{0,0.25,0.0625},{0,0,0}},

            {{1,0,0.25},{1.16563,0.163056,0.258301},{1.32951,0.330468,0.0661414},{1.5,0.5,0},{0.666561,0.332965,0.259344},
             {0.886507,0.495189,0.258799},{0.889831,0.490079,0.263339},{1.10471,0.661768,0.0577948},
             {1.10602,0.662846,0.0604107},{1.32951,0.830468,-0.017192},{0.333439,0.667035,0.240656},
             {0.611396,0.830849,0.259686},{0.612033,0.824062,0.25533},{0.895292,1.0049,0.108872},
             {0.893983,1.00382,0.106256},{1.17049,1.16953,0.017192},{0,1,0.25},{0.332274,1.16298,0.260185},
             {0.670487,1.3362,0.100525},{1,1.5,0}}

          };

          int nfaces = sizeof(faces)/sizeof(int);
          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  // to generate the data above, use these settings:
                  bool convert=false, printHeader=false;
                  //bool convert=true, printHeader=true;
                  std::cout << fitter.printForMathematica(face_i, convert, printHeader, 80) << std::endl;
                }
            }
          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  //eMesh.print_entity(std::cout, face_i);
                  double *Cp = (eMesh.entity_rank(face_i) == eMesh.side_rank()
                                ? stk::mesh::field_data( *eMesh.m_gregory_control_points_field, face_i)
                                : stk::mesh::field_data( *eMesh.m_gregory_control_points_field_shell, face_i));

                  unsigned nn = eMesh.get_bulk_data()->num_nodes(face_i);
                  for (unsigned k=0; k < (nn == 3 ? 18 : 20); ++k)
                    {
                      for (unsigned jc=0; jc < 3; ++jc)
                        {
                          EXPECT_NEAR(expected_data[jf][k][jc], Cp[k+jc*FitGregoryPatches::MaxControlPoints()], 1.e-5);
                        }
                    }

                }
            }
          eMesh.save_as("stk_geom_3d_test3.e");

          // now, reproduce the FarinHansford test case

          // set normals
          TriQuadSurfaceMesh3D::Point normals[] = { {{-0.639602, -0.639602, 0.426401}}, {{0.229416, -0.688247, 0.688247}},
                                 {{-0.688247, 0.229416, 0.688247}}, {{0.801784, -0.267261, 0.534522}},
                                 {{0.267261, 0.801784, 0.534522}} };

          eMesh.save_as("farin-hansford-normals-0.e");
          const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.node_rank() );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  double *ndata = stk::mesh::field_data( *eMesh.m_node_normals , node );
                  stk::mesh::EntityId id = eMesh.identifier(node);
                  Math::normalize_3d(normals[id-1].c_array());
                  for (unsigned jc=0; jc < 3; ++jc)
                    {
                      ndata[jc] = normals[id-1][jc];
                    }
                }
            }
          eMesh.save_as("farin-hansford-normals.e");

          // re-fit the data
          // Note: these are flattened and sorted to make comparison easy with results from FarinHansford
          double expected_tri_data[] =
            {-0.1028708133971292, -0.1028708133971292, -0.08522727272727271,
             -0.08522727272727271, -0.06907894736842109, -0.06907894736842109, 0,
             0, 0, 0., 0., 0.004175964256766342, 0.02422248803827741,
             0.02422248803827741, 0.026913875598086154, 0.026913875598086154,
             0.07602853386275349, 0.09210526315789466, 0.09210526315789466,
             0.09708931419457738, 0.09708931419457738, 0.1193181818181818,
             0.1193181818181818, 0.1647727272727273, 0.1647727272727273,
             0.20893141945773525, 0.20893141945773525, 0.25, 0.25,
             0.2565789473684211, 0.2565789473684211, 0.291267942583732,
             0.291267942583732, 0.40380781499202556, 0.40380781499202556,
             0.4078947368421053, 0.4078947368421053, 0.4298245614035087,
             0.4298245614035087, 0.4299774108048654, 0.4306842766789487,
             0.4585326953748006, 0.4585326953748006, 0.4605263157894738,
             0.5049483426508528, 0.5373803827751196, 0.5373803827751196,
             0.5771382342907696, 0.7730263157894737, 0.7730263157894737,
             0.8026315789473684, 0.8026315789473684, 1., 1.};

          double expected_quad_data[] =
            {0, 0, 0, 0., 0.07456140350877187, 0.10714285714285716,
             0.11904761904761905, 0.11904761904761908, 0.12280701754385956,
             0.12280701754385956, 0.1622807017543859, 0.1785714285714286,
             0.2477025898078528, 0.25, 0.25, 0.25000000000000006,
             0.2566833751044278, 0.25877192982456143, 0.2611668832396433,
             0.27380952380952384, 0.3195488721804509, 0.3214285714285714,
             0.32811194653299924, 0.3377192982456141, 0.3577955482936576,
             0.43880534670008364, 0.4605263157894738, 0.4605263157894738,
             0.46458149838700574, 0.49160415373477423, 0.5, 0.5177527151211363,
             0.5449039264828737, 0.6468253968253967, 0.7142857142857143,
             0.7368421052631579, 0.7368421052631579, 0.7738095238095237,
             0.8642784429202794, 0.9060150375939848, 0.9232327341017901,
             0.992063492063492, 0.9966583124477861, 1, 1, 1., 1.0258980785296572,
             1.1973684210526314, 1.2236842105263157, 1.226190476190476,
             1.26984126984127, 1.2852965747702592, 1.289264828738513,
             1.3293650793650793, 1.3452380952380953, 1.369047619047619,
             1.4761904761904763, 1.5, 1.5, 1.5119047619047619};

          std::cout << "re-fit to FarinHansford data" << std::endl;
          bool doGetNormals = false;
          fitter.computeControlPoints(doGetNormals);
          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  std::cout << fitter.printForMathematica(face_i, false) << std::endl;
                }
            }
          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  //eMesh.print_entity(std::cout, face_i);
                  double *Cp = (eMesh.entity_rank(face_i) == eMesh.side_rank()
                                ? stk::mesh::field_data( *eMesh.m_gregory_control_points_field, face_i)
                                : stk::mesh::field_data( *eMesh.m_gregory_control_points_field_shell, face_i));

                  unsigned nn = eMesh.get_bulk_data()->num_nodes(face_i);
                  unsigned npts = (nn == 3 ? 18 : 20);
                  std::vector<double> Cpv;
                  for (unsigned k=0; k < npts; ++k)
                    {
                      for (unsigned jc=0; jc < 3; ++jc)
                        {
                          Cpv.push_back( Cp[k+jc*FitGregoryPatches::MaxControlPoints()]);
                        }
                    }
                  std::sort(Cpv.begin(), Cpv.end());
                  for (unsigned k=0; k < Cpv.size(); ++k)
                    {
                      double exp_data = (nn == 3 ? expected_tri_data[k] : expected_quad_data[k]);
                      if (std::fabs(exp_data - Cpv[k]) > 1.e-5)
                        std::cout << "k = " << k << " nn= " << nn << " expected_data= " << exp_data << " Cp= " << Cpv[k] << " diff= " << (exp_data-Cpv[k]) << std::endl;
                      //EXPECT_NEAR(expected_data, Cpv[k], 1.e-6);
                    }
                }
            }
          PerceptMesh::create_refined_mesh(eMesh, "ref_test_3.e", 10);
        }
    }

    static void do_test(const std::string& prefix,
                        unsigned nNodes , unsigned nTriFaces , unsigned nQuadFaces,
                        TriQuadSurfaceMesh3D::Point coords[],
                        TriQuadSurfaceMesh3D::TriIds tris[],
                        TriQuadSurfaceMesh3D::QuadIds quads[],
                        double expected_data[][20][3],
                        int faces[],
                        int map[],
                        TriQuadSurfaceMesh3D::Point normals[],
                        double expected_tri_data[][18*3],
                        double expected_quad_data[][20*3],
                        int num_divisions = 10
                        )
    {
      stk::ParallelMachine pm = MPI_COMM_WORLD ;

      //const unsigned p_rank = stk::parallel_machine_rank( pm );
      const unsigned p_size = stk::parallel_machine_size( pm );
      if (p_size <= 1)
        {
          TriQuadSurfaceMesh3D mesh(pm, false);
          bool isCommitted = false;
          percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
          FitGregoryPatches fitter(eMesh);
          fitter.register_or_set_fields();
          eMesh.commit();
          eMesh.setProperty("FitGregoryPatches::noEdgeFitting", "true");

          mesh.populate(nNodes, nTriFaces, nQuadFaces, coords, tris, quads);

          eMesh.save_as(prefix+"_0.e");

          fitter.computeControlPoints();

          int nfaces = nTriFaces + nQuadFaces;
          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  // to generate the data above, use these settings:
                  bool convert=false, printHeader=false;
                  //bool convert=true, printHeader=true;
                  std::cout << fitter.printForMathematica(face_i, convert, printHeader, 80) << std::endl;
                }
            }
          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  //eMesh.print_entity(std::cout, face_i);
                  double *Cp = (eMesh.entity_rank(face_i) == eMesh.side_rank()
                                ? stk::mesh::field_data( *eMesh.m_gregory_control_points_field, face_i)
                                : stk::mesh::field_data( *eMesh.m_gregory_control_points_field_shell, face_i));

                  unsigned nn = eMesh.get_bulk_data()->num_nodes(face_i);
                  for (unsigned k=0; k < (nn == 3 ? 18 : 20); ++k)
                    {
                      for (unsigned jc=0; jc < 3; ++jc)
                        {
                          EXPECT_NEAR(expected_data[jf][k][jc], Cp[k+jc*FitGregoryPatches::MaxControlPoints()], 1.e-5);
                        }
                    }

                }
            }
          eMesh.save_as(prefix+"_1.e");

          // now, reproduce the FarinHansford test case
          const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.node_rank() );
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();
              for (unsigned iNode = 0; iNode < num_nodes_in_bucket; iNode++)
                {
                  stk::mesh::Entity node = bucket[iNode];
                  double *ndata = stk::mesh::field_data( *eMesh.m_node_normals , node );
                  stk::mesh::EntityId id = eMesh.identifier(node);
                  Math::normalize_3d(normals[id-1].c_array());
                  for (unsigned jc=0; jc < 3; ++jc)
                    {
                      ndata[jc] = normals[id-1][jc];
                    }
                }
            }

          // re-fit the data

          std::cout << "re-fit to FarinHansford data" << std::endl;
          bool doGetNormals = false;
          fitter.computeControlPoints(doGetNormals);
          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  std::cout << fitter.printForMathematica(face_i, false) << std::endl;
                }
            }

          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  //eMesh.print_entity(std::cout, face_i);
                  double *Cp = (eMesh.entity_rank(face_i) == eMesh.side_rank()
                                ? stk::mesh::field_data( *eMesh.m_gregory_control_points_field, face_i)
                                : stk::mesh::field_data( *eMesh.m_gregory_control_points_field_shell, face_i));

                  unsigned nn = eMesh.get_bulk_data()->num_nodes(face_i);
                  unsigned npts = (nn == 3 ? 18 : 20);
                  std::vector<double> Cpv;
                  for (unsigned k=0; k < npts; ++k)
                    {
                      for (unsigned jc=0; jc < 3; ++jc)
                        {
                          Cpv.push_back( Cp[k+jc*FitGregoryPatches::MaxControlPoints()]);
                        }
                    }
                  std::sort(Cpv.begin(), Cpv.end());

                  for (unsigned k=0; k < Cpv.size(); ++k)
                    {
                      double exp_data = (nn == 3 ? expected_tri_data[map[jf]][k] : expected_quad_data[map[jf]][k]);
                      if (std::fabs(exp_data - Cpv[k]) > 1.e-5)
                        std::cout << "k = " << k << " jf= " << jf << " nn= " << nn << " expected_data= " << exp_data << " Cp= " << Cpv[k] << " diff= " << (exp_data-Cpv[k]) << std::endl;
                      EXPECT_NEAR(exp_data, Cpv[k], 1.e-5);
                    }
                }
            }
          PerceptMesh::create_refined_mesh(eMesh, "out1.e", num_divisions);
        }
    }

    TEST(unit_stk_geom_3d, test_4)
    {
      unsigned nNodes = 7, nTriFaces = 3, nQuadFaces = 2;
      TriQuadSurfaceMesh3D::Point coords[] = {
        {{0., 2., 0.}}, {{-1.90211, 0.618034, 0.}}, {{-1.17557, -1.61803, 0.}},
        {{1.17557, -1.61803, 0.}}, {{1.90211, 0.618034, 0.}}, {{0, 0, 1.76556}}, {{0., -1.80278, 0.651388}}
      };
      // note: 1-based
      TriQuadSurfaceMesh3D::TriIds tris[] = { {{3,4,7}}, {{7,4,6}}, {{4,5,6}} };
      TriQuadSurfaceMesh3D::QuadIds quads[] = { {{1,2,6,5}}, {{2,3,7,6}} };
      //const unsigned p_rank = stk::parallel_machine_rank( pm );

      int faces[] = {1,2,3,4,5};  // always list tris first
      int map[] = {0,1,2,0,1};

      double expected_data[][20][3]={

        {{-1.17557,-1.61803,0},{-0.695365,-1.83313,0.0734427},{0.0295007,-1.90226,0.146865},
         {0.739616,-1.82928,0.146855},{1.17557,-1.61803,0},{-0.898218,-1.69729,0.174139},{-0.442149,-1.9817,0.274305},
         {-0.41047,-1.93777,0.274751},{0.50324,-1.97418,0.37443},{0.370577,-1.99016,0.386574},
         {0.900307,-1.69014,0.180866},{-0.604015,-1.79198,0.37373},{0.000502539,-2.07067,0.571443},
         {-0.0707369,-2.03235,0.511574},{0.592966,-1.8105,0.394064},{-0.301697,-1.84588,0.549303},
         {0.283034,-1.88081,0.573077},{0,-1.80278,0.651388}},

        {{0,-1.80278,0.651388},{0.283034,-1.88081,0.573077},{0.592966,-1.8105,0.394064},{0.900307,-1.69014,0.180866},
         {1.17557,-1.61803,0},{-0.0155206,-1.52964,1.05076},{0.299878,-1.64699,0.999855},{0.267687,-1.67243,1.03417},
         {0.658923,-1.51812,0.781796},{0.685504,-1.38917,0.734224},{0.992692,-1.36798,0.548766},
         {0.0170775,-1.03649,1.4325},{0.378912,-0.936364,1.33048},{0.355802,-1.02073,1.49855},
         {0.702817,-0.937006,1.16898},{0.0411369,-0.475784,1.70223},{0.355426,-0.442036,1.64609},
         {0,0,1.76556}},

        {{1.17557,-1.61803,0},{1.51119,-1.27326,0.148942},{1.70528,-0.603355,0.228942},{1.81614,0.118231,0.194472},
         {1.90211,0.618034,0},{0.992692,-1.36798,0.548766},{1.32293,-1.05534,0.671499},{1.28401,-0.93635,0.792514},
         {1.37009,-0.0924225,0.656351},{1.27957,-0.156777,0.774432},{1.36746,0.42693,0.3212},
         {0.702817,-0.937006,1.16898},{0.949716,-0.346862,1.83676},{0.940562,-0.425111,1.39213},
         {0.956355,0.257347,1.0366},{0.355426,-0.442036,1.64609},{0.542603,0.113599,1.67509},
         {0,0,1.76556}},

        {{0,2,0},{-0.634037,1.7411,0.228549},{-1.32587,1.11071,0.132386},{-1.90211,0.618034,0},
         {0.634037,1.7411,0.228549},{0.0305102,1.48293,0.507201},{-0.0262779,1.48357,0.52367},
         {-0.630811,0.853262,0.46114},{-0.467625,0.670502,1.08613},{-1.17654,0.361305,0.378857},
         {1.33133,1.11785,0.128598},{0.672244,0.691489,0.838195},{0.644742,0.861676,0.490291},
         {0.21869,0.357063,1.99904},{0.274317,0.308052,2.90582},{-0.500483,0.124558,1.87575},
         {1.90211,0.618034,0},{1.18924,0.363228,0.428266},{0.723471,0.151466,1.64493},{0,0,1.76556}},

        {{-1.90211,0.618034,0},{-1.7599,-0.0719288,0.228987},{-1.64618,-1.32941,0.155947},
         {-1.17557,-1.61803,0},{-1.17654,0.361305,0.378857},{-1.15292,-0.278307,0.558953},
         {-1.00367,-0.360553,0.900971},{-1.15779,-1.48544,0.437023},{-1.13494,-1.42677,0.520134},
         {-0.805767,-1.72371,0.232185},{-0.500483,0.124558,1.87575},{-0.434912,-0.619079,1.65502},
         {-0.33993,-0.590736,2.70052},{-0.466182,-1.45285,1.09477},{-0.579509,-1.54185,0.933058},
         {-0.402263,-1.86025,0.515275},{0,0,1.76556},{0.0548492,-0.634379,1.68112},{-0.0206941,-1.43859,1.18389},
         {0,-1.80278,0.651388}}

      };

      TriQuadSurfaceMesh3D::Point normals[] = {
        {{0., 0.992278, 0.124035}}, {{-0.943712, 0.306631, 0.124035}},
        {{-0.583246, -0.80277, 0.124035}}, {{0.583246, -0.80277, 0.124035}},
        {{0.943712, 0.306631, 0.124035}}, {{0., 0., 1.}}, {{0., -0.894427, 0.447214}}
      };

      double expected_tri_data[][18*3] =
        {{-2.122789499065, -2.1227894990649996, -2.038847897678856, -2.032970216150994, -1.9849784773652446,
          -1.9789908316858547, -1.9750574703353603, -1.8932423552114073, -1.8932423552114073, -1.8586773371694894,
          -1.8586773371694894, -1.8395466273324446, -1.8395466273324446, -1.8027756377319946, -1.75584501044965,
          -1.7558450104496497, -1.618033988749895, -1.6180339887498947, -1.1755705045849463, -0.9482477803616179,
          -0.7877358345896461, -0.6321651869077453, -0.5826833698501034, -0.5192858052898774, -0.29389262614623657,
          -0.041404378013531314, -1.1102230246251565e-16, 0., 0., 0., 0.0022527969474931686, 0.04252203111538782,
          0.04252203111538783, 0.05669604148718377, 0.17700388993715876, 0.1770038899371588, 0.23369993142434253,
          0.23369993142434256, 0.23476051218377852, 0.2645055821619979, 0.2938926261462365, 0.36916090347444486,
          0.36916090347444486, 0.41108248173295614, 0.45694896937180474, 0.46513107466165005, 0.5395844199910078,
          0.5395844199910078, 0.5826833698501032, 0.6321651869077451, 0.6513878188659973, 0.7877358345896459,
          0.9482477803616176, 1.175570504584946},
         {-1.8586773371694894, -1.8395466273324446, -1.819692181389536,
          -1.8027756377319946, -1.75584501044965, -1.7088772605779516, -1.6873372464870886, -1.618033988749895,
          -1.601219194024531, -1.5678608846621667, -1.566146758187251, -1.0674794626830206, -1.0452405897747779,
          -1.0343466973584206, -1.0117809653603458, -0.4506939094329986, -0.4045084971874737, 0., 0., 0., 0., 0., 0., 0.,
          0.17700388993715876, 0.2938926261462365, 0.2938926261462365, 0.3042743854350562, 0.30475193425802377,
          0.34398462139802366, 0.36916090347444486, 0.40765360804454037, 0.4961389383568337, 0.5395844199910078,
          0.6321651869077451, 0.6513878188659973, 0.721892419408445, 0.7340613408381139, 0.7594117404687813,
          0.7899999772193829, 0.8104239326046958, 0.9447553317146793, 0.9482477803616176, 1.0260818493717214,
          1.0545007062809248, 1.139117610703172, 1.175570504584946, 1.2135415107752077, 1.3937759901378093,
          1.4772180529136023, 1.621331608805313, 1.7655644370746373, 1.7655644370746373, 1.7655644370746373},
         {-1.618033988749895, -1.5678608846621667, -1.33422536083646, -1.2125618889536383, -1.1839159203685716,
          -1.0452405897747779, -0.6133920829743678, -0.4119496880703892, -0.4045084971874737, -0.27573988160444396, 0., 0.,
          0., 0., 0.016786406652205127, 0.04252203111538782, 0.042522031115387854, 0.05669604148718378, 0.14169214241306152,
          0.15450849718747364, 0.16413723637490832, 0.2938926261462365, 0.3992463788729792, 0.4755282581475768,
          0.4961389383568337, 0.49613893835683376, 0.5528349798440175, 0.5528349798440175, 0.5988695683094688,
          0.6180339887498946, 0.664541461532378, 0.7222449788100959, 0.7594117404687813, 0.7819916390125279,
          0.887084220074153, 1.139117610703172, 1.175570504584946, 1.2135415107752077, 1.2135415107752079, 1.2287540075342027,
          1.4479505180228482, 1.539750517614548, 1.5571567188834592, 1.6403895206321133, 1.7655644370746373,
          1.7655644370746373, 1.7655644370746373, 1.843131011301304, 1.8878267156038482, 1.9021130325903073,
          2.006797650428053, 2.0076331733040544, 2.0398289954048057, 2.044004238816126}};

      double expected_quad_data[][20*3] =
        {{-1.902113032590307, -1.823470337538303, -1.699444472696417, -1.493569855436121, -1.4356701168183803,
          -0.6340376775301023, -0.6340376775301023, -0.18513166082603272, -0.02455431808193187, 0, 0, 0., 0., 0., 0.,
          0.05669604148718381, 0.056696041487183824, 0.05669604148718383, 0.05669604148718384, 0.11738246763047006,
          0.18513166082603283, 0.20601132958329818, 0.2060113295832983, 0.31499959730389915, 0.31499959730389915,
          0.44806366077446985, 0.4789786274124744, 0.5166071116334306, 0.5166071116334306, 0.5924814281626603,
          0.5924814281626607, 0.6180339887498946, 0.6180339887498949, 0.6340376775301024, 0.6340376775301024,
          0.6615185844757784, 0.6615185844757784, 0.8000638721319938, 0.8937013384728383, 1.0444748465480451,
          1.056204937157409, 1.1994522730464223, 1.199452273046423, 1.2188496484998779, 1.2188496484998783,
          1.4356701168183803, 1.4544551993882389, 1.6994444726964173, 1.7655644370746373, 1.7655644370746373,
          1.7655644370746373, 1.823470337538303, 1.9021130325903073, 1.9796708044944256, 1.9796708044944256,
          1.9929129948141022, 1.9929129948141022, 2., 2.0209807362995877, 2.233356763407489},
         {-2.0913013075580658, -2.010704954806625, -1.9378405223130282, -1.902113032590307, -1.8773112369819875,
          -1.823470337538303, -1.8027756377319946, -1.8017820176829014, -1.6843521236496306, -1.6180339887498947,
          -1.588663917966741, -1.5732838601403372, -1.5340337127887098, -1.4731274877837723, -1.4560732482115597,
          -1.3706386910163972, -1.2396224848653148, -1.2220783752454434, -1.1755705045849463, -0.8724735389538418,
          -0.6340376775301023, -0.6009252125773314, -0.5782204672481287, -0.5471984030316877, -0.5432890603311795,
          -0.48047957929618756, -0.4755969026901859, -0.4659939066525691, -0.39185683486164874, -0.06544606445257828,
          -0.00048055544813625783, 0, 0, 0., 0., 0., 0., 0., 0.012838318916580077, 0.056696041487183775,
          0.056696041487183844, 0.2060113295832983, 0.2360051865828784, 0.43453902736769545, 0.4519637521726792,
          0.5023166203660113, 0.5763768266653289, 0.5924814281626607, 0.6180339887498949, 0.6513878188659973,
          0.6615185844757784, 0.8171280848025331, 0.8940387274471406, 1.0425347800105367, 1.1888716687525673,
          1.6502071393191313, 1.7655644370746373, 1.7655644370746373, 1.7655644370746373, 2.1567835097371857}};

      do_test("sphere_cap", nNodes ,  nTriFaces ,  nQuadFaces, coords,
               tris, quads, expected_data, faces, map, normals, expected_tri_data, expected_quad_data);
    }

    TEST(unit_stk_geom_3d, test_5)
    {
      unsigned nNodes = 18, nTriFaces = 8, nQuadFaces = 8;

      // always list tris first
      int faces[] = {
        1,2,3,4,5,6,7,8,
        9,10,11,12,13,14,15,16
      };
      int map[] = {
        0,1,2,3,4,5,6,7,
        0,1,2,3,4,5,6,7
      };

      double r2 = 1.7;

      TriQuadSurfaceMesh3D::Point coords[] = {
        {{-1,1,r2}},
        {{-1,-1,2}},
        {{1,-1,2}},
        {{1,1,r2}},

        {{-r2,1,1}},
        {{-2,-1,1}},
        {{-2,-1,-1}},
        {{-r2,1,-1}},

        {{-1,-1,-2}},
        {{1,-1,-2}},
        {{1,1,-r2}},
        {{-1,1,-r2}},

        {{2,-1,-1}},
        {{r2,1,-1}},
        {{r2,1,1}},
        {{2,-1,1}},

        {{0, -4, 0}},
        {{0, 4, 0}}
      };

      //cc={0,0,0};
      TriQuadSurfaceMesh3D::Point normals[18]= {{{}}};
      for (unsigned ii = 0; ii < 18; ++ii)
        {
          for (unsigned jc=0; jc < 3; ++jc)
            normals[ii][jc] = coords[ii][jc];
          Math::normalize_3d(normals[ii].c_array());
        }

      TriQuadSurfaceMesh3D::QuadIds quads[] = {
        {{1,2,3,4}},
        {{3,16,15,4}},
        {{13,14,15,16}},
        {{11,14,13,10}},
        {{11,10,9,12}},
        {{12,9,7,8}},
        {{7,6,5,8}},
        {{6,2,1,5}}
      };

      TriQuadSurfaceMesh3D::TriIds tris[] = {
        {{17,3,2}},
        {{17,16,3}},
        {{17,13,16}},
        {{17,10,13}},
        {{17,9,10}},
        {{17,7,9}},
        {{17,6,7}},
        {{17,2,6}}
      };

      double expected_data[][20][3]={

        {{0,-4,0},{0.25,-4,0.5},{0.586018,-3.059,1.20834},{0.879027,-1.83849,1.81251},{1,-1,2},
         {-0.25,-4,0.5},{0.0809156,-3.40229,1.19598},{-0.0809156,-3.40229,1.19598},{0.427616,-2.04002,1.9302},
         {0.385177,-1.90141,2.0347},{0.568149,-1.04674,2.16506},{-0.586018,-3.059,1.20834},
         {-0.385177,-1.90141,2.0347},{-0.427616,-2.04002,1.9302},{0,-1.06232,2.22008},{-0.879027,-1.83849,1.81251},
         {-0.568149,-1.04674,2.16506},{-1,-1,2}},

        {{0,-4,0},{0.5,-4,0.25},{1.20834,-3.059,0.586018},{1.81251,-1.83849,0.879027},{2,-1,1},
         {0.25,-4,0.5},{0.935039,-3.27203,0.826397},{0.826397,-3.27203,0.935039},{1.65288,-1.90618,1.16634},
         {1.65543,-1.88221,1.25424},{1.86737,-1.03323,1.29846},{0.586018,-3.059,1.20834},{1.25424,-1.88221,1.65543},
         {1.16634,-1.90618,1.65288},{1.61055,-1.04431,1.61055},{0.879027,-1.83849,1.81251},
         {1.29846,-1.03323,1.86737},{1,-1,2}},

        {{0,-4,0},{0.5,-4,-0.25},{1.20834,-3.059,-0.586018},{1.81251,-1.83849,-0.879027},{2,-1,-1},
         {0.5,-4,0.25},{1.19598,-3.40229,-0.0809156},{1.19598,-3.40229,0.0809156},{1.9302,-2.04002,-0.427616},
         {2.0347,-1.90141,-0.385177},{2.16506,-1.04674,-0.568149},{1.20834,-3.059,0.586018},
         {2.0347,-1.90141,0.385177},{1.9302,-2.04002,0.427616},{2.22008,-1.06232,0},{1.81251,-1.83849,0.879027},
         {2.16506,-1.04674,0.568149},{2,-1,1}},

        {{0,-4,0},{0.25,-4,-0.5},{0.586018,-3.059,-1.20834},{0.879027,-1.83849,-1.81251},{1,-1,-2},
         {0.5,-4,-0.25},{0.826397,-3.27203,-0.935039},{0.935039,-3.27203,-0.826397},{1.16634,-1.90618,-1.65288},
         {1.25424,-1.88221,-1.65543},{1.29846,-1.03323,-1.86737},{1.20834,-3.059,-0.586018},
         {1.65543,-1.88221,-1.25424},{1.65288,-1.90618,-1.16634},{1.61055,-1.04431,-1.61055},
         {1.81251,-1.83849,-0.879027},{1.86737,-1.03323,-1.29846},{2,-1,-1}},

        {{0,-4,0},{-0.25,-4,-0.5},{-0.586018,-3.059,-1.20834},{-0.879027,-1.83849,-1.81251},
         {-1,-1,-2},{0.25,-4,-0.5},{-0.0809156,-3.40229,-1.19598},{0.0809156,-3.40229,-1.19598},
         {-0.427616,-2.04002,-1.9302},{-0.385177,-1.90141,-2.0347},{-0.568149,-1.04674,-2.16506},
         {0.586018,-3.059,-1.20834},{0.385177,-1.90141,-2.0347},{0.427616,-2.04002,-1.9302},
         {0,-1.06232,-2.22008},{0.879027,-1.83849,-1.81251},{0.568149,-1.04674,-2.16506},{1,-1,-2}},

        {{0,-4,0},{-0.5,-4,-0.25},{-1.20834,-3.059,-0.586018},{-1.81251,-1.83849,-0.879027},
         {-2,-1,-1},{-0.25,-4,-0.5},{-0.935039,-3.27203,-0.826397},{-0.826397,-3.27203,-0.935039},
         {-1.65288,-1.90618,-1.16634},{-1.65543,-1.88221,-1.25424},{-1.86737,-1.03323,-1.29846},
         {-0.586018,-3.059,-1.20834},{-1.25424,-1.88221,-1.65543},{-1.16634,-1.90618,-1.65288},
         {-1.61055,-1.04431,-1.61055},{-0.879027,-1.83849,-1.81251},{-1.29846,-1.03323,-1.86737},
         {-1,-1,-2}},

        {{0,-4,0},{-0.5,-4,0.25},{-1.20834,-3.059,0.586018},{-1.81251,-1.83849,0.879027},{-2,-1,1},
         {-0.5,-4,-0.25},{-1.19598,-3.40229,0.0809156},{-1.19598,-3.40229,-0.0809156},{-1.9302,-2.04002,0.427616},
         {-2.0347,-1.90141,0.385177},{-2.16506,-1.04674,0.568149},{-1.20834,-3.059,-0.586018},
         {-2.0347,-1.90141,-0.385177},{-1.9302,-2.04002,-0.427616},{-2.22008,-1.06232,0},{-1.81251,-1.83849,-0.879027},
         {-2.16506,-1.04674,-0.568149},{-2,-1,-1}},

        {{0,-4,0},{-0.25,-4,0.5},{-0.586018,-3.059,1.20834},{-0.879027,-1.83849,1.81251},{-1,-1,2},
         {-0.5,-4,0.25},{-0.826397,-3.27203,0.935039},{-0.935039,-3.27203,0.826397},{-1.16634,-1.90618,1.65288},
         {-1.25424,-1.88221,1.65543},{-1.29846,-1.03323,1.86737},{-1.20834,-3.059,0.586018},
         {-1.65543,-1.88221,1.25424},{-1.65288,-1.90618,1.16634},{-1.61055,-1.04431,1.61055},
         {-1.81251,-1.83849,0.879027},{-1.86737,-1.03323,1.29846},{-2,-1,1}},

        {{-1,1,1.7},{-1,0.333333,1.8},{-1.09533,-0.398717,2.1309},{-1,-1,2},{-0.429657,1.03477,1.93179},
         {-0.424585,0.335589,2.02683},{-0.429657,0.368102,2.03179},{-0.529626,-0.436874,2.38216},
         {-0.457796,-0.461869,2.35392},{-0.424199,-1.06232,2.22008},{0.429657,1.03477,1.93179},
         {0.424585,0.335589,2.02683},{0.429657,0.368102,2.03179},{0.529626,-0.436874,2.38216},
         {0.457796,-0.461869,2.35392},{0.424199,-1.06232,2.22008},{1,1,1.7},{1,0.333333,1.8},
         {1.09533,-0.398717,2.1309},{1,-1,2}},

        {{1,-1,2},{1.39794,-1.04431,1.82315},{1.82315,-1.04431,1.39794},{2,-1,1},{1.09533,-0.398717,2.1309},
         {1.50435,-0.442201,1.93885},{1.46284,-0.432641,2.0061},{1.93885,-0.442201,1.50435},
         {2.0061,-0.432641,1.46284},{2.1309,-0.398717,1.09533},{1,0.333333,1.8},{1.31408,0.350448,1.64743},
         {1.31691,0.327934,1.66607},{1.64743,0.350448,1.31408},{1.66607,0.327934,1.31691},
         {1.8,0.333333,1},{1,1,1.7},{1.28075,1.01711,1.58077},{1.58077,1.01711,1.28075},{1.7,1,1}},

        {{2,-1,-1},{2.1309,-0.398717,-1.09533},{1.8,0.333333,-1},{1.7,1,-1},{2.22008,-1.06232,-0.424199},
         {2.38216,-0.436874,-0.529626},{2.35392,-0.461869,-0.457796},{2.02683,0.335589,-0.424585},
         {2.03179,0.368102,-0.429657},{1.93179,1.03477,-0.429657},{2.22008,-1.06232,0.424199},
         {2.38216,-0.436874,0.529626},{2.35392,-0.461869,0.457796},{2.02683,0.335589,0.424585},
         {2.03179,0.368102,0.429657},{1.93179,1.03477,0.429657},{2,-1,1},{2.1309,-0.398717,1.09533},
         {1.8,0.333333,1},{1.7,1,1}},

        {{1,1,-1.7},{1.28075,1.01711,-1.58077},{1.58077,1.01711,-1.28075},{1.7,1,-1},{1,0.333333,-1.8},
         {1.31408,0.350448,-1.64743},{1.31691,0.327934,-1.66607},{1.64743,0.350448,-1.31408},
         {1.66607,0.327934,-1.31691},{1.8,0.333333,-1},{1.09533,-0.398717,-2.1309},{1.50435,-0.442201,-1.93885},
         {1.46284,-0.432641,-2.0061},{1.93885,-0.442201,-1.50435},{2.0061,-0.432641,-1.46284},
         {2.1309,-0.398717,-1.09533},{1,-1,-2},{1.39794,-1.04431,-1.82315},{1.82315,-1.04431,-1.39794},
         {2,-1,-1}},

        {{1,1,-1.7},{1,0.333333,-1.8},{1.09533,-0.398717,-2.1309},{1,-1,-2},{0.429657,1.03477,-1.93179},
         {0.424585,0.335589,-2.02683},{0.429657,0.368102,-2.03179},{0.529626,-0.436874,-2.38216},
         {0.457796,-0.461869,-2.35392},{0.424199,-1.06232,-2.22008},{-0.429657,1.03477,-1.93179},
         {-0.424585,0.335589,-2.02683},{-0.429657,0.368102,-2.03179},{-0.529626,-0.436874,-2.38216},
         {-0.457796,-0.461869,-2.35392},{-0.424199,-1.06232,-2.22008},{-1,1,-1.7},{-1,0.333333,-1.8},
         {-1.09533,-0.398717,-2.1309},{-1,-1,-2}},

        {{-1,1,-1.7},{-1,0.333333,-1.8},{-1.09533,-0.398717,-2.1309},{-1,-1,-2},{-1.28075,1.01711,-1.58077},
         {-1.31691,0.327934,-1.66607},{-1.31408,0.350448,-1.64743},{-1.46284,-0.432641,-2.0061},
         {-1.50435,-0.442201,-1.93885},{-1.39794,-1.04431,-1.82315},{-1.58077,1.01711,-1.28075},
         {-1.66607,0.327934,-1.31691},{-1.64743,0.350448,-1.31408},{-2.0061,-0.432641,-1.46284},
         {-1.93885,-0.442201,-1.50435},{-1.82315,-1.04431,-1.39794},{-1.7,1,-1},{-1.8,0.333333,-1},
         {-2.1309,-0.398717,-1.09533},{-2,-1,-1}},

        {{-2,-1,-1},{-2.22008,-1.06232,-0.424199},{-2.22008,-1.06232,0.424199},{-2,-1,1},{-2.1309,-0.398717,-1.09533},
         {-2.35392,-0.461869,-0.457796},{-2.38216,-0.436874,-0.529626},{-2.35392,-0.461869,0.457796},
         {-2.38216,-0.436874,0.529626},{-2.1309,-0.398717,1.09533},{-1.8,0.333333,-1},{-2.03179,0.368102,-0.429657},
         {-2.02683,0.335589,-0.424585},{-2.03179,0.368102,0.429657},{-2.02683,0.335589,0.424585},
         {-1.8,0.333333,1},{-1.7,1,-1},{-1.93179,1.03477,-0.429657},{-1.93179,1.03477,0.429657},
         {-1.7,1,1}},

        {{-2,-1,1},{-1.82315,-1.04431,1.39794},{-1.39794,-1.04431,1.82315},{-1,-1,2},{-2.1309,-0.398717,1.09533},
         {-1.93885,-0.442201,1.50435},{-2.0061,-0.432641,1.46284},{-1.50435,-0.442201,1.93885},
         {-1.46284,-0.432641,2.0061},{-1.09533,-0.398717,2.1309},{-1.8,0.333333,1},{-1.64743,0.350448,1.31408},
         {-1.66607,0.327934,1.31691},{-1.31408,0.350448,1.64743},{-1.31691,0.327934,1.66607},
         {-1,0.333333,1.8},{-1.7,1,1},{-1.58077,1.01711,1.28075},{-1.28075,1.01711,1.58077},
         {-1,1,1.7}}

      };


      double expected_tri_data[][18*3] =
        {{-4., -4., -4., -3.444444444444444, -3.4444444444444438, -3.0555555555555554, -3.0555555555555554,
          -2.0233333333333334, -2.0233333333333325, -1.9494085552648646, -1.9494085552648646, -1.8333333333333333,
          -1.8333333333333333, -1.1111111111111112, -1.0833333333333333, -1.0833333333333333, -1., -1., -1.,
          -0.8333333333333334, -0.5833333333333334, -0.5555555555555556, -0.3963350551195186, -0.39333333333333337, -0.25,
          -0.05944444444444444, 0., 0., 0., 0.05944444444444434, 0.25, 0.39333333333333337, 0.3963350551195189, 0.5, 0.5,
          0.5555555555555556, 0.5833333333333334, 0.8333333333333334, 1., 1.1111111111111112, 1.1111111111111112,
          1.1188888888888886, 1.1188888888888888, 1.6666666666666667, 1.6666666666666667, 1.7866666666666668,
          1.786666666666667, 1.8988171105297285, 1.898817110529729, 2., 2., 2.1666666666666665, 2.1666666666666665,
          2.2222222222222223},
         {-4., -4., -4., -3.2361111111111107, -3.2361111111111107, -3.0555555555555554,
          -3.0555555555555554, -1.9477777777777774, -1.9477777777777767, -1.8855405762622262, -1.885540576262226,
          -1.8333333333333333, -1.8333333333333333, -1.0555555555555556, -1.0416666666666667, -1.0416666666666667, -1., -1.,
          0., 0., 0.25, 0.25, 0.5, 0.5, 0.5555555555555556, 0.5555555555555556, 0.8061111111111111, 0.8061111111111111,
          0.8333333333333334, 0.8333333333333334, 0.8622222222222223, 0.8622222222222226, 1., 1., 1.1077777777777778,
          1.1077777777777782, 1.1111111111111112, 1.1111111111111112, 1.164155432196669, 1.1641554321966694,
          1.2916666666666665, 1.2916666666666665, 1.465555555555556, 1.4655555555555566, 1.4924662965900077,
          1.4924662965900082, 1.5833333333333333, 1.5833333333333333, 1.6666666666666667, 1.6666666666666667,
          1.8333333333333335, 1.8333333333333335, 2., 2.},
         {-4., -4., -4., -3.4444444444444446, -3.444444444444444,
          -3.0555555555555554, -3.0555555555555554, -2.0233333333333334, -2.023333333333333, -1.9494085552648641,
          -1.949408555264864, -1.8333333333333333, -1.8333333333333333, -1.1111111111111112, -1.0833333333333333,
          -1.0833333333333333, -1., -1., -1., -0.8333333333333334, -0.5833333333333334, -0.5555555555555556,
          -0.3963350551195189, -0.3933333333333335, -0.25, -0.05944444444444425, 0., 0., 0., 0.059444444444444376, 0.25,
          0.39333333333333337, 0.3963350551195189, 0.5, 0.5, 0.5555555555555556, 0.5833333333333334, 0.8333333333333334, 1.,
          1.1111111111111112, 1.1111111111111112, 1.1188888888888886, 1.1188888888888888, 1.6666666666666667,
          1.6666666666666667, 1.7866666666666662, 1.7866666666666668, 1.898817110529729, 1.898817110529729, 2., 2.,
          2.1666666666666665, 2.1666666666666665, 2.2222222222222223},
         {-4., -4., -4., -3.236111111111111,
          -3.2361111111111107, -3.0555555555555554, -3.0555555555555554, -2., -1.9477777777777783, -1.9477777777777774,
          -1.8855405762622262, -1.8855405762622257, -1.8333333333333335, -1.8333333333333333, -1.8333333333333333,
          -1.6666666666666667, -1.5833333333333333, -1.4924662965900082, -1.465555555555556, -1.2916666666666665,
          -1.1641554321966692, -1.1111111111111112, -1.1077777777777778, -1.0555555555555556, -1.0416666666666667,
          -1.0416666666666667, -1., -1., -1., -0.862222222222222, -0.8333333333333334, -0.8061111111111112,
          -0.5555555555555556, -0.5, -0.25, 0., 0., 0.25, 0.5, 0.5555555555555556, 0.806111111111111, 0.8333333333333334,
          0.8622222222222223, 1., 1.1077777777777778, 1.1111111111111112, 1.1641554321966692, 1.2916666666666665,
          1.465555555555556, 1.4924662965900082, 1.5833333333333333, 1.6666666666666667, 1.8333333333333335, 2.},
         {-4., -4.,
          -4., -3.444444444444444, -3.4444444444444438, -3.0555555555555554, -3.0555555555555554, -2.2222222222222223,
          -2.1666666666666665, -2.1666666666666665, -2.0233333333333334, -2.0233333333333325, -2., -2., -1.9494085552648646,
          -1.9494085552648646, -1.898817110529729, -1.8988171105297285, -1.8333333333333333, -1.8333333333333333,
          -1.786666666666667, -1.7866666666666668, -1.6666666666666667, -1.6666666666666667, -1.1188888888888888,
          -1.1188888888888886, -1.1111111111111112, -1.1111111111111112, -1.1111111111111112, -1.0833333333333333,
          -1.0833333333333333, -1., -1., -1., -0.8333333333333334, -0.5833333333333334, -0.5555555555555556, -0.5, -0.5,
          -0.3963350551195189, -0.39333333333333337, -0.25, -0.05944444444444434, 0., 0., 0., 0.05944444444444444, 0.25,
          0.39333333333333337, 0.3963350551195186, 0.5555555555555556, 0.5833333333333334, 0.8333333333333334, 1.},
         {-4.,
          -4., -4., -3.2361111111111107, -3.2361111111111107, -3.0555555555555554, -3.0555555555555554, -2., -2.,
          -1.9477777777777774, -1.9477777777777767, -1.8855405762622262, -1.885540576262226, -1.8333333333333335,
          -1.8333333333333335, -1.8333333333333333, -1.8333333333333333, -1.6666666666666667, -1.6666666666666667,
          -1.5833333333333333, -1.5833333333333333, -1.4924662965900082, -1.4924662965900077, -1.4655555555555566,
          -1.465555555555556, -1.2916666666666665, -1.2916666666666665, -1.1641554321966694, -1.164155432196669,
          -1.1111111111111112, -1.1111111111111112, -1.1077777777777782, -1.1077777777777778, -1.0555555555555556,
          -1.0416666666666667, -1.0416666666666667, -1., -1., -1., -1., -0.8622222222222226, -0.8622222222222223,
          -0.8333333333333334, -0.8333333333333334, -0.8061111111111111, -0.8061111111111111, -0.5555555555555556,
          -0.5555555555555556, -0.5, -0.5, -0.25, -0.25, 0., 0.},
         {-4., -4., -4., -3.4444444444444446, -3.444444444444444,
          -3.0555555555555554, -3.0555555555555554, -2.2222222222222223, -2.1666666666666665, -2.1666666666666665,
          -2.0233333333333334, -2.023333333333333, -2., -2., -1.9494085552648641, -1.949408555264864, -1.898817110529729,
          -1.898817110529729, -1.8333333333333333, -1.8333333333333333, -1.7866666666666668, -1.7866666666666662,
          -1.6666666666666667, -1.6666666666666667, -1.1188888888888888, -1.1188888888888886, -1.1111111111111112,
          -1.1111111111111112, -1.1111111111111112, -1.0833333333333333, -1.0833333333333333, -1., -1., -1.,
          -0.8333333333333334, -0.5833333333333334, -0.5555555555555556, -0.5, -0.5, -0.3963350551195189,
          -0.39333333333333337, -0.25, -0.059444444444444376, 0., 0., 0., 0.05944444444444425, 0.25, 0.3933333333333335,
          0.3963350551195189, 0.5555555555555556, 0.5833333333333334, 0.8333333333333334, 1.},
         {-4., -4., -4.,
          -3.236111111111111, -3.2361111111111107, -3.0555555555555554, -3.0555555555555554, -2., -1.9477777777777783,
          -1.9477777777777774, -1.8855405762622262, -1.8855405762622257, -1.8333333333333335, -1.8333333333333333,
          -1.8333333333333333, -1.6666666666666667, -1.5833333333333333, -1.4924662965900082, -1.465555555555556,
          -1.2916666666666665, -1.1641554321966692, -1.1111111111111112, -1.1077777777777778, -1.0555555555555556,
          -1.0416666666666667, -1.0416666666666667, -1., -1., -1., -0.8622222222222223, -0.8333333333333334,
          -0.806111111111111, -0.5555555555555556, -0.5, -0.25, 0., 0., 0.25, 0.5, 0.5555555555555556, 0.8061111111111112,
          0.8333333333333334, 0.862222222222222, 1., 1.1077777777777778, 1.1111111111111112, 1.1641554321966692,
          1.2916666666666665, 1.465555555555556, 1.4924662965900082, 1.5833333333333333, 1.6666666666666667,
          1.8333333333333335, 2.}};

      double expected_quad_data[][20*3] =
        {{-1.1444444444444446, -1.1111111111111112, -1.1111111111111112, -1.1015678254942058, -1., -1., -1., -1.,
          -0.6184230472606548, -0.5969387983274075, -0.5969387983274075, -0.5715590608060064, -0.5526892433882773,
          -0.5526892433882766, -0.5046674567503712, -0.5035219268348102, -0.47777777777777786, -0.47777777777777786,
          -0.4696659850034083, -0.4444444444444444, 0.4349011588275392, 0.4349011588275392, 0.4444444444444444,
          0.4696659850034083, 0.5035219268348102, 0.5046674567503711, 0.5102230702533861, 0.510223070253387,
          0.5712338104976142, 0.5712338104976142, 0.5715590608060063, 0.6184230472606549, 1., 1., 1., 1., 1.1015678254942058,
          1.136332651670075, 1.136332651670075, 1.1444444444444446, 1.7, 1.7, 1.9317655078391274, 1.9317655078391274,
          1.9726653033401498, 1.9726653033401498, 2., 2., 2.188888888888889, 2.188888888888889, 2.204430811179277,
          2.204430811179277, 2.220108000413337, 2.2201080004133407, 2.2222222222222223, 2.2222222222222223,
          2.4272109299881484, 2.4272109299881484, 2.4665216287357734, 2.466521628735777},
         {-1.0555555555555556,
          -1.0555555555555556, -1., -1., -0.5494253647660937, -0.5494253647660936, -0.5280465239228126, -0.5280465239228124,
          -0.47777777777777786, -0.47777777777777786, 0.4349011588275392, 0.4349011588275392, 0.4606580463534165,
          0.4606580463534167, 0.4683026584867076, 0.4683026584867076, 1., 1., 1., 1., 1., 1., 1.0334014996591683,
          1.0334014996591683, 1.1015678254942058, 1.1015678254942058, 1.1444444444444446, 1.1444444444444446,
          1.2667348329925017, 1.2667348329925017, 1.3888888888888888, 1.3888888888888888, 1.4245724133787394,
          1.4245724133787396, 1.4253351511020222, 1.4253351511020222, 1.5205615986518382, 1.5205615986518388,
          1.5234492160872528, 1.5234492160872528, 1.5441830410902573, 1.544183041090258, 1.7, 1.7, 1.739082026812088,
          1.739082026812088, 1.777777777777778, 1.777777777777778, 1.8155669371217904, 1.8155669371217915, 1.93995653067818,
          1.93995653067818, 1.9726653033401498, 1.9726653033401498, 2., 2., 2.0499378478800403, 2.049937847880041,
          2.188888888888889, 2.188888888888889},
         {-1.1444444444444446, -1.1111111111111112, -1.1111111111111112,
          -1.1015678254942058, -1., -1., -1., -1., -0.618423047260656, -0.596938798327408, -0.5969387983274078,
          -0.5715590608060053, -0.5526892433882772, -0.5526892433882771, -0.5046674567503708, -0.5035219268348102,
          -0.47777777777777786, -0.47777777777777786, -0.4696659850034083, -0.4444444444444444, 0.4349011588275392,
          0.4349011588275392, 0.4444444444444444, 0.4696659850034083, 0.5035219268348102, 0.5046674567503707,
          0.5102230702533866, 0.5102230702533868, 0.5712338104976142, 0.5712338104976142, 0.5715590608060058,
          0.6184230472606552, 1., 1., 1., 1., 1.1015678254942058, 1.136332651670075, 1.136332651670075, 1.1444444444444446,
          1.7, 1.7, 1.9317655078391274, 1.9317655078391274, 1.9726653033401498, 1.9726653033401498, 2., 2.,
          2.188888888888889, 2.188888888888889, 2.204430811179277, 2.204430811179277, 2.220108000413338, 2.2201080004133393,
          2.2222222222222223, 2.2222222222222223, 2.427210929988148, 2.4272109299881484, 2.4665216287357747,
          2.466521628735776},
         {-2.188888888888889, -2.0499378478800394, -2., -1.9726653033401498, -1.93995653067818,
          -1.8155669371217924, -1.777777777777778, -1.739082026812088, -1.7, -1.5441830410902573, -1.5234492160872528,
          -1.520561598651839, -1.4253351511020222, -1.424572413378739, -1.3888888888888888, -1.2667348329925017,
          -1.1444444444444446, -1.1015678254942058, -1.0555555555555556, -1.0555555555555556, -1., -1., -1., -1.,
          -0.5494253647660937, -0.5494253647660936, -0.5280465239228126, -0.5280465239228124, -0.47777777777777786,
          -0.47777777777777786, 0.4349011588275392, 0.4349011588275392, 0.46065804635341595, 0.46065804635341684,
          0.4683026584867076, 0.4683026584867076, 1., 1., 1., 1., 1.0334014996591683, 1.0334014996591683, 1.1015678254942058,
          1.1444444444444446, 1.2667348329925017, 1.3888888888888888, 1.4245724133787383, 1.4253351511020222,
          1.5205615986518402, 1.5234492160872528, 1.544183041090258, 1.7, 1.739082026812088, 1.777777777777778,
          1.8155669371217924, 1.9399565306781799, 1.9726653033401498, 2., 2.0499378478800394, 2.188888888888889},

         {-2.466521628735777, -2.4665216287357734, -2.4272109299881484, -2.4272109299881484, -2.2222222222222223,
          -2.2222222222222223, -2.2201080004133407, -2.220108000413337, -2.204430811179277, -2.204430811179277,
          -2.188888888888889, -2.188888888888889, -2., -2., -1.9726653033401498, -1.9726653033401498, -1.9317655078391274,
          -1.9317655078391274, -1.7, -1.7, -1.1444444444444446, -1.1111111111111112, -1.1111111111111112,
          -1.1015678254942058, -1., -1., -1., -1., -0.6184230472606549, -0.5969387983274075, -0.5969387983274075,
          -0.5715590608060063, -0.5526892433882773, -0.5526892433882766, -0.5046674567503711, -0.5035219268348102,
          -0.47777777777777786, -0.47777777777777786, -0.4696659850034083, -0.4444444444444444, 0.4349011588275392,
          0.4349011588275392, 0.4444444444444444, 0.4696659850034083, 0.5035219268348102, 0.5046674567503712,
          0.5102230702533861, 0.510223070253387, 0.5712338104976142, 0.5712338104976142, 0.5715590608060064,
          0.6184230472606548, 1., 1., 1., 1., 1.1015678254942058, 1.136332651670075, 1.136332651670075, 1.1444444444444446},

         {-2.188888888888889, -2.188888888888889, -2.049937847880041, -2.0499378478800403, -2., -2., -1.9726653033401498,
          -1.9726653033401498, -1.93995653067818, -1.93995653067818, -1.8155669371217915, -1.8155669371217904,
          -1.777777777777778, -1.777777777777778, -1.739082026812088, -1.739082026812088, -1.7, -1.7, -1.544183041090258,
          -1.5441830410902573, -1.5234492160872528, -1.5234492160872528, -1.5205615986518388, -1.5205615986518382,
          -1.4253351511020222, -1.4253351511020222, -1.4245724133787396, -1.4245724133787394, -1.3888888888888888,
          -1.3888888888888888, -1.2667348329925017, -1.2667348329925017, -1.1444444444444446, -1.1444444444444446,
          -1.1015678254942058, -1.1015678254942058, -1.0555555555555556, -1.0555555555555556, -1., -1., -1., -1., -1., -1.,
          -0.5494253647660937, -0.5494253647660936, -0.5280465239228126, -0.5280465239228124, -0.47777777777777786,
          -0.47777777777777786, 0.4349011588275392, 0.4349011588275392, 0.4606580463534165, 0.4606580463534167,
          0.4683026584867076, 0.4683026584867076, 1., 1., 1.0334014996591683, 1.0334014996591683},
         {-2.466521628735776,
          -2.4665216287357747, -2.4272109299881484, -2.427210929988148, -2.2222222222222223, -2.2222222222222223,
          -2.2201080004133393, -2.220108000413338, -2.204430811179277, -2.204430811179277, -2.188888888888889,
          -2.188888888888889, -2., -2., -1.9726653033401498, -1.9726653033401498, -1.9317655078391274, -1.9317655078391274,
          -1.7, -1.7, -1.1444444444444446, -1.1111111111111112, -1.1111111111111112, -1.1015678254942058, -1., -1., -1., -1.,
          -0.6184230472606552, -0.596938798327408, -0.5969387983274078, -0.5715590608060058, -0.5526892433882772,
          -0.5526892433882771, -0.5046674567503707, -0.5035219268348102, -0.47777777777777786, -0.47777777777777786,
          -0.4696659850034083, -0.4444444444444444, 0.4349011588275392, 0.4349011588275392, 0.4444444444444444,
          0.4696659850034083, 0.5035219268348102, 0.5046674567503708, 0.5102230702533866, 0.5102230702533868,
          0.5712338104976142, 0.5712338104976142, 0.5715590608060053, 0.618423047260656, 1., 1., 1., 1., 1.1015678254942058,
          1.136332651670075, 1.136332651670075, 1.1444444444444446},
         {-2.188888888888889, -2.0499378478800394, -2.,
          -1.9726653033401498, -1.9399565306781799, -1.8155669371217924, -1.777777777777778, -1.739082026812088, -1.7,
          -1.544183041090258, -1.5234492160872528, -1.5205615986518402, -1.4253351511020222, -1.4245724133787383,
          -1.3888888888888888, -1.2667348329925017, -1.1444444444444446, -1.1015678254942058, -1.0555555555555556,
          -1.0555555555555556, -1., -1., -1., -1., -0.5494253647660937, -0.5494253647660936, -0.5280465239228126,
          -0.5280465239228124, -0.47777777777777786, -0.47777777777777786, 0.4349011588275392, 0.4349011588275392,
          0.46065804635341595, 0.46065804635341684, 0.4683026584867076, 0.4683026584867076, 1., 1., 1., 1.,
          1.0334014996591683, 1.0334014996591683, 1.1015678254942058, 1.1444444444444446, 1.2667348329925017,
          1.3888888888888888, 1.424572413378739, 1.4253351511020222, 1.520561598651839, 1.5234492160872528,
          1.5441830410902573, 1.7, 1.739082026812088, 1.777777777777778, 1.8155669371217924, 1.93995653067818,
          1.9726653033401498, 2., 2.0499378478800394, 2.188888888888889}};

      do_test("cone", nNodes ,  nTriFaces ,  nQuadFaces, coords,
               tris, quads, expected_data, faces, map, normals, expected_tri_data, expected_quad_data);
    }


    TEST(unit_stk_geom_3d, test_6_hex_fit_closest_point)
    {
      const bool debug_print = false;
      //SplineFit::s_debug_print = debug_print;
      EXCEPTWATCH;
      stk::ParallelMachine pm = MPI_COMM_WORLD ;
      MPI_Barrier( MPI_COMM_WORLD );

      const unsigned p_rank = stk::parallel_machine_rank( pm );
      const unsigned p_size = stk::parallel_machine_size( pm );

      if (p_size <= 1)
        {
          unsigned n = 2;
          if (debug_print) std::cout << "P["<<p_rank<<"] " << "tmp srk test_6 n = " << n << std::endl;
          std::string gmesh_spec = toString(n)+"x"+toString(n)+"x"+toString(n)+std::string("|bbox:0,0,0,1,1,1|sideset:xXyYzZ");
          PerceptMesh eMesh(3);
          eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));

          FitGregoryPatches fitter(eMesh);
          fitter.register_or_set_fields();

          eMesh.commit();
          eMesh.setProperty("FitGregoryPatches::noEdgeFitting", "true");

          fitter.computeControlPoints();

          bool debug = false;
          EvaluateGregoryPatch eval(eMesh, debug);

          PerceptMesh::create_refined_mesh(eMesh, "ref_test_6.e", 10);

          int faces[] = {11,14,15,21};

          int nfaces = sizeof(faces)/sizeof(int);

          for (int jf=0; jf < nfaces; ++jf)
            {
              stk::mesh::Entity face_i = eMesh.get_bulk_data()->get_entity(eMesh.side_rank(),faces[jf]);
              if (eMesh.is_valid(face_i))
                {
                  //eMesh.print_entity(std::cout, face_i);
                  // to generate the data above, use these settings:
                  //bool convert=false, printHeader=false;
                  bool convert=true, printHeader=true;
                  std::cout << fitter.printForMathematica(face_i, convert, printHeader, 80) << std::endl;

                  if (1)
                    {
                      //double uv[2] = {.2001,.2001};
                      double xyz[3] = {0.5,0.4,0.3};
                      double c_xyz[3], f_uv[2];
                      bool bad = eval.findClosestPoint(xyz, face_i, c_xyz, f_uv);
                      EXPECT_FALSE(bad);
                      if (debug)
                        std::cout << "bad= " << bad << " err= " << eval.error_message()
                                << " f_uv= " << f_uv[0] << ", " << f_uv[1]
                                << " xyz= " << xyz[0] << ", " << xyz[1] << ", " << xyz[2]
                                << " c_xyz= " << c_xyz[0] << ", " << c_xyz[1] << ", " << c_xyz[2]
                                << std::endl;
                    }

                  if (1)
                    {
                      double ndiv = 3;
                      for (double u=0.0; u <= 1.0; u += 1./ndiv)
                        {
                          for (double v=0.0; v <= 1.0; v += 1./ndiv)
                            {
                              double uv[2] = {u,v};
                              double xyz[3], c_xyz[3], f_uv[2];
                              eval.evaluateGregoryPatch(uv, face_i, xyz);
                              bool bad = eval.findClosestPoint(xyz, face_i, c_xyz, f_uv);
                              EXPECT_FALSE(bad);
                              EXPECT_NEAR(xyz[0], c_xyz[0], 1.e-6);
                              EXPECT_NEAR(xyz[1], c_xyz[1], 1.e-6);
                              EXPECT_NEAR(xyz[2], c_xyz[2], 1.e-6);
                              if (debug)
                                std::cout << "bad= " << bad << " err= " << eval.error_message()
                                          << " uv= " << uv[0] << ", " << uv[1]
                                          << " f_uv= " << f_uv[0] << ", " << f_uv[1]
                                          << " xyz= " << xyz[0] << ", " << xyz[1] << ", " << xyz[2]
                                          << " c_xyz= " << c_xyz[0] << ", " << c_xyz[1] << ", " << c_xyz[2]
                                          << std::endl;
                            }
                        }
                    }
                }
            }
          eMesh.save_as("stk_geom_3d_test_6.e");
        }
    }

  }
}
