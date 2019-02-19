// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <math.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <gtest/gtest.h>

#include <stk_util/environment/WallTime.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/PerceptUtils.hpp>

#if !STK_PERCEPT_LITE
#include <percept/structured/StructuredGridRefiner.hpp>
#include <percept/structured/StructuredGridSnap.hpp>
#include <percept/structured/PGeomAssocStructured.hpp>
#endif

#include <PGeom.hpp>
#ifdef HAVE_ACIS
#include <PGeomACIS.hpp>
#endif

// Ioss
#if defined(STK_BUILT_IN_SIERRA)
#include <main/io_info.h>
#endif

namespace percept
{
  namespace regression_tests
  {

#if defined(STK_BUILT_IN_SIERRA)
#if !STK_PERCEPT_LITE


    struct CylinderTestData
    {
        int surface_id;
        double expected_xy_radius;

        CylinderTestData(int surf_id, double radius) : surface_id(surf_id), expected_xy_radius(radius) {};
    };

    void do_refine_cgns_structured_with_geom(const std::string& file,
                                             const std::string& geom_file,
                                             std::vector<CylinderTestData> &cylinders,
                                             int print_level=0,
                                             bool allow_parallel=false)
    {
      stk::diag::Timer     my_timer(file, rootTimerStructured());
      stk::diag::TimeBlock my_timeblock(my_timer);

      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1 && !allow_parallel) return;

      eMesh.open(file, "cgns_structured");
      eMesh.readBulkData();

      // snap the refined mesh to the geometry
      if (geom_file.length() != 0.0)
        {
          std::shared_ptr<PGeom> pgeom;

          if      (geom_file.find(".3dm") != std::string::npos) {
              pgeom = std::shared_ptr<PGeom>(new PGeom());
              pgeom->initialize(OPENNURBS_GEOMETRY_ENGINE);
              pgeom->import_open_nurbs_file(geom_file.c_str());
          }
          // bcarnes: this code does not work for some reason
#ifdef HAVE_ACIS
          else if (geom_file.find(".sat") != std::string::npos) {
              pgeom = std::shared_ptr<PGeom>(new PGeomACIS());
              pgeom->initialize(     ACIS_GEOMETRY_ENGINE);
              pgeom->import_acis_file(      geom_file.c_str());
          }
#endif
          else
              throw std::runtime_error("Could not read geom file = " + geom_file);

          {
              bool use_local_tolerance = false;
              double default_tolerance = 1e-4;

              PGeomAssocStructured assoc_structured;
              assoc_structured.find_mesh_association(pgeom,
                                                     eMesh.get_block_structured_grid(),
                                                     use_local_tolerance,
                                                     default_tolerance,
                                                     print_level);

              StructuredGridRefiner sgr(eMesh.get_block_structured_grid(), print_level);

              size_t pos = file.find(".cgns");
              VERIFY_OP_ON(pos, !=, std::string::npos, "bad filename: need .cgns extension, file= "+file);
              std::string prefix = file.substr(0,pos);
              if (print_level)
                sgr.m_input->print(std::cout, print_level-1);
              sgr.m_input->dump_vtk(prefix);

              sgr.do_refine();

              assoc_structured.make_refined_mesh_association(eMesh.get_block_structured_grid(),
                                                             sgr.m_output);

              StructuredGridSnap str_snap(sgr.m_output);
              str_snap.snap_to_geometry(pgeom, assoc_structured);

              sgr.m_output->dump_vtk(prefix+"_ref");
              if (print_level)
                sgr.m_output->print(std::cout, print_level-1);

              for (CylinderTestData cyl_test : cylinders)
              {
                  assoc_structured.check_cylinder_surface_projection(pgeom, sgr.m_output, cyl_test.surface_id,
                                                                     cyl_test.expected_xy_radius,
                                                                     default_tolerance);
              }
            }
        }
    }

    TEST(regr_PerceptCGNS, test_sgrid_refine_and_snap1)
    {
        std::vector<CylinderTestData> cylinder_tests;
        int cylinder_id = 6;
        double cylinder_radius = 7.0;
        cylinder_tests.push_back(CylinderTestData(cylinder_id, cylinder_radius));

        do_refine_cgns_structured_with_geom("bulge.cgns", "bulge.3dm", cylinder_tests, 0, true);
    }

    TEST(regr_PerceptCGNS, test_sgrid_refine_and_snap2)
    {
        std::vector<CylinderTestData> cylinder_tests;
        int cylinder_id = 4;
        double cylinder_radius = 2.0;
        cylinder_tests.push_back(CylinderTestData(cylinder_id, cylinder_radius));
        cylinder_id = 5;
        cylinder_radius = 5.0;
        cylinder_tests.push_back(CylinderTestData(cylinder_id, cylinder_radius));

        do_refine_cgns_structured_with_geom("half_cyl.cgns", "half_cyl.3dm", cylinder_tests, 0, true);
    }

    TEST(regr_PerceptCGNS, test_sgrid_refine_and_snap3)
    {
        std::vector<CylinderTestData> cylinder_tests;
        do_refine_cgns_structured_with_geom("quarter_sphere.cgns", "quarter_sphere.3dm", cylinder_tests, 1, true);
    }
//    TEST(regr_PerceptCGNS, test_sgrid_refine_and_snap4)
//    {
//      //do_refine_cgns_structured_with_geom("TransitionDuct.cgns", "TransitionDuct.3dm", 0, true);
//      do_refine_cgns_structured_with_geom("TransitionDuct.cgns", "TransitionDuct.sat", 0, true);
//    }

    void do_associate_to_geom(const std::string& file, const std::string& geom_file, int print_level=0, bool allow_parallel=false)
    {
      stk::diag::Timer     my_timer(file, rootTimerStructured());
      stk::diag::TimeBlock my_timeblock(my_timer);

      PerceptMesh eMesh(3u);
      if (!allow_parallel && eMesh.get_parallel_size() > 1) return;

      eMesh.open(file, "cgns_structured");
      eMesh.readBulkData();

      // snap the refined mesh to the geometry
      if (geom_file.length() != 0.0)
        {
          std::shared_ptr<PGeom> pgeom(new PGeom());

          if (pgeom->initialize(OPENNURBS_GEOMETRY_ENGINE))
            {
              pgeom->import_open_nurbs_file(geom_file.c_str());

              bool use_local_tolerance = false;
              double default_tolerance = 1e-4;

              PGeomAssocStructured assoc_structured;
              assoc_structured.find_mesh_association(pgeom,
                                                     eMesh.get_block_structured_grid(),
                                                     use_local_tolerance,
                                                     default_tolerance,
                                                     print_level);
             }
        }
    }

    TEST(regr_PerceptCGNS, test_sgrid_associate_geometry)
    {
#if !defined(KOKKOS_ENABLE_CUDA)
        do_associate_to_geom("bulge.cgns", "bulge.3dm", 0, false);
#endif
    }

#endif // STK_PERCEPT_LITE
#endif // STK_BUILT_IN_SIERRA

  } // namespace regression_tests
} // namespace percept

