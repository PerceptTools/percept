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

#include <adapt/AdaptedMeshVerifier.hpp>

#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/PerceptUtils.hpp>

#include <percept/norm/Norm.hpp>
#include <percept/RebalanceMesh.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>

#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <percept/RunEnvironment.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/PyramidFixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>
#include <percept/fixtures/SingleTetFixture.hpp>
#include <percept/util/Loops.hpp>

#include <adapt/IEdgeBasedAdapterPredicate.hpp>
#include <adapt/ElementRefinePredicate.hpp>
#include <adapt/PredicateBasedEdgeAdapter.hpp>

#include <adapt/FindValidCentroid.hpp>

#include <adapt/PredicateBasedElementAdapter.hpp>
#include <adapt/QualityElementAdapter.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/TEA_SpecialWedgeRefinement.hpp>
#include <adapt/HangingNodeAdapter.hpp>

#include <adapt/PredicateTemplateAdapter.hpp>

#include <percept/function/ElementOp.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MeshUtils.hpp>

#include <regression_tests/RegressionTestLocalRefiner.hpp>
#include <regression_tests/BlockRefiner.hpp>

#include <adapt/markers/MarkerUsingErrIndFraction.hpp>
#if defined(STK_BUILT_IN_SIERRA)
#include <percept/mesh/geometry/volume/sierra_only/FiniteVolumeMesh.hpp>
#endif

//TEMP
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>
#include <adapt/UniformRefinerPattern_Tri3_Tri3_4_sierra.hpp>
#include <adapt/UniformRefinerPattern_def.hpp>

  namespace percept
  {
    namespace regression_tests
    {

#include "RegressionTestFileLoc.hpp"

//       const std::string input_files_loc="./input_files/";
//       const std::string output_files_loc="./output_files/";

      static std::string post_fix(int i) {
        return "np" + toString(i);
      }

#define EXTRA_PRINT 0

      void normalize(double input_normal[3], double normal[3])
      {
        double sum = std::sqrt(input_normal[0]*input_normal[0]+
                               input_normal[1]*input_normal[1]+
                               input_normal[2]*input_normal[2]);
        normal[0] = input_normal[0] / sum;
        normal[1] = input_normal[1] / sum;
        normal[2] = input_normal[2] / sum;
      }

      void normalize(double input_output_normal[3])
      {
        normalize(input_output_normal, input_output_normal);
      }

      double distance(double c0[3], double c1[3])
      {
        return std::sqrt((c0[0]-c1[0])*(c0[0]-c1[0]) + (c0[1]-c1[1])*(c0[1]-c1[1]) + (c0[2]-c1[2])*(c0[2]-c1[2]) );
      }

      void difference(double v01[3], double c0[3], double c1[3])
      {
        v01[0] = c0[0] - c1[0];
        v01[1] = c0[1] - c1[1];
        v01[2] = c0[2] - c1[2];
      }

      double dot(double c0[3], double c1[3])
      {
        return c0[0]*c1[0] + c0[1]*c1[1] + c0[2]*c1[2];
      }

      double plane_dot_product(double plane_point[3], double plane_normal[3], double point[3])
      {
        double normal[3]={0,0,0};
        normalize(plane_normal, normal);
        double dot = 0.0;
        for (int i = 0; i < 3; i++)
          {
            dot += (point[i] - plane_point[i])*normal[i];
          }
        return dot;
      }

      PlaneShock::PlaneShock()
      {
        plane_point_init[0]=0;
        plane_point_init[1]=0;
        plane_point_init[2]=0;
        plane_point[0]=0;
        plane_point[1]=0;
        plane_point[2]=0;
        plane_normal[0]=1;
        plane_normal[1]=0;
        plane_normal[2]=0;
        shock_width = 0.0;
      }

      void PlaneShock::setCurrentPlanePoint(double shock_displacement)
      {
        normalize(plane_normal);
        plane_point[0] = plane_point_init[0] + shock_displacement*plane_normal[0];
        plane_point[1] = plane_point_init[1] + shock_displacement*plane_normal[1];
        plane_point[2] = plane_point_init[2] + shock_displacement*plane_normal[2];
      }
      
      double PlaneShock::shock_function(double x)
      {
        // normalize by width
        return std::tanh(x/shock_width);
      }

      double shock_diff(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh,
                        stk::mesh::Entity node0, stk::mesh::Entity node1, double *coord0, double *coord1, PlaneShock& shock, double shock_displacement)
      {
        shock.setCurrentPlanePoint(shock_displacement);
        double *plane_point = shock.plane_point;
        double *plane_normal = shock.plane_normal;

        double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
        double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

        double v01[3] = {0,0,0};
        difference(v01, coord1, coord0);
        normalize(v01);
        normalize(plane_normal);
        double v01dotn = std::abs(dot(v01, plane_normal));

        double d01p = std::abs(dot_0)+std::abs(dot_1);
        dot_0 = shock.shock_function(dot_0);
        dot_1 = shock.shock_function(dot_1);

        if (nodal_refine_field)
          {
            double *fd0 = eMesh.field_data(nodal_refine_field, node0);
            double *fd1 = eMesh.field_data(nodal_refine_field, node1);
            fd0[0] = dot_0;
            fd1[0] = dot_1;
          }

        double d01 = distance(coord0, coord1);

        return (1 + 0*d01 + 0*d01p + 0*v01dotn)*std::abs(dot_0 - dot_1);
      }

      int shock_diff1(stk::mesh::FieldBase* nodal_refine_field, percept::PerceptMesh& eMesh,
                      stk::mesh::Entity node0, stk::mesh::Entity node1, 
                      double *coord0, double *coord1, PlaneShock& shock, double shock_displacement)
      {
        shock.setCurrentPlanePoint(shock_displacement);
        double *plane_point = shock.plane_point;
        double *plane_normal = shock.plane_normal;

        if (distance(coord0, coord1) < 1./200.)
          return DO_UNREFINE;

        double dot_0 = plane_dot_product(plane_point, plane_normal, coord0);
        double dot_1 = plane_dot_product(plane_point, plane_normal, coord1);

        if (nodal_refine_field)
          {
            double *fd0 = eMesh.field_data(nodal_refine_field, node0);
            double *fd1 = eMesh.field_data(nodal_refine_field, node1);
            fd0[0] = dot_0 < 0? -1 : 1;
            fd1[0] = dot_1 < 0? -1 : 1;
          }

        if (dot_0*dot_1 < 1.e-6)
          return DO_REFINE;
        else
          return DO_UNREFINE;

      }

      TEST(regr_localRefiner, break_tet_to_tet_N_5_ElementBased)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes_small.e");

            Local_Tet4_Tet4_N_HangingNode break_tet_to_tet_N(eMesh);
            
            eMesh.register_and_set_refine_fields();

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            eMesh.commit();

            SetRefineField set_ref_field(eMesh);
            elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, eMesh.m_refine_field);

            eMesh.save_as( output_files_loc+"local_tet_N_5_ElementBased_0_"+post_fix(p_size)+".e");

            ElementRefinePredicate erp(eMesh, 0, eMesh.m_refine_field, 0.0);

            TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            for (int ipass=0; ipass < 2; ipass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, eMesh.m_refine_field);

                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_ipass_"+toString(ipass)+"_"+post_fix(p_size)+".e");
              }

            breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_ElementBased_1_"+post_fix(p_size)+".e");
          }
      }

      TEST(regr_localRefiner, break_tet_to_tet_N_5_ElementBased_shock)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes_small.e");

            Local_Tet4_Tet4_N_HangingNode break_tet_to_tet_N(eMesh);

            eMesh.register_and_set_refine_fields();

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            eMesh.save_as( output_files_loc+"local_tet_N_5_EdgeBased_shock_0_"+post_fix(p_size)+".e");

            PlaneShock shock;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0,  eMesh.m_refine_field, 0.0, shock, 0.0);

            TransitionElementAdapter<ShockBasedRefinePredicate> breaker(srp, eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            for (int ipass=0; ipass < 3; ipass++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                breaker.doBreak();
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << std::endl;
                eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_ipass_"+toString(ipass)+"_"+post_fix(p_size)+".e");
              }

            breaker.deleteParentElements();
            eMesh.save_as(output_files_loc+"local_tet_N_5_EdgeBased_shock_1_"+post_fix(p_size)+".e");
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================


      TEST(regr_localRefiner, break_tet_to_tet_N_5_ElementBased_moving_shock)
      {
        const double rebalThreshold = 1.2;
        const int num_time_steps = 2;
        const bool save_intermediate=true;

        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+"cylinder_with_5_holes_small.e");

            Local_Tet4_Tet4_N_HangingNode break_tet_to_tet_N(eMesh);

            eMesh.register_and_set_refine_fields();

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.commit();

            eMesh.set_ioss_write_options("large");

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"moving_shock_"+post_fix(p_size)+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = 2.0;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;
            shock.shock_width = 1./5.0;

            ShockBasedRefinePredicate srp(nodal_refine_field, eMesh, 0, eMesh.m_refine_field, 0.0, shock, 0.0, 0.4);

            TransitionElementAdapter<ShockBasedRefinePredicate> breaker(srp, eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            if (rebalThreshold != 0.0)
              breaker.setDoRebalance(true, rebalThreshold);

            double delta_shock_displacement = 0.2;
            double shock_displacement = -2.0;
            int num_ref_passes = 1;
            int num_unref_passes = 1;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass << std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " moving_shock number elements= " << eMesh.get_number_elements()
                              << std::endl;
                    if (save_intermediate)
                      eMesh.save_as(output_files_loc+"tmp_moving_shock_ref_istep_ipass_"+toString(istep)+"_"+toString(ipass)+"_"+post_fix(p_size)+".e");
                  }

                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
                    breaker.buildUnrefineList(elements_to_unref);

                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << " unref list size= " << elements_to_unref.size() << std::endl;
                    breaker.unrefineTheseElements(elements_to_unref);
                    if (save_intermediate)
                      eMesh.save_as(output_files_loc+"tmp_moving_shock_unref_istep_ipass_"+toString(istep)+"_"+toString(iunref_pass)+"_"+post_fix(p_size)+".e");
                  }

                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"moving_shock_"+post_fix(p_size)+".e."+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"final_moving_shock_"+post_fix(p_size)+".e."+toString(num_time_steps) );
            for (int iunref=0; iunref < 2; iunref++)
              {
                std::cout << "moving shock unrefineAll pass= " << iunref << std::endl;
                breaker.unrefineAll();
              }
            std::cout << "moving shock deleteParentElements = "  << std::endl;
            breaker.deleteParentElements();
            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"final_unrefed_moving_shock_"+post_fix(p_size)+".e."+toString(num_time_steps) );
          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(regr_localRefiner, break_tet_to_tet_N_5_ElementBased_moving_shock_cyl_sidesets)
      {
        const std::string filename = "cylinder_tet4_0.e";
        const int num_time_steps = 4;  // 10 for stress testing
        
        const bool save_intermediate=false;
        const bool delete_parents=false;

        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        double shock_diff_criterion = 0.04;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            percept::PerceptMesh eMesh;
            eMesh.open(input_files_loc+filename);

            Local_Tet4_Tet4_N_HangingNode break_tet_to_tet_N(eMesh);

            eMesh.register_and_set_refine_fields();

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field       = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            stk::mesh::FieldBase* nodal_refine_field    = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);

            eMesh.commit();
            eMesh.delete_side_sets();
            eMesh.output_active_children_only(true);

            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"cyl_sidesets_moving_shock_"+post_fix(p_size)+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = -0.8;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;
            shock.shock_width = 1./25.0;

            ShockBasedRefinePredicate1 srp(nodal_refine_field, eMesh, 0, eMesh.m_refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            TransitionElementAdapter<ShockBasedRefinePredicate1> breaker(srp, eMesh, break_tet_to_tet_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.setDoLevelBasedUnrefinement(true);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.02;
            double shock_displacement = 0.0;
            int num_ref_passes = 1;
            int num_unref_passes = 1;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                srp.m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak();
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                SetRefineField set_ref_field(eMesh);
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, eMesh.m_refine_field);

                eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_ref_"+post_fix(p_size)+"_step_"+toString(istep)+".e");

                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << " number elements before unrefine= " << eMesh.get_number_elements() <<  std::endl;
                    ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
                    breaker.buildUnrefineList(elements_to_unref);

                    size_t unref_sz = elements_to_unref.size();
                    breaker.unrefinePass2(elements_to_unref);

                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " unref list size = " << unref_sz << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                eMesh.save_as(output_files_loc+"tmp_cyl_sidesets_moving_shock_unref_"+post_fix(p_size)+"_step_"+toString(istep)+".e");

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+"cyl_sidesets_moving_shock_"+post_fix(p_size)+".e."+toString(istep+1) );

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", istep);
                    if (istep==0)
                      eMesh.save_as("cyl_tet_anim.e");
                    else
                      eMesh.save_as("cyl_tet_anim.e-s"+std::string(buf));
                  }

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+"cyl_sidesets_final_moving_shock_"+post_fix(p_size)+".e."+toString(num_time_steps) );
            for (int iunref=0; iunref < 4; iunref++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker.unrefineAll();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                eMesh.save_as(output_files_loc+"cyl_sidesets_final_moving_shock_"+post_fix(p_size)+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
              }

            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"cyl_sidesets_final_unrefed_moving_shock_"+post_fix(p_size)+".e."+toString(num_time_steps) );

            // end_demo
          }
      }

      void tet_transition(PerceptMesh& eMesh)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        // 121
        std::string filename = "cylinder_tet4_0.e";

        eMesh.open(input_files_loc+filename);
        eMesh.output_active_children_only(true);

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Local_Tet4_Tet4_N_HangingNode localBreakPattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

            // for plotting, use doubles, for internal use, use int
            RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
            stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);
            eMesh.set_refine_level_field(&refine_level);

            {
              TransitionElementType& transition_element       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>(stk::topology::ELEMENT_RANK, "transition_element_3");
              stk::mesh::put_field_on_mesh( transition_element , eMesh.get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(transition_element, Ioss::Field::TRANSIENT);
              eMesh.set_transition_element_field(&transition_element);
            }
            {
              TransitionElementType& transition_element       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>(stk::topology::FACE_RANK, "transition_element");
              stk::mesh::put_field_on_mesh( transition_element , eMesh.get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(transition_element, Ioss::Field::TRANSIENT);
              eMesh.set_transition_element_field_2d(&transition_element);
            }

            eMesh.commit();
            std::cout << "tet_transition number elements after pass " << -1 << " = " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"tet_"+post_fix(p_size)+".0.e");
            eMesh.save_as("tet.e");

            AdaptedMeshVerifier amv(0);
            VERIFY_OP_ON(amv.isValid(eMesh, true), ==, true, "Invalid initial mesh");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, localBreakPattern, proc_rank_field);

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);


            //stk::mesh::Entity element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 121);
            int nref=3;
            for (int iref=0; iref < nref; ++iref)
              {
                stk::mesh::Entity element = eMesh.get_element(-.25, -.6, .98);
                stk::mesh::FieldTraits<RefineFieldType>::data_type *e_data = 0;
                if (eMesh.is_valid(element))
                  {
                    e_data = stk::mesh::field_data(*refine_field, element);
                    e_data[0] = 1.0;
                    //std::cout << "P[" << eMesh.get_rank() << "] found element = " << eMesh.identifier(element) << std::endl;
                  }

                breaker.refine();
                eMesh.save_as("tet.e-s000"+toString(iref+1));
                VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid refined mesh at step: "+toString(iref));
                std::cout << "tet_transition number elements after pass " << iref << " = " << eMesh.get_number_elements() << std::endl;

                std::cout << "tet_transition number elements after pass " << iref << " = " << eMesh.get_number_elements() << std::endl;

                if (1)
                  {
                    SetElementRefineFieldValue srf(eMesh, 0);
                    elementOpLoop(*eMesh.get_bulk_data(), srf, refine_field);

                  }
              }

            SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1);

            for (int iunref=0; iunref < 4; iunref++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field_val_unref_all, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker.unrefine();

                eMesh.save_as("tet.e-s000"+toString(nref+iunref+1));
                VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid unrefined mesh at step: "+toString(iunref));

                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref
                          << " tet_local number elements= " << eMesh.get_number_elements() << std::endl;
              }

          }
      }

      TEST(regr_localRefiner, tet_transition)
      {
        PerceptMesh eMesh(3);
        tet_transition(eMesh);
      }

      TEST(regr_localRefiner, tet_elem_based)
      {

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        PerceptMesh eMesh(3);

        std::string filename = "cylinder_tet4_0.e";

        eMesh.open(input_files_loc+filename);
        eMesh.output_active_children_only(true);

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Local_Tet4_Tet4_N_HangingNode localBreakPattern(eMesh);

            eMesh.register_and_set_refine_fields();

            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);

            eMesh.commit();

            eMesh.save_as( output_files_loc+"tet_"+post_fix(p_size)+".0.e");
            eMesh.save_as("tet.e");

            AdaptedMeshVerifier amv(0);
            VERIFY_OP_ON(amv.isValid(eMesh, true), ==, true, "Invalid initial mesh");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);
            TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, localBreakPattern, proc_rank_field);

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            int nref=3;
            for (int iref=0; iref < nref; ++iref)
              {
                stk::mesh::Entity element = eMesh.get_element(-.25, -.6, .98);
                RefineFieldType_type *e_data = 0;
                if (eMesh.is_valid(element))
                  {
                    e_data = stk::mesh::field_data(*eMesh.m_refine_field, element);
                    e_data[0] = 1.0;
                  }

                breaker.doBreak();
                eMesh.save_as("tet.e-s000"+toString(iref+1));
                VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid refined mesh at step: "+toString(iref));

                if (1)
                  {
                    SetElementRefineFieldValue srf(eMesh, 0);
                    elementOpLoop(*eMesh.get_bulk_data(), srf, eMesh.m_refine_field);
                  }
              }
          }
      }

      class SetElementField : public percept::ElementOp
      {
        IEdgeAdapter *m_iea;
      public:
        SetElementField(percept::PerceptMesh& eMesh, IEdgeAdapter* iea) : m_iea(iea) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          int unref = m_iea->markUnrefine(element);
          int ref_count = m_iea->markCountRefinedEdges(element);
          f_data[0] = ((unref & DO_UNREFINE ? -1 : 0));
          f_data[0] = (ref_count ? 1 : f_data[0]);


          return false;  // don't terminate the loop
        }
        virtual void init_elementOp() {}
        virtual void fini_elementOp() {}

      };



      template< class Breaker >
      static void do_moving_shock_test_square_sidesets(PerceptMesh& eMesh, int num_time_steps, bool save_intermediate=false, bool delete_parents=false,
                                                       std::string prefix="")
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        double shock_diff_criterion = 0.04;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Breaker break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);

            // just some debugging fields, not necessary
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));
            stk::mesh::FieldBase* nodal_refine_field = eMesh.add_field("nodal_refine_field", eMesh.node_rank(), scalarDimension);
            eMesh.add_field("normal_kept_deleted", eMesh.node_rank(), scalarDimension);

            RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
            stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);
            eMesh.set_refine_level_field(&refine_level);

            eMesh.commit();
            eMesh.output_active_children_only(true);

            //eMesh.delete_side_sets();
            std::cout << "moving_shock initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+prefix+"square_sidesets_moving_shock_"+post_fix(p_size)+".e.0");

            PlaneShock shock;
            shock.plane_point_init[0] = -0.8;
            shock.plane_point_init[1] = 0.0;
            shock.plane_point_init[2] = 0.0;
            shock.plane_normal[0] = 1;
            shock.plane_normal[1] = 0;
            shock.plane_normal[2] = 0;
            shock.shock_width = 1./25.0;

            ShockBasedRefinePredicate1 srp(nodal_refine_field, eMesh, 0, refine_field, 0.0, shock, 0.0, shock_diff_criterion);

            PredicateBasedEdgeAdapter<ShockBasedRefinePredicate1>
              breaker(srp,
                      eMesh, break_tri_to_tri_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            //breaker.setDoLevelBasedUnrefinement(true);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            double delta_shock_displacement = 0.02;
            double shock_displacement = 0.0;
            int num_ref_passes = 3;
            int num_unref_passes = 1;

            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak();

                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_moving_shock_ref_"+post_fix(p_size)+".e.s-"+toString(istep+1));

                if (1)
                  {
                    nodalOpLoop(*eMesh.get_bulk_data(), srp, nodal_refine_field);
                  }

                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;

                    if (istep == 12)
                      {
                        Util::setFlag(1256,true);
                      }
                    if (istep==0)
                      {
                        Util::setFlag(1257,true);
                      }

                    if (1)
                      {
                        //ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                        ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
                        breaker.buildUnrefineList(elements_to_unref);

                        breaker.unrefinePass2(elements_to_unref);
                      }

                    if (istep==11)
                      {
                        char buf1[1000];
                        sprintf(buf1, "%04d", iunref_pass);
                        eMesh.save_as("test1.e"+(iunref_pass==0?"":"-s"+std::string(buf1)));
                        if (iunref_pass==15) {
                          std::cout << "shock_displacement= " << shock_displacement << std::endl;
                          //exit(1);
                        }
                      }

                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
                  }
                if (1)
                  {
                    nodalOpLoop(*eMesh.get_bulk_data(), srp, nodal_refine_field);
                  }


                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", istep);
                    if (istep==0)
                      eMesh.save_as(prefix+"square_anim.e");
                    else
                      eMesh.save_as(prefix+"square_anim.e-s"+std::string(buf));
                  }

                eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_moving_shock_unref_"+post_fix(p_size)+".e");

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+prefix+"square_sidesets_moving_shock_"+post_fix(p_size)+".e.s-"+toString(istep+1) );

                shock_displacement += delta_shock_displacement;

              }

            eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_moving_shock_"+post_fix(p_size)+".e.s-"+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_moving_shock_"+post_fix(p_size)+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
                breaker.unrefineAll();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " moving_shock number elements= " << eMesh.get_number_elements() << std::endl;
              }

            std::cout << "moving_shock final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_unrefed_moving_shock_"+post_fix(p_size)+".e."+toString(num_time_steps) );

            // end_demo
          }
      }

      TEST(regr_localRefiner, break_tri_to_tri_N_5_EdgeBased_moving_shock_square_sidesets)
      {
        const int num_time_steps = 10;
        
        {
          // structured mesh
          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"square_tri3_0.e");
            do_moving_shock_test_square_sidesets<Local_Tri3_Tri3_N>(eMesh, num_time_steps, false, true, "str-");
          }

          // unstructured mesh
          {
            PerceptMesh eMesh2;
            eMesh2.open(input_files_loc+"square_tri3_uns_xformed.e");
            do_moving_shock_test_square_sidesets<Local_Tri3_Tri3_N>(eMesh2, num_time_steps, false, false, "uns-");
          }
        }
      }

      void do_point_source_square_sidesets(PerceptMesh& eMesh,
                                                  int num_time_steps,
                                                  double bad_angle, double bad_angle_transition_elements,
                                                  int nq_iter,
                                                  bool save_intermediate=false, bool delete_parents=false,
                                                  std::string prefix="")
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        std::vector<double> point_source_position(2);
        point_source_position[0] = 0.0;
        point_source_position[1] = 0.1;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);

            // just some debugging fields, not necessary
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));
            RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
            stk::mesh::put_field_on_mesh( refine_level,  eMesh.get_fem_meta_data()->universal_part(), nullptr);
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            eMesh.commit();
            //eMesh.output_active_children_only(true);

            //eMesh.delete_side_sets();
            std::cout << "point_source initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+prefix+"square_sidesets_point_source_"+post_fix(p_size)+".e.0");

            PointSourceRefinePredicate srp(eMesh, 0, refine_field, 0.0, point_source_position);

            QualityElementAdapter<PointSourceRefinePredicate>
              breaker(srp,
                      eMesh, break_tri_to_tri_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.setQualityCriteriaAngle(bad_angle, bad_angle_transition_elements);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            int num_ref_passes = 5;
            int num_unref_passes = 1;
            int cnt1=0;

            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", cnt1);
                if (cnt1==0)
                  eMesh.save_as(prefix+"ref.e");
                else
                  eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                ++cnt1;
              }

            num_time_steps = 1; // FIXME
            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                //breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {

                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak(nq_iter);
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " point_source number elements= " << eMesh.get_number_elements() << std::endl;

                    if (1)
                      {
                        char buf[1000];
                        sprintf(buf, "%04d", cnt1);
                        if (cnt1==0)
                          eMesh.save_as(prefix+"ref.e");
                        else
                          eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                        ++cnt1;
                      }

                  }

                eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_point_source_ref_"+post_fix(p_size)+".e.s-"+toString(istep+1));


                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;
                    //ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
                    breaker.buildUnrefineList(elements_to_unref);

                    breaker.unrefinePass2(elements_to_unref);
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " point_source number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", cnt1);
                    if (cnt1==0)
                      eMesh.save_as(prefix+"ref.e");
                    else
                      eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                    ++cnt1;
                  }

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", istep);
                    if (istep==0)
                      eMesh.save_as(prefix+"square_anim.e");
                    else
                      eMesh.save_as(prefix+"square_anim.e-s"+std::string(buf));
                  }

                eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_point_source_unref_"+post_fix(p_size)+".e");

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if (istep == num_time_steps-1 || save_intermediate)
                  eMesh.save_as(output_files_loc+prefix+"square_sidesets_point_source_"+post_fix(p_size)+".e.s-"+toString(istep+1) );

                //shock_displacement += delta_shock_displacement;
              }

            eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_point_source_"+post_fix(p_size)+".e.s-"+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_point_source_"+post_fix(p_size)+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
                breaker.unrefineAll();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " point_source number elements= " << eMesh.get_number_elements() << std::endl;

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", cnt1);
                    if (cnt1==0)
                      eMesh.save_as(prefix+"ref.e");
                    else
                      eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                    ++cnt1;
                  }

              }

            std::cout << "point_source final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_unrefed_point_source_"+post_fix(p_size)+".e."+toString(num_time_steps) );

            // end_demo
          }
        else
          {
            eMesh.commit();
          }
        eMesh.close();
      }

      void do_point_source_square_sidesets_template(PerceptMesh& eMesh,
                                                    int num_time_steps,
                                                    double bad_angle, double bad_angle_transition_elements,
                                                    int nq_iter,
                                                    bool save_intermediate=false, bool delete_parents=false,
                                                    std::string prefix="")
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        std::vector<double> point_source_position(2);
        point_source_position[0] = 0.0;
        point_source_position[1] = 0.1;

        bool extra_output = false;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);

            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

            {
              RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
              stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);
              eMesh.set_refine_level_field(&refine_level);
            }

            {
              TransitionElementType& transition_element       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>(stk::topology::ELEMENT_RANK, "transition_element");
              stk::mesh::put_field_on_mesh( transition_element , eMesh.get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(transition_element, Ioss::Field::TRANSIENT);
              eMesh.set_transition_element_field(&transition_element);
            }

            Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
            eMesh.commit();
            //eMesh.output_active_children_only(true);

            //eMesh.delete_side_sets();
            std::cout << "point_source initial number elements= " << eMesh.get_number_elements() << std::endl;

            if (extra_output) eMesh.save_as( output_files_loc+prefix+"square_sidesets_point_source_"+post_fix(p_size)+".e.0");

            PointSourceRefinePredicate srp(eMesh, 0, refine_field, 0.0, point_source_position);

            PredicateTemplateAdapter<PointSourceRefinePredicate>
              breaker(srp,
                      eMesh, break_tri_to_tri_N, proc_rank_field);

            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            // x,y,z: [-.08,.08],[-.08,.08],[-.5,-.2]
            int num_ref_passes = 5;
            int num_unref_passes = 1;
            int cnt1=0;

            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", cnt1);
                if (cnt1==0)
                  eMesh.save_as(prefix+"ref.e");
                else
                  eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                ++cnt1;
              }

            num_time_steps = 1; // FIXME
            for (int istep = 0; istep < num_time_steps; istep++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] istep= " << istep << std::endl;

                //breaker.getRefinePredicate().m_shock_displacement = shock_displacement;

                for (int ipass=0; ipass < num_ref_passes; ipass++)
                  {

                    std::cout << "P[" << eMesh.get_rank() << "] ipass= " << ipass <<  std::endl;
                    breaker.doBreak(nq_iter);
                    std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " point_source number elements= " << eMesh.get_number_elements() << std::endl;

                    if (1)
                      {
                        char buf[1000];
                        sprintf(buf, "%04d", cnt1);
                        if (cnt1==0)
                          eMesh.save_as(prefix+"ref.e");
                        else
                          eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                        ++cnt1;
                      }

                  }

                if (extra_output) eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_point_source_ref_"+post_fix(p_size)+".e.s-"+toString(istep+1));


                for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
                  {
                    std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass <<  std::endl;
                    //ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
                    ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
                    breaker.buildUnrefineList(elements_to_unref);

                    breaker.unrefinePass2(elements_to_unref);
                    std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " point_source number elements= " << eMesh.get_number_elements() << std::endl;
                  }

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", cnt1);
                    if (cnt1==0)
                      eMesh.save_as(prefix+"ref.e");
                    else
                      eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                    ++cnt1;
                  }

                if (extra_output)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", istep);
                    if (istep==0)
                      eMesh.save_as(prefix+"square_anim.e");
                    else
                      eMesh.save_as(prefix+"square_anim.e-s"+std::string(buf));
                  }

                if (extra_output) eMesh.save_as(output_files_loc+prefix+"tmp_square_sidesets_point_source_unref_"+post_fix(p_size)+".e");

                if (delete_parents && istep == num_time_steps-1)
                  {
                    breaker.deleteParentElements();
                  }
                if ((istep == num_time_steps-1 || save_intermediate) && extra_output)
                  eMesh.save_as(output_files_loc+prefix+"square_sidesets_point_source_"+post_fix(p_size)+".e.s-"+toString(istep+1) );

                //shock_displacement += delta_shock_displacement;
              }

            if (extra_output) eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_point_source_"+post_fix(p_size)+".e.s-"+toString(num_time_steps) );
            for (int iunref=0; iunref < 10; iunref++)
              {
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                if (extra_output) eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_point_source_"+post_fix(p_size)+"_unrefAll_pass_"+toString(iunref)+".e."+toString(num_time_steps) );
                breaker.unrefineAll();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " point_source number elements= " << eMesh.get_number_elements() << std::endl;

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", cnt1);
                    if (cnt1==0)
                      eMesh.save_as(prefix+"ref.e");
                    else
                      eMesh.save_as(prefix+"ref.e-s"+std::string(buf));
                    ++cnt1;
                  }

              }

            std::cout << "point_source final number elements= " << eMesh.get_number_elements() << std::endl;
            if (extra_output) eMesh.save_as(output_files_loc+prefix+"square_sidesets_final_unrefed_point_source_"+post_fix(p_size)+".e."+toString(num_time_steps) );

            // end_demo
          }
        else
          {
            eMesh.commit();
          }
        eMesh.close();
      }

      TEST(regr_localRefiner, point_source_tri)
      {
        if (1)
          {
            int num_time_steps = 1;
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"square_tri3_uns_xformed.e");
            do_point_source_square_sidesets(eMesh, num_time_steps, 30., 30., 5, false, false, "point-source-uns-");
          }
      }

      TEST(regr_localRefiner, point_source_tri_template)
      {
        if (1)
          {
            int num_time_steps = 1;
            PerceptMesh eMesh;
            int nq_iter = 10;
            eMesh.open(input_files_loc+"square_tri3_uns_xformed.e");
            do_point_source_square_sidesets_template(eMesh, num_time_steps, 30., 30., nq_iter, false, false, "point-source-uns-");
          }
      }

      static void do_tet_edge_test(unsigned ntets, bool do_part_for_edges)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        std::string ntetstr = toString(ntets);

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1)
          {
            {
              unsigned npts=5;
              //unsigned ntets=2;
              if (ntets == 1) npts=4;
              static  SingleTetFixture::Point node_coord_data[  ] = {
                { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 }, {1, 1, 1} };

              // Hard coded tetra node ids for all the tetra nodes in the entire mesh
              static  SingleTetFixture::TetIds tetra_node_ids[] = {
                { 1, 2, 3, 4}, {2, 3, 4, 5} };

              percept::SingleTetFixture mesh(pm, false, npts, node_coord_data, ntets, tetra_node_ids);

              stk::io::put_io_part_attribute(  mesh.m_block_tet );

              //stk::mesh::Part& part = mesh.m_metaData.declare_part_with_topology("elems", stk::topology::TET_4);

              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              eMesh.save_as("tmp.e");
            }

            {
              percept::PerceptMesh eMesh;
              eMesh.open("tmp.e");

              int scalarDimension = 0; // a scalar
              stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
              RefineFieldType *refine_field = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

              const std::string part_for_edges_name = "tmp_part_for_edges";
              stk::mesh::Part * part_for_edges = 0;
              if (do_part_for_edges) {
                part_for_edges = &eMesh.get_fem_meta_data()->declare_part_with_topology(part_for_edges_name, stk::topology::LINE_2);
              }

              if (part_for_edges) stk::io::put_io_part_attribute(*part_for_edges);
              Local_Tet4_Tet4_N break_tet_to_tet_N(eMesh);

              eMesh.commit();

              {
                std::vector<size_t> counts;
                stk::mesh::comm_mesh_counts(*eMesh.get_bulk_data() , counts);
                std::cout << "Initial mesh has  "
                          << counts[0] << " nodes, "
                          << counts[1] << " edges, "
                          << counts[2] << " faces, "
                          << counts[3] << " elements" << std::endl;
              }

              if (part_for_edges)
                {
                  stk::mesh::create_edges(*eMesh.get_bulk_data());

                  {
                    std::vector<size_t> counts;
                    stk::mesh::comm_mesh_counts(*eMesh.get_bulk_data() , counts);
                    std::cout << "Mesh after create_edges has  "
                              << counts[0] << " nodes, "
                              << counts[1] << " edges, "
                              << counts[2] << " faces, "
                              << counts[3] << " elements" << std::endl;
                  }

                  if (1)
                    {
                      stk::mesh::PartVector add_parts(1, part_for_edges), remove_parts;
                      std::vector<stk::mesh::Entity> edges;
                      const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( eMesh.edge_rank() );
                      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                        {
                          stk::mesh::Bucket & bucket = **k ;
                          const unsigned num_edges_in_bucket = bucket.size();
                          for (unsigned iEdge = 0; iEdge < num_edges_in_bucket; iEdge++)
                            {
                              stk::mesh::Entity edge = bucket[iEdge];
                              edges.push_back(edge);
                            }
                        }

                      eMesh.get_bulk_data()->modification_begin();
                      for (unsigned iEdge = 0; iEdge < edges.size(); iEdge++)
                        {
                          eMesh.get_bulk_data()->change_entity_parts( edges[iEdge], add_parts, remove_parts );
                        }
                      stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
                      eMesh.get_bulk_data()->modification_end();
                    }
                  //eMesh.print_info("after create edges", 2);
                }

              std::cout << "nele= " << eMesh.get_number_elements() << std::endl;
              std::cout << "nnode= " << eMesh.get_number_nodes() << std::endl;
              stk::mesh::Entity element_0 = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);
              RefineFieldType_type *fd = stk::mesh::field_data(*refine_field, element_0);
              fd[0] = 1.0;
              std::cout << "print_entity= " << std::endl;
              eMesh.print_entity(std::cout, element_0, refine_field);
              eMesh.print_entity_compact(element_0, refine_field);

              eMesh.save_as( output_files_loc+ntetstr+"_tet_with_edges_0.e");

              ElementRefinePredicate erp(eMesh, 0, refine_field, 0.0);

              PredicateBasedElementAdapter<ElementRefinePredicate>
                breaker(erp, eMesh, break_tet_to_tet_N, proc_rank_field);

              breaker.setRemoveOldElements(false);
              breaker.setAlwaysInitializeNodeRegistry(false);

              breaker.doBreak();
              eMesh.save_as(output_files_loc+ntetstr+"_tet_with_edges_1.e");


              breaker.deleteParentElements();
              eMesh.save_as(output_files_loc+ntetstr+"_tet_with_edges_2.e");
              std::cout << "nele= " << eMesh.get_number_elements() << std::endl;
              std::cout << "nnode= " << eMesh.get_number_nodes() << std::endl;

              eMesh.dump_vtk(ntetstr+"-tet.vtk",true);

              if (1)
                {
                  // sanity on mesh counts; overall time
                  std::vector<size_t> counts;
                  stk::mesh::comm_mesh_counts(*eMesh.get_bulk_data() , counts);

                  std::cout << "Mesh has  "
                            << counts[0] << " nodes, "
                            << counts[1] << " edges, "
                            << counts[2] << " faces, "
                            << counts[3] << " elements" << std::endl;
                }
            }
          }
      }

      TEST(regr_localRefiner, tet_with_edges)
      {
        do_tet_edge_test(1, true);
        // do_tet_edge_test(2, false);
        // do_tet_edge_test(1, true);
        // do_tet_edge_test(2, true);

      }

      //  ====================================================================================================================================
      //  ====================================================================================================================================
      //  Local quad hanging-node refinement
      //  ====================================================================================================================================
      //  ====================================================================================================================================

      class SetElementFieldQuadCorner : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase *m_field;
        std::vector<double> m_bbox;
      public:
        SetElementFieldQuadCorner(percept::PerceptMesh& eMesh, stk::mesh::FieldBase *field, std::vector<double> bbox)
          : m_eMesh(eMesh), m_field(field), m_bbox(bbox)
        {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          stk::mesh::FieldBase* coord_field = m_eMesh.get_coordinates_field();

          const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );

          unsigned num_node = elem_nodes.size();
          double c[2] = {0,0};
          for (unsigned inode=0; inode < num_node; inode++)
            {
              stk::mesh::Entity node = elem_nodes[ inode ].entity();
              double *c_data = m_eMesh.field_data(coord_field, node);
              c[0] += c_data[0]/double(num_node);
              c[1] += c_data[1]/double(num_node);
            }
          if ((m_bbox[0] < c[0] && c[0] < m_bbox[1]) &&
              (m_bbox[2] < c[1] && c[1] < m_bbox[3]) )
            {
              f_data[0] = 1;
            }
          else
            {
              f_data[0] = -1;
            }
          // FIXME tmp
          //f_data[0] = 1.0;
          return false;  // don't terminate the loop
        }
        virtual void init_elementOp()
        {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }
        virtual void fini_elementOp() {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }

      };

      class SetElementFieldQuadPointSource : public percept::ElementOp
      {
        percept::PerceptMesh& m_eMesh;
        stk::mesh::FieldBase *m_field;

      public:
        std::vector<double> m_pointSourcePosition;

        SetElementFieldQuadPointSource(percept::PerceptMesh& eMesh, stk::mesh::FieldBase *field, std::vector<double>& pointSourcePosition)
          : m_eMesh(eMesh), m_field(field), m_pointSourcePosition(pointSourcePosition) {
        }

        virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
        {
          RefineFieldType_type *f_data = stk::mesh::field_data(*static_cast<RefineFieldType *>(field), element);
          //stk::mesh::FieldBase* coord_field = m_eMesh.get_coordinates_field();

          const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);
          //int spatialDimension = m_eMesh.get_spatial_dim();
          CellTopology cell_topo(cell_topo_data);
          const percept::MyPairIterRelation elem_nodes (m_eMesh, element,stk::topology::NODE_RANK);
          //VERIFY_OP_ON(elem_nodes.size(), ==, 3, "only for tris");

          CoordinatesFieldType* coordField = m_eMesh.get_coordinates_field();

          unsigned numSubDimNeededEntities = cell_topo_data->edge_count;

          //int count=0;
          unsigned ref_count=0, unref_count=0;
          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              stk::mesh::Entity node0 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[0]].entity();
              stk::mesh::Entity node1 = elem_nodes[cell_topo_data->edge[iSubDimOrd].node[1]].entity();
              double * const coord0 = m_eMesh.field_data( *coordField , node0 );
              double * const coord1 = m_eMesh.field_data( *coordField , node1 );
              double cross =
                (coord1[0] - coord0[0])*(m_pointSourcePosition[1] - coord0[1])
                - (coord1[1] - coord0[1])*(m_pointSourcePosition[0] - coord0[0]);

              if (cross > 0)
                ++ref_count;
              else
                ++unref_count;
            }
          if (ref_count == numSubDimNeededEntities)
            f_data[0] = 1;
          else
            f_data[0] = -1;
          return false;  // don't terminate the loop
        }
        virtual void init_elementOp()
        {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }
        virtual void fini_elementOp() {
          std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }

      };


      template<class LocalBreakPattern>
      static void do_quad_local_corner_refine_sidesets(PerceptMesh& eMesh, int num_time_steps,
                                                       const std::string& prefix="",
                                                       std::vector<double> bbox=std::vector<double>(),
                                                       int num_ref_passes = 8, int num_unref_passes = 4)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            LocalBreakPattern localBreakPattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

            // for plotting, use doubles, for internal use, use int
            RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
            stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);
            eMesh.set_refine_level_field(&refine_level);

            {
              TransitionElementType& transition_element       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>(stk::topology::ELEMENT_RANK, "transition_element");
              stk::mesh::put_field_on_mesh( transition_element , eMesh.get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(transition_element, Ioss::Field::TRANSIENT);
              eMesh.set_transition_element_field(&transition_element);
            }

            eMesh.commit();

            std::cout << "qual_local initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"quad_local_"+post_fix(p_size)+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, localBreakPattern, proc_rank_field);

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            std::vector<double> pts(2); pts[0] = -0.97; pts[1] = 0.21;
            SetElementFieldQuadPointSource set_ref_field(eMesh, refine_field, pts);

            for (int ipass=0; ipass < num_ref_passes; ipass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);
                breaker.refine();
              }

            for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                breaker.unrefine();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " quad_local number elements= " << eMesh.get_number_elements() << std::endl;
              }
            SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1);

            eMesh.save_as(output_files_loc+"quad_tmp_square_sidesets_quad_local_unref_"+post_fix(p_size)+".e");

            eMesh.save_as(output_files_loc+"quad_square_sidesets_final_quad_local_"+post_fix(p_size)+".e.s-"+toString(num_time_steps) );

            for (int iunref=0; iunref < (num_unref_passes ? 10 : 0); iunref++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field_val_unref_all, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker.unrefine();

                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " quad_local number elements= " << eMesh.get_number_elements() << std::endl;
              }

            std::cout << "quad_local final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"quad_square_sidesets_final_unrefed_quad_local_"+post_fix(p_size)+".e."+toString(num_time_steps) );
          }
      }

      // quad transition element

      TEST(regr_localRefiner, quad_transition)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        {
          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(-1,1, -1, 1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();
            eMesh.save_as(input_files_loc+"quad_square_quad4_0.e");
          }

          {
            PerceptMesh eMesh;
            eMesh.set_ioss_write_options("large");
            eMesh.open(input_files_loc+"quad_square_quad4_0.e");

            eMesh.output_active_children_only(true);

            double bb[4] = {-0.1,0.1,-0.1,0.1};
            std::vector<double> bbox(&bb[0], &bb[0]+4);

            do_quad_local_corner_refine_sidesets<Local_Quad4_Quad4_N_Transition>(eMesh, 10, "quad-het-", bbox);
          }
        }
      }

      //  ====================================================================================================================================
      //  ====================================================================================================================================
      //  Local tri hanging-node refinement
      //  ====================================================================================================================================
      //  ====================================================================================================================================

      TEST(regr_localRefiner, point_source_tri_hanging_node)
      {
        if (1)
          {
            PerceptMesh eMesh;
            eMesh.set_ioss_write_options("large");
            eMesh.open(input_files_loc+"square_tri3_uns_xformed.e");

            do_quad_local_corner_refine_sidesets<Local_Tri3_Tri3_N_HangingNode>(eMesh, 10, "point-source-uns-");
          }
      }


      TEST(regr_localRefiner, tri_transition_unit)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (1)
          {
            const unsigned n = 2;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(-1,1, -1, 1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();
            eMesh.save_as(input_files_loc+"tri_transition_0.e");
          }

        PerceptMesh eMesh;
        eMesh.open(input_files_loc+"tri_transition_0.e");
        eMesh.output_active_children_only(true);

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 4)
          {
            Local_Tri3_Tri3_N_HangingNode localBreakPattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

            // for plotting, use doubles, for internal use, use int
            RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
            stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            {
              TransitionElementType& transition_element       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>(stk::topology::ELEMENT_RANK, "transition_element");
              stk::mesh::put_field_on_mesh( transition_element , eMesh.get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(transition_element, Ioss::Field::TRANSIENT);
              eMesh.set_transition_element_field(&transition_element);
            }

            eMesh.commit();

            eMesh.save_as( output_files_loc+"tri_"+post_fix(p_size)+".0.e");
            eMesh.save_as("tri.e");

            AdaptedMeshVerifier amv(0);
            VERIFY_OP_ON(amv.isValid(eMesh, true), ==, true, "Invalid initial mesh");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            HangingNodeAdapter<ElementRefinePredicate> *breaker_p=0;
            TransitionElementAdapter<ElementRefinePredicate> *trea_p=0;

            trea_p = new TransitionElementAdapter<ElementRefinePredicate>(erp, eMesh, localBreakPattern, proc_rank_field);
            breaker_p = trea_p;

            HangingNodeAdapter<ElementRefinePredicate>& breaker = *breaker_p;

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            stk::mesh::Entity element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);
            RefineFieldType_type *e_data = 0;
            if (eMesh.is_valid(element))
              {
                e_data = stk::mesh::field_data(*refine_field, element);
                e_data[0] = 1.0;
              }

            breaker.refine();
            eMesh.save_as("tri.e-s0001");
            VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid refined mesh at step: 0");

            if (1)
              {
                SetElementRefineFieldValue srf(eMesh, 0);
                elementOpLoop(*eMesh.get_bulk_data(), srf, refine_field);

                stk::mesh::EntityId elemId = (p_size == 4 ? 34 : 10);
                element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elemId);
                if (eMesh.is_valid(element))
                  {
                    e_data = stk::mesh::field_data(*refine_field, element);
                    e_data[0] = 1.0;
                  }

                if (1)
                  {
                    element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 4);
                    if (eMesh.is_valid(element))
                      {
                        e_data = stk::mesh::field_data(*refine_field, element);
                        e_data[0] = 1.0;
                      }
                  }

                breaker.refine();
                eMesh.save_as("tri.e-s0002");
                VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid refined mesh at step: 1");
              }

            if (1 && p_size == 4)
              {
                SetElementRefineFieldValue srf(eMesh, 1);
                elementOpLoop(*eMesh.get_bulk_data(), srf, refine_field);

                breaker.refine();
                eMesh.save_as("tri.e-s0003");
                VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid refined mesh at step: 3");
              }

            eMesh.save_as( output_files_loc+"tri_transition.e");

            delete breaker_p;

            // end_demo
          }
      }

      TEST(regr_localRefiner, tri_transition_demo)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (1)
          {
            const unsigned n = 2;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = false;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(-1,1, -1, 1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();
            eMesh.save_as(input_files_loc+"tri_transition_0.e");
          }

        PerceptMesh eMesh;
        eMesh.open(input_files_loc+"tri_transition_0.e");
        eMesh.output_active_children_only(true);

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 4)
          {
            Local_Tri3_Tri3_N_HangingNode localBreakPattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

            // for plotting, use doubles, for internal use, use int
            RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
            stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            {
              TransitionElementType& transition_element       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>(stk::topology::ELEMENT_RANK, "transition_element");
              stk::mesh::put_field_on_mesh( transition_element , eMesh.get_fem_meta_data()->universal_part(), nullptr);
              stk::io::set_field_role(transition_element, Ioss::Field::TRANSIENT);
              eMesh.set_transition_element_field(&transition_element);
            }

            eMesh.commit();

            eMesh.save_as( output_files_loc+"tri_"+post_fix(p_size)+".0.e");
            eMesh.save_as("tri.e");

            AdaptedMeshVerifier amv(0);
            VERIFY_OP_ON(amv.isValid(eMesh, true), ==, true, "Invalid initial mesh");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            HangingNodeAdapter<ElementRefinePredicate> *breaker_p=0;
            TransitionElementAdapter<ElementRefinePredicate> *trea_p=0;

            trea_p = new TransitionElementAdapter<ElementRefinePredicate>(erp, eMesh, localBreakPattern, proc_rank_field, true);
            breaker_p = trea_p;

            HangingNodeAdapter<ElementRefinePredicate>& breaker = *breaker_p;

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            stk::mesh::Entity element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);
            RefineFieldType_type *e_data = 0;
            if (eMesh.is_valid(element))
              {
                e_data = stk::mesh::field_data(*refine_field, element);
                e_data[0] = 1.0;
              }

            breaker.refine();
            eMesh.save_as("tri.e-s0001");
            VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid refined mesh at step: 0");

            if (1)
              {
                SetElementRefineFieldValue srf(eMesh, 0);
                elementOpLoop(*eMesh.get_bulk_data(), srf, refine_field);

                if (1)
                  {
                    stk::mesh::EntityId elemId = 11;
                    element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elemId);
                    if (eMesh.is_valid(element))
                      {
                        e_data = stk::mesh::field_data(*refine_field, element);
                        e_data[0] = 1.0;
                      }
                  }
                if (1)
                  {
                    stk::mesh::EntityId elemId = 8;
                    element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elemId);
                    if (eMesh.is_valid(element))
                      {
                        e_data = stk::mesh::field_data(*refine_field, element);
                        e_data[0] = 1.0;
                      }
                  }

                breaker.refine();
                eMesh.save_as("tri.e-s0002");
                VERIFY_OP_ON(amv.isValid(eMesh, false), ==, true, "Invalid refined mesh at step: 1");
              }

            //if (delete_parents)
            //  breaker.deleteParentElements();
            eMesh.save_as( output_files_loc+"tri_transition.e");


            delete breaker_p;

            // end_demo
          }
      }

      //  ====================================================================================================================================
      //  ====================================================================================================================================
      //  Local hybrid transition element refinement
      //  ====================================================================================================================================
      //  ====================================================================================================================================

      template<class LocalBreakPattern>
      static void do_hybrid_quad_tri_local_corner_refine_sidesets(PerceptMesh& eMesh, int num_time_steps,
                                                                  bool save_intermediate=false, bool delete_parents=false, const std::string& prefix="",
                                                                  std::vector<double> bbox=std::vector<double>())
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            LocalBreakPattern localBreakPattern(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            eMesh.register_and_set_refine_fields();
            RefineFieldType *refine_field       = eMesh.m_refine_field;
            eMesh.commit();

            std::cout << "hybrid_quad_tri initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"hybrid_quad_local_"+post_fix(p_size)+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            //erp.m_mark_centroid_always = true;
            TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, localBreakPattern, proc_rank_field);

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            SetElementFieldQuadCorner set_ref_field_qc(eMesh, refine_field, bbox);
            percept::ElementOp* set_ref_field_ptr = &set_ref_field_qc;
            percept::ElementOp& set_ref_field = *set_ref_field_ptr;

            int num_ref_passes = 4;
            int num_unref_passes = 3;
            int iplot=0;

            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", iplot);
                if (iplot == 0)
                  eMesh.save_as(prefix+"quad_tri_anim.e");
                else
                  eMesh.save_as(prefix+"quad_tri_anim.e-s"+std::string(buf));
                ++iplot;
              }

            for (int ipass=0; ipass < num_ref_passes; ipass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                //eMesh.save_as(output_files_loc+"quad_anim_set_field_"+post_fix(p_size)+".e.s-"+toString(ipass+1));

                breaker.refine();

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as(prefix+"quad_tri_anim.e");
                    else
                      eMesh.save_as(prefix+"quad_tri_anim.e-s"+std::string(buf));
                    ++iplot;
                  }
              }

            for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                breaker.unrefine();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " quad_local number elements= " << eMesh.get_number_elements() << std::endl;
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as(prefix+"quad_tri_anim.e");
                    else
                      eMesh.save_as(prefix+"quad_tri_anim.e-s"+std::string(buf));
                    ++iplot;
                  }
              }

            SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1);

            eMesh.save_as(output_files_loc+"quad_tmp_tri_sidesets_quad_local_unref_"+post_fix(p_size)+".e");

            eMesh.save_as(output_files_loc+"quad_tri_sidesets_final_quad_local_"+post_fix(p_size)+".e.s-"+toString(num_time_steps) );

            for (int iunref=0; iunref < 10; iunref++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field_val_unref_all, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker.unrefine();

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as(prefix+"quad_tri_anim.e");
                    else
                      eMesh.save_as(prefix+"quad_tri_anim.e-s"+std::string(buf));
                    ++iplot;
                  }

                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " quad_local number elements= " << eMesh.get_number_elements() << std::endl;
              }

            if (delete_parents)
              breaker.deleteParentElements();
            std::cout << "quad_local final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"quad_tri_sidesets_final_unrefed_quad_local_"+post_fix(p_size)+".e."+toString(num_time_steps) );

            // end_demo
          }
      }

      TEST(regr_localRefiner, local_hybrid_quad4_tri3_n_transition)
      {
        bool do_test = true;
        //stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"lshape_hybrid.e");

            //eMesh.output_active_children_only(true);

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)

            double bb[4] = {-0.3,0.3,-0.3,0.3};
            //double bb[4] = {-1.0,-0.78,-1.0,-0.78};
            //double bb[4] = {-1.0,-0.78,-0.1,0.1};
            std::vector<double> bbox(&bb[0], &bb[0]+4);

            do_hybrid_quad_tri_local_corner_refine_sidesets<Local_Quad4_Tri3_Hybrid_Transition>(eMesh, 10, false, false, "", bbox);
          }
        }
      }

      //  ====================================================================================================================================
      //  ====================================================================================================================================
      //  Local hex hanging-node refinement
      //  ====================================================================================================================================
      //  ====================================================================================================================================

      template<class LocalBreakPattern>
      static void do_hex_local_corner_refine_sidesets(PerceptMesh& eMesh, int num_time_steps, bool save_intermediate=false, bool delete_parents=false, int num_ref_passes = 6,  int num_unref_passes = 4,
                                                      std::vector<double> bbox_new = std::vector<double>(), std::string block_name = "")
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        // for some reason, writing the animation files creates an Ioss exception, turn off for now
        bool do_anim = false;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3)
          {
            LocalBreakPattern break_hex_to_hex_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("refine_field", stk::topology::ELEMENT_RANK, scalarDimension));

            // for plotting, use doubles, for internal use, use int
            RefineLevelType& refine_level       = eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
            stk::mesh::put_field_on_mesh( refine_level , eMesh.get_fem_meta_data()->universal_part(), nullptr);
            stk::io::set_field_role(refine_level, Ioss::Field::TRANSIENT);

            eMesh.register_and_set_refine_fields();
            eMesh.commit();

            std::cout << "hex_local initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"hex_local_"+post_fix(p_size)+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            HangingNodeAdapter<ElementRefinePredicate> *breaker = 0;
            if (typeid(LocalBreakPattern) == typeid(Local_Hex8_Hex8_N_Transition) || typeid(LocalBreakPattern) == typeid(Local_Hybrid_3D))
              {
                breaker = new TransitionElementAdapter<ElementRefinePredicate>(erp, eMesh, break_hex_to_hex_N, proc_rank_field);
              }
            else
              breaker = new HangingNodeAdapter<ElementRefinePredicate>(erp, eMesh, break_hex_to_hex_N, proc_rank_field);

            // special for local adaptivity
            breaker->setRemoveOldElements(false);
            breaker->setAlwaysInitializeNodeRegistry(false);

            double bbox[] = {-0.1, 0.1, -0.1, 0.1, -0.1, 0.1};
            std::vector<double> vbbox(&bbox[0], &bbox[0]+6);
            if (bbox_new.size() == 6)
              vbbox = bbox_new;
            std::shared_ptr<percept::ElementOp> set_ref_field;
            if (!block_name.size())
              {
                set_ref_field.reset(new SetElementFieldBoundingBox(eMesh, vbbox, refine_field));
              }
            else
              {
                stk::mesh::Part *part = eMesh.get_fem_meta_data()->get_part(block_name);
                VERIFY_OP_ON(part, !=, 0, "bad block_name= " + block_name);
                stk::mesh::Selector sel = *part;
                set_ref_field.reset(new SetElementFieldSelector(eMesh, sel, refine_field));
              }

            int iplot=0;
            if (do_anim)
              {
                char buf[1000];
                sprintf(buf, "%04d", iplot);
                if (iplot == 0)
                  eMesh.save_as("hex_cube_anim.e");
                else
                  eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                ++iplot;
              }
            for (int ipass=0; ipass < num_ref_passes; ipass++)
              {

                int nel = eMesh.get_number_elements();
                if (eMesh.get_rank() == 0)
                  std::cout << "P[" << eMesh.get_rank() << "] before refine ipass= " << ipass << " hex_local number elements= "
                            << nel << std::endl;

                elementOpLoop(*eMesh.get_bulk_data(), *set_ref_field, refine_field);
                //eMesh.save_as(output_files_loc+"hex_anim_set_field_"+post_fix(p_size)+".e.s-"+toString(ipass+1));
                if (do_anim)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", ipass);
                    if (ipass == 0)
                      eMesh.save_as("hex_set_field.e");
                    else
                      eMesh.save_as("hex_set_field.e-s"+std::string(buf));
                  }

                // {node-, edge-, face-neighors}
                breaker->refine();
                //breaker.doBreak();

                MPI_Barrier( MPI_COMM_WORLD );
                //VERIFY_OP_ON(true,==,false,"here");
                std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " hex_local number elements= "
                          << eMesh.get_number_elements() << std::endl;
                //eMesh.save_as("square_anim."+toString(ipass+1)+".e");
                if (do_anim)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("hex_cube_anim.e");
                    else
                      eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                    ++iplot;
                  }

              }

            for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), *set_ref_field, refine_field);

                breaker->unrefine();
                std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " hex_local number elements= " << eMesh.get_number_elements() << std::endl;
                if (do_anim)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("hex_cube_anim.e");
                    else
                      eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                    ++iplot;
                  }
              }

            SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1);

            eMesh.save_as(output_files_loc+"hex_tmp_square_sidesets_hex_local_unref_"+post_fix(p_size)+".e");

            eMesh.save_as(output_files_loc+"hex_cube_sidesets_final_hex_local_"+post_fix(p_size)+".e.s-"+toString(num_time_steps) );

            for (int iunref=0; iunref < (num_unref_passes? 10 : 0); iunref++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field_val_unref_all, refine_field);
                std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker->unrefine();

                if (do_anim)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    if (iplot == 0)
                      eMesh.save_as("hex_cube_anim.e");
                    else
                      eMesh.save_as("hex_cube_anim.e-s"+std::string(buf));
                    ++iplot;
                  }

                std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " hex_local number elements= " << eMesh.get_number_elements() << std::endl;
              }

            if (delete_parents)
              breaker->deleteParentElements();
            std::cout << "hex_local final number elements= " << eMesh.get_number_elements() << std::endl;
            eMesh.save_as(output_files_loc+"hex_cube_sidesets_final_unrefed_hex_local_"+post_fix(p_size)+".e."+toString(num_time_steps) );
            //exit(123);

            // end_demo
            delete breaker;
          }
      }

      TEST(regr_localRefiner, break_hex_to_hex_N_5_ElementBased_hex_local_square_sidesets)
      {
        {

          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n, nz = n;

            percept::PerceptMesh eMesh(3u);
            std::string gmesh_spec = toString(nx)+"x"+toString(ny)+"x"+toString(nz)+"|bbox:-1,-1,-1,1,1,1";
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();

            eMesh.save_as(input_files_loc+"hex_cube_hex4_0.e");
          }

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"hex_cube_hex4_0.e");

            ScalarFieldType & weight_field =  eMesh.get_fem_meta_data()->declare_field<ScalarFieldType>(stk::topology::ELEMENT_RANK, "rebalance_weights");
            stk::mesh::put_field_on_mesh( weight_field , eMesh.get_fem_meta_data()->universal_part() , nullptr);
            stk::io::set_field_role(weight_field, Ioss::Field::TRANSIENT);

            eMesh.output_active_children_only(true);

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            bool delete_parents = true;
            do_hex_local_corner_refine_sidesets<Local_Hex8_Hex8_N>(eMesh, 10, false, delete_parents, 6, 0);

            percept::RebalanceMesh rb(eMesh, &weight_field, true);
            double imb_before = rb.compute_imbalance();
            double imb_after = rb.rebalance();

            if (eMesh.get_rank() == 0) std::cout << "imbalance before= " << imb_before << " imbalance after= " << imb_after << std::endl;

            EXPECT_NEAR(1.0, imb_after, 1.0e-6);
          }
        }
      }

      TEST(regr_localRefiner, tri_rebalance)
      {
        if (1)
          {
            bool debug = false;
            PerceptMesh eMesh;
            eMesh.set_ioss_write_options("large");
            eMesh.open(input_files_loc+"square_tri3_uns_xformed.e");
            eMesh.register_and_set_refine_fields();
            ScalarFieldType & weight_field =  eMesh.get_fem_meta_data()->declare_field<ScalarFieldType>(stk::topology::ELEMENT_RANK, "rebalance_weights");
            stk::mesh::put_field_on_mesh( weight_field , eMesh.get_fem_meta_data()->universal_part() , nullptr);
            stk::io::set_field_role(weight_field, Ioss::Field::TRANSIENT);

            eMesh.output_active_children_only(false);

            int num_time_steps = 10;
            const std::string prefix="tri_rebal_anim";
            int nref=8;
            int nunref=0;
            std::vector<double> bbox;
            do_quad_local_corner_refine_sidesets<Local_Tri3_Tri3_N_HangingNode>(eMesh, num_time_steps, prefix, bbox, nref, nunref);
            eMesh.save_as("tri_adapted.e");

            if (1)
              {
                percept::RebalanceMesh rb(eMesh, &weight_field, debug);
                double imb_before = rb.compute_imbalance();
                double imb_after = rb.rebalance();

                if (eMesh.get_rank() == 0) std::cout << "imbalance before= " << imb_before << " imbalance after= " << imb_after << std::endl;

                EXPECT_NEAR(1.0, imb_after, 0.33);
              }

            if (1)
              {
                Tri3_Tri3_4 tri3_tri3_4(eMesh);
                Refiner rf(eMesh, tri3_tri3_4);
                rf.deleteParentElements();
                std::vector<double> load_factor;
                eMesh.get_load_factor(load_factor, true, "after rebalance and delete_parents");
                eMesh.save_as("tri_rebal.e");
              }

          }
      }

      TEST(regr_localRefiner, tri_rebalance_sidesets)
      {
        if (1)
          {
            bool debug = false;
            PerceptMesh eMesh;
            eMesh.set_ioss_write_options("large");
            eMesh.open(input_files_loc+"square_tri3_0.e");
            eMesh.register_and_set_refine_fields();
            ScalarFieldType & weight_field =  *eMesh.m_weights_field;

            eMesh.output_active_children_only(false);

            int num_time_steps = 10;
            const std::string prefix="tri_rebal_ss_anim";
            int nref=8;
            int nunref=0;
            std::vector<double> bbox;
            do_quad_local_corner_refine_sidesets<Local_Tri3_Tri3_N_HangingNode>(eMesh, num_time_steps, prefix, bbox, nref, nunref);
            eMesh.save_as("tri_ss_adapted.e");

            if (1)
              {
                percept::RebalanceMesh rb(eMesh, &weight_field, debug);
                double imb_before = rb.compute_imbalance();
                double imb_after = rb.rebalance();

                if (eMesh.get_rank() == 0) std::cout << "imbalance before= " << imb_before << " imbalance after= " << imb_after << std::endl;

                EXPECT_NEAR(1.0, imb_after, 0.70);
              }
          }
      }

      TEST(regr_localRefiner, hex_rebalance)
      {
        bool do_test = true;
        //stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          bool do_bootstrap_mesh = true;
          if (do_bootstrap_mesh)
          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n, nz = n;

            percept::PerceptMesh eMesh(3u);
            std::string gmesh_spec = toString(nx)+"x"+toString(ny)+"x"+toString(nz)+"|bbox:-1,-1,-1,1,1,1";
            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();

            eMesh.save_as(input_files_loc+"hex_cube_hex4_0.e");
            //eMesh.print_info("test1",2);
            //eMesh.dump_vtk("sqhex3.vtk",false);
          }

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"hex_cube_hex4_0.e");

            ScalarFieldType & weight_field =  eMesh.get_fem_meta_data()->declare_field<ScalarFieldType>(stk::topology::ELEMENT_RANK, "rebalance_weights");
            stk::mesh::put_field_on_mesh( weight_field , eMesh.get_fem_meta_data()->universal_part() , nullptr);
            stk::io::set_field_role(weight_field, Ioss::Field::TRANSIENT);

            eMesh.output_active_children_only(true);

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            bool delete_parents = false;
            do_hex_local_corner_refine_sidesets<Local_Hex8_Hex8_N>(eMesh, 10, false, delete_parents, 6, 0);
            eMesh.save_as("hex_adapted.e");

            bool debug_print = false;
            percept::RebalanceMesh rb(eMesh, &weight_field, debug_print);
            double imb_before = rb.compute_imbalance();
            double imb_after = rb.rebalance();

            if (eMesh.get_rank() == 0) std::cout << "imbalance before= " << imb_before << " imbalance after= " << imb_after << std::endl;

            EXPECT_NEAR(1.0, imb_after, 1.86);
          }
        }
      }

      TEST(regr_localRefiner, hex_transition)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 3) return;

        //const unsigned n = (p_size <= 3 ? 3 : p_size);
        const unsigned n = 5;
        const unsigned nx = n , ny = n, nz = n;
        //const unsigned nx = 1 , ny = 1, nz = n;
        int nref = 2;
        int nunref = 1;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          bool do_bootstrap_mesh = true;
          if (do_bootstrap_mesh)
          {
            percept::PerceptMesh eMesh(3u);
            //std::string gmesh_spec = toString(nx)+"x"+toString(ny)+"x"+toString(nz)+"|bbox:-1,-1,-1,1,1,1";
            std::string gmesh_spec = toString(nx)+"x"+toString(ny)+"x"+toString(nz)+std::string("|bbox:-1,-1,-1,1,1,1|sideset:xXyYzZ");

            eMesh.new_mesh(percept::GMeshSpec(gmesh_spec));
            eMesh.commit();

            eMesh.save_as(input_files_loc+"hex_cube_hex8_0.e");
          }

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"hex_cube_hex8_0.e");

            eMesh.output_active_children_only(false);

            //int num_time_steps = 10;  // 10 for stress testing
            //for (int istep=1; istep <= num_time_steps; istep++)
            double bbox[] = {-0.9, -0.7, -0.9, -0.7, -0.9, -0.7, -0.9, -0.7};
            std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

            do_hex_local_corner_refine_sidesets<Local_Hex8_Hex8_N_Transition>(eMesh, 1, false, true, nref, nunref, vbbox);
          }
        }
      }

      TEST(regr_localRefiner, hex_transition_7block)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;

        int nref = 3;
        int nunref = 1;

        if (do_test) {

          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"hex-7block.g");

            eMesh.output_active_children_only(false);

            double bl = 0.25;
            double br = 0.75;
            double bbox[] = {bl,br,bl,br,bl,br};
            std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

            do_hex_local_corner_refine_sidesets<Local_Hex8_Hex8_N_Transition>(eMesh, 1, false, true, nref, nunref, vbbox, "block_1");
            //do_hex_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, false, true, nref, nunref, vbbox);
          }
        }
      }

      //  ====================================================================================================================================
      //  ====================================================================================================================================
      //  Local pyramid transition refine
      //  ====================================================================================================================================
      //  ====================================================================================================================================

      template<class LocalBreakPattern>
      static void do_het_local_corner_refine_sidesets(PerceptMesh& eMesh, int num_time_steps, std::vector<double>& vbbox, const std::string& file,
                                                      bool save_intermediate=false, bool delete_parents=false, int num_ref_passes = 6,  int num_unref_passes = 4,
                                                      std::vector<std::string> do_wedge_bl= std::vector<std::string>(), int num_loops = 1)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3 || p_size >= 8)
          {
            LocalBreakPattern break_het_to_het_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = 0; 

            eMesh.register_and_set_refine_fields();
            refine_field = eMesh.m_refine_field;

            eMesh.commit();

            std::cout << "het_local initial number elements= " << eMesh.get_number_elements() << std::endl;

            eMesh.save_as( output_files_loc+"het_local_"+post_fix(p_size)+".e.0");

            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            TransitionElementAdapter<ElementRefinePredicate> *breaker_p;
            stk::mesh::Selector wedge_selector;
            if (do_wedge_bl.size())
              {
                stk::mesh::PartVector pv;
                for (unsigned ii = 0; ii < do_wedge_bl.size(); ++ii)
                  {
                    stk::mesh::Part* part = eMesh.get_fem_meta_data()->get_part(do_wedge_bl[ii]);
                    VERIFY_OP_ON(part, !=, 0, "can't find wedge part "+do_wedge_bl[ii]);
                    VERIFY_OP_ON(part->topology(), ==, stk::topology::WEDGE_6, "block has wrong topology: "+do_wedge_bl[ii]);
                    if (eMesh.get_rank() == 0) std::cout << "TEA_SpecialWedgeRefinement:: found part = " << part->name() << std::endl;
                    pv.push_back(part);
                  }
                wedge_selector = stk::mesh::selectUnion(pv);
                breaker_p = new TEA_SpecialWedgeRefinement<ElementRefinePredicate>(erp, eMesh, break_het_to_het_N, proc_rank_field, &wedge_selector);
              }
            else
              breaker_p = new TransitionElementAdapter<ElementRefinePredicate>(erp, eMesh, break_het_to_het_N, proc_rank_field);
            TransitionElementAdapter<ElementRefinePredicate>& breaker = *breaker_p;
            if (eMesh.get_rank() == 0)
              std::cout << "type= " << eMesh.demangle(typeid(breaker).name()) << std::endl;

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            // FIXME
            //breaker.setAddChildrenToParts(false);

            SetElementFieldBoundingBox set_ref_field(eMesh, vbbox, refine_field);

            int iplot=0;
            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", iplot);
                if (iplot == 0)
                  eMesh.save_as(file+"_anim.e", double(iplot));
                else
                  eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                ++iplot;
              }
            for (int loop=0; loop < num_loops; ++loop)
              {
            for (int ipass=0; ipass < num_ref_passes; ipass++)
              {

                int nel = eMesh.get_number_elements();
                if (eMesh.get_rank() == 0)
                  std::cout << "P[" << eMesh.get_rank() << "] before refine ipass= " << ipass << " local number elements= "
                            << nel << std::endl;

                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", ipass);
                    eMesh.save_as(file+"_set_field.e-s"+std::string(buf), double(iplot));
                    ++iplot;
                  }

                // {node-, edge-, face-neighors}
                breaker.refine();

                //eMesh.print_info("refined", 2);

                nel = eMesh.get_number_elements();
                if (eMesh.get_rank() == 0)
                  std::cout << "P[" << eMesh.get_rank() << "] done... ipass= " << ipass << " local number elements= "
                            << nel << std::endl;
                //eMesh.dump_elements_compact();
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                    ++iplot;
                  }

              }

            for (int iunref_pass=0; iunref_pass < num_unref_passes; iunref_pass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                breaker.unrefine();
                int nel = eMesh.get_number_elements();
                if (eMesh.get_rank() == 0)
                  std::cout << "P[" << eMesh.get_rank() << "] done... iunref_pass= " << iunref_pass << " local number elements= " << nel << std::endl;
                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                    ++iplot;
                  }
              }
              }

            SetElementRefineFieldValue set_ref_field_val_unref_all(eMesh, -1);

            eMesh.save_as(output_files_loc+file+"_tmp_square_sidesets_local_unref_"+post_fix(p_size)+".e");

            eMesh.save_as(output_files_loc+file+"_sidesets_final_local_"+post_fix(p_size)+".e.s-"+toString(num_time_steps) );

            for (int iunref=0; iunref < (num_unref_passes? 10 : 0); iunref++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field_val_unref_all, refine_field);
                if (eMesh.get_rank() == 0)
                  std::cout << "P[" << eMesh.get_rank() << "] iunrefAll_pass= " << iunref <<  std::endl;
                breaker.unrefine();

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                    ++iplot;
                  }

                int nel = eMesh.get_number_elements();
                if (eMesh.get_rank() == 0)
                  std::cout << "P[" << eMesh.get_rank() << "] done... iunrefAll_pass= " << iunref << " local number elements= " << nel << std::endl;
              }


            if (delete_parents)
              breaker.deleteParentElements();
            int nel = eMesh.get_number_elements();
            if (eMesh.get_rank() == 0)
              std::cout << "local final number elements= " << nel << std::endl;
            eMesh.save_as(output_files_loc+file+"_sidesets_final_unrefed_local_"+post_fix(p_size)+".e."+toString(num_time_steps) );

            delete breaker_p;

            // end_demo
          }
      }

      template<class LocalBreakPattern>
      static void do_refine_seq_spec_elems(PerceptMesh& eMesh, const std::string& file,
                                           std::vector<stk::mesh::EntityId>& elem_ids)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size==1 || p_size == 3 || p_size >= 8)
          {
            LocalBreakPattern break_het_to_het_N(eMesh);
            int scalarDimension = 0; // a scalar
            stk::mesh::FieldBase* proc_rank_field    = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
            RefineFieldType *refine_field       = 0; 

            // for plotting, use doubles, for internal use, use int
            eMesh.register_and_set_refine_fields();
            refine_field = eMesh.m_refine_field;

            eMesh.commit();
            stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

            ElementRefinePredicate erp(eMesh, &univ_selector, refine_field, 0.0);
            TransitionElementAdapter<ElementRefinePredicate> *breaker_p;
            breaker_p = new TransitionElementAdapter<ElementRefinePredicate>(erp, eMesh, break_het_to_het_N, proc_rank_field);
            TransitionElementAdapter<ElementRefinePredicate>& breaker = *breaker_p;

            // special for local adaptivity
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);

            SetElementRefineFieldValue set_ref_field(eMesh, 0.0);

            int iplot=0;
            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", iplot);
                if (iplot == 0)
                  eMesh.save_as(file+"_anim.e", double(iplot));
                else
                  eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                ++iplot;
              }
            int num_ref_passes = elem_ids.size();

            for (int ipass=0; ipass < num_ref_passes; ipass++)
              {
                elementOpLoop(*eMesh.get_bulk_data(), set_ref_field, refine_field);

                stk::mesh::Entity element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), elem_ids[ipass]);
                VERIFY_OP_ON(eMesh.is_valid(element), ==, true, "bad element");
                RefineFieldType_type *ref_field = stk::mesh::field_data(*refine_field, element);
                ref_field[0] = 1;

                // {node-, edge-, face-neighors}
                breaker.refine();

                if (1)
                  {
                    char buf[1000];
                    sprintf(buf, "%04d", iplot);
                    eMesh.save_as(file+"_anim.e-s"+std::string(buf), double(iplot));
                    ++iplot;
                  }

              }

            delete breaker_p;

            // end_demo
          }
      }


      TEST(regr_localRefiner, pyr_transition)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 3) return;

        int nref = 1;
        int nunref = 1;

        if (do_test) {

          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"heterogeneous_sideset_0.e");
          eMesh.output_active_children_only(false);

          //int num_time_steps = 10;  // 10 for stress testing
          //for (int istep=1; istep <= num_time_steps; istep++)

          //double bbox[] = {0,2,0,1,-2,-1};  // all prisms
          double bbox[] = {0,1,0,1,-2,-1};  // just one prism
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          //do_het_local_corner_refine_sidesets<Local_Pyr5_Pyr5_N_Transition>(eMesh, 1, vbbox, "pyr", false, true, nref, nunref);
          do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "pyr", false, true, nref, nunref);
        }
      }

      TEST(regr_localRefiner, wedge_transition)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 3) return;

        int nref = 1;
        int nunref = 1;

        if (do_test) {

          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"heterogeneous_sideset_0.e");
          eMesh.output_active_children_only(false);

          //int num_time_steps = 10;  // 10 for stress testing
          //for (int istep=1; istep <= num_time_steps; istep++)

          double bbox[] = {0,2,1,2,-1,0}; // all
          //double bbox[] = {1,2,1,2,-1,0};  // just one wedge
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          do_het_local_corner_refine_sidesets<Local_Wedge6_Wedge6_N_Transition>(eMesh, 1, vbbox, "wedge", false, true, nref, nunref);
        }

      }

      TEST(regr_localRefiner, hybrid_transition)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 3) return;

        int nref = 3;
        int nunref = 1;

        if (!do_test) return;

        // test 1
        if (1)
        {
          PerceptMesh eMesh;
          //eMesh.open(input_files_loc+"heterogeneous_sideset_0.e");
          eMesh.open(input_files_loc+"heterogeneous_0.e");
          eMesh.output_active_children_only(false);

          double bbox[] = {0,3,0,2,-2,0}; // all
          //double bbox[] = {1,2,1,2,-1,0};  // just one wedge
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "hybrid", false, true, nref, nunref);
        }

        // test 2
        bool tmp = true;
        if (tmp)
        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"heterogeneous_sideset_0.e");
          //eMesh.open(input_files_loc+"heterogeneous_0.e");
          eMesh.output_active_children_only(false);

          //double bbox[] = {0,3,0,2,-2,0}; // all
          double bbox[] = {1,2,1,2,-1,0};  // just one wedge
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "hybrid1", false, true, nref, nunref);
        }


        // test 3
        if (tmp)
        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"heterogeneous_sideset_0.e");
          //eMesh.open(input_files_loc+"heterogeneous_0.e");
          eMesh.output_active_children_only(false);

          double bbox[] = {0,1,0,1,-2,-1};  // just one prism
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "hybrid2", false, true, nref, nunref);
        }

        // FIXME - repeat with no sidesets

        // FIXME - random marks?
      }

      TEST(regr_localRefiner, hybrid_transition_multi)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 3) return;

        int nref = 3;
        int nunref = 1;

        if (!do_test) return;

        bool tmp = true;
        if (tmp)
        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"heterogeneous_sideset_multi_0.e");
          //eMesh.open(input_files_loc+"heterogeneous_0.e");
          eMesh.output_active_children_only(false);

          //double bbox[] = {0,3,0,2,-2,0}; // all
          double bbox[] = {1,2,1,2,-1,0};  // just one wedge
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "hybrid1", false, true, nref, nunref);
        }
      }


      TEST(regr_localRefiner, hybrid_tet_wedge_transition)
      {
        int nref = 2;
        int nunref = 2;

        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"tet4_wedge6_mesh.g");
          eMesh.output_active_children_only(true);

          double bbox[] = {-0.5,0.5,-0.5,0.5,0,1};
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "hybrid_tet4_wedge6", false, true, nref, nunref);
        }
      }

      TEST(regr_localRefiner, wedge_boundary_layer_special)
      {
        int nref = 1; //2;
        int nunref = 1; //2;

        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"tet4_wedge6_mesh.g");
          //eMesh.output_active_children_only(true);
          eMesh.commit();
          eMesh.delete_side_sets();
          eMesh.save_as(input_files_loc+"tet4_wedge6_mesh_no_sidesets.g");
        }

        // partial refinement
        if (1)
          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"tet4_wedge6_mesh_no_sidesets.g");
            eMesh.output_active_children_only(true);

            double bbox[] = {-0.5,0.5,-0.5,0.5,0,1};
            std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

            std::vector<std::string> do_wedge_bl;
            do_wedge_bl.push_back("block_1");
            do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "wedge_bl", false, true, nref, nunref, do_wedge_bl);
          }

      }

      TEST(regr_localRefiner, wedge_boundary_layer_special_1)
      {
        int nref = 1; //2;
        int nunref = 1; //2;
        int num_loops = 3; //2;

        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"tet4_wedge6_mesh.g");
          //eMesh.output_active_children_only(true);
          eMesh.commit();
          eMesh.delete_side_sets();
          eMesh.save_as(input_files_loc+"tet4_wedge6_mesh_no_sidesets.g");
        }

        // partial refinement
        if (1)
          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"tet4_wedge6_mesh_no_sidesets.g");
            eMesh.output_active_children_only(true);

            double bbox[] = {-0.5,0.5,-0.5,0.5,0,1};
            std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

            std::vector<std::string> do_wedge_bl;
            do_wedge_bl.push_back("block_1");
            do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "wedge_bl_2nd_ref", false, true, nref, nunref, do_wedge_bl, num_loops);
          }
      }

      TEST(regr_localRefiner, wedge_boundary_layer_special_sidesets)
      {
        int nref = 3; //3;
        int nunref = 3; //3;

        // uniformly refine all elements (wedges get excluded)
        if (0)
          {
            PerceptMesh eMesh;
            eMesh.open(input_files_loc+"tet4_wedge6_mesh.g");
            eMesh.output_active_children_only(true);

            double bbox[] = {-1,1,-1,1,0,1};
            std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

            std::vector<std::string> do_wedge_bl;
            do_wedge_bl.push_back("block_1");
            do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "wedge_bl_unif_sidesets", false, true, nref, nunref, do_wedge_bl);
          }
      }

      TEST(regr_localRefiner, pyr_tet_special_case_transition)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;

        int nref = 1;
        int nunref = 0;

        double coord[][3] = {

          {3.56232, -0.0407922, 0.0426477},
          {3.56252, -0.041082, 0.0424578},
          {3.56236, -0.0412391, 0.0421968},
          {3.56268, -0.0409035, 0.0428227},
          {3.56287, -0.0412162, 0.0426341},

          {-3.56252 +2*3.56232, -1*(-0.041082)+2*(-0.0407922), -1*0.0424578 +2*0.0426477},

          {-1*3.56287 +2*3.56268, -1*(-0.0412162) + 2*(-0.0409035), -1*0.0426341 + 2*0.0428227},

          {3.56236+0.001, -0.0412391-0.001, 0.0421968+0.001}

        };
        (void) coord;

        stk::mesh::EntityIdVector pyramid_node_ids_0[2] {
          {2, 5, 4, 1, 3},
            {1, 4, 7, 6, 8}
        };

        unsigned nnodes=8, npyr=2;

        {
          bool doCommit = false;
          percept::PyramidFixture mesh(MPI_COMM_WORLD, doCommit, false);
          mesh.m_metaData.commit();
          if (1) mesh.init(nnodes, npyr, coord, &pyramid_node_ids_0[0]);
          mesh.populate();

          bool isCommitted = true;
          percept::PerceptMesh em1(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);

          em1.save_as("pyr_0.e");
          bool isBad = em1.check_mesh_volumes(true, 0.0, 1);
          EXPECT_FALSE(isBad);

          em1.close();
        }

        std::string input_mesh = "pyr_0.e";

        if (do_test) {

          PerceptMesh eMesh;
          eMesh.open_read_only(input_mesh);
          eMesh.output_active_children_only(false);

          double C[3]={0,0,0};
          unsigned ielem = 2;
          stk::mesh::Entity pyr = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), ielem);
          stk::mesh::Entity pyr1 = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 1);

          VERIFY_OP_ON(eMesh.is_valid(pyr), ==, true, "bad pyr");
          percept::computeCentroid(pyr, &C[0], *eMesh.get_coordinates_field());
          if (1)
            {
              VolumeUtil jacA;
              shards::CellTopology cell_topo(eMesh.get_cell_topology(pyr));
              double volScale = jacA.getJacobianToVolumeScale(cell_topo);

              double jacobian = 0.0;
              jacA(jacobian, eMesh, pyr1, eMesh.get_coordinates_field(), eMesh.get_cell_topology(pyr));
              double cellVol = jacobian*volScale;
              for (unsigned i=0; i < 5; ++i)
                std::cout << "cellVol = " << cellVol << " i= " << i << " det= " << jacA.m_detJ[i] << std::endl;

            }
          eMesh.reopen();


          double ep=1.e-5;
          double bbox[] = {C[0]-ep,C[0]+ep,C[1]-ep,C[1]+ep,C[2]-ep,C[2]+ep};
          std::vector<double> vbbox(&bbox[0], &bbox[0]+6);

          bool delParents = false;
          do_het_local_corner_refine_sidesets<Local_Hybrid_3D>(eMesh, 1, vbbox, "pyr", false, delParents, nref, nunref);

          bool isBadRef = eMesh.check_mesh_volumes(true, 0.0, 1);
          EXPECT_FALSE(isBadRef);
          eMesh.dump_vtk("mesh.vtk");

        }
      }

      TEST(regr_localRefiner, pyr_find_valid_centroid)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;

        double coord[][3] = {
          { 3.56324, 0.0406481, 0.0438404 },
          { 3.56313, 0.0406429, 0.0437659 },
          { 3.56354, 0.0408522, 0.0440267 },
          { 3.56362, 0.0408083, 0.0441883 },
          { 3.56333, 0.0408, 0.0438086 },
          {3.56337, 0.0407503, 0.043926} // centroid of 1st 5 points (pyramid)
        };
        (void) coord;

        stk::mesh::EntityIdVector pyramid_node_ids_0[2] {
          {1, 2, 3, 4, 5},
          {1, 2, 3, 4, 6}
        };
        stk::mesh::EntityIdVector tet_node_ids_0[4] {
          {6, 1, 2, 5},
            {6, 2, 3, 5},
              {6, 3, 4, 5},
                {6, 4, 1, 5}
        };

        unsigned nnodes=6, npyr=2, ntet=4;

        {
          bool doCommit = false;
          percept::HeterogeneousFixture mesh(MPI_COMM_WORLD, doCommit, false);
          mesh.m_metaData.commit();
          //mesh.init(nnodes, npyr, coord, &pyramid_node_ids_0[0]);
          mesh.init(nnodes, coord,
                    0, 0,
                    ntet, &tet_node_ids_0[0],
                    0, 0,
                    npyr, &pyramid_node_ids_0[0]);
          mesh.populate();

          bool isCommitted = true;
          percept::PerceptMesh em1(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);

          em1.save_as("pyr_fvc_0.e");
          bool isBad = em1.check_mesh_volumes(true, 0.0, 1);
          EXPECT_FALSE(isBad);

          em1.close();
        }

        std::string input_mesh = "pyr_fvc_0.e";

        if (do_test) {

          PerceptMesh eMesh;
          eMesh.open_read_only(input_mesh);
          bool isBad = eMesh.check_mesh_volumes(true, 0.0, 1);
          EXPECT_FALSE(isBad);

          eMesh.output_active_children_only(false);
          unsigned ielem = 5;
          stk::mesh::Entity pyr = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), ielem);
          if (1)
            {
              double sc_volume[5];
              FiniteVolumeMesh3D fvm(*eMesh.get_bulk_data());
              fvm.elementVolume(pyr, sc_volume);
              for (unsigned ii=0; ii < 5; ++ii)
                {
                  std::cout << "sc_volume[" << ii << "] = " << sc_volume[ii] << std::endl;
                }
            }
          unsigned childIds[] = {1,2,3,4,6};

          const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(eMesh.element_rank() + 1u);
          std::vector<stk::mesh::Entity> fts;
          eMesh.get_bulk_data()->modification_begin();
          eMesh.createEntities(FAMILY_TREE_RANK, 5, fts);
          eMesh.get_bulk_data()->modification_end();
          eMesh.get_bulk_data()->modification_begin();
          for (unsigned ii=0; ii < 5; ++ii)
            {
              stk::mesh::Entity child = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), childIds[ii]);
              UniformRefinerPatternBase::set_parent_child_relations(eMesh, pyr, child, fts[ii], ii);
            }
          eMesh.get_bulk_data()->modification_end();
          int ndiv=5;
          double C[3]={0,0,0};
          percept::computeCentroid(pyr, &C[0], *eMesh.get_coordinates_field());
          FindValidCentroid findValidCentroid(eMesh, ndiv, false, true);
          std::vector<stk::mesh::Entity> nodes, pyrs;
          pyrs.push_back(pyr);
          stk::mesh::get_entities_through_relations(*eMesh.get_bulk_data(), pyrs, stk::topology::NODE_RANK, nodes);
          stk::mesh::Entity cnode = eMesh.get_bulk_data()->get_entity(eMesh.node_rank(), 6);
          bool didChange = findValidCentroid.findCentroid(pyr, &C[0], nodes, cnode);
          std::cout << "didChange= " << didChange << std::endl;
        }

      }

      TEST(regr_localRefiner, wedge_recurse_tea)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size > 1) return;

        double coord[][3] = {

          {0,0,0}, {1,0,0}, {0,1,0},
          {0,0,1}, {1,0,1}, {0,1,1},
          {1,1,0}, {1,1,1},
          {0,0,2}, {1,0,2}, {0,1,2}
        };
        (void) coord;

        stk::mesh::EntityIdVector wedge_node_ids_0[3] {
          {1,2,3,4,5,6},
            {2,7,3,5,8,6},
              {4,5,6,9,10,11}
        };


        unsigned nnodes=11, nwedge=3;

        {
          bool doCommit = false;
          percept::HeterogeneousFixture mesh(MPI_COMM_WORLD, doCommit, false);
          mesh.m_metaData.commit();
          if (1) mesh.init(nnodes, coord,
                           0, 0,
                           0, 0,
                           nwedge, wedge_node_ids_0,
                           0, 0);
          mesh.populate();

          bool isCommitted = true;
          percept::PerceptMesh em1(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);

          em1.save_as("wedge_0.e");
          bool isBad = em1.check_mesh_volumes(true, 0.0, 1);
          EXPECT_FALSE(isBad);

          em1.close();
        }

        std::string input_mesh = "wedge_0.e";

        if (do_test) {

          PerceptMesh eMesh;
          eMesh.open_read_only(input_mesh);
          eMesh.output_active_children_only(false);

          if (1)
            {
              eMesh.reopen();

              std::vector<stk::mesh::EntityId> ids;
              ids.push_back(1);
              ids.push_back(7);

              do_refine_seq_spec_elems<Local_Hybrid_3D>(eMesh, "wedge", ids);
            }

          bool isBadRef = eMesh.check_mesh_volumes(true, 0.0, 1);
          EXPECT_FALSE(isBadRef);
          eMesh.dump_vtk("mesh.vtk");

        }
      }

      //=============================================================================
      // Encore/Percept mimic
      //=============================================================================
      //=============================================================================

      class Mimic_Encr_Entity_Marker {
      public:
        Mimic_Encr_Entity_Marker(percept::PerceptMesh & pMesh, const stk::mesh::FieldBase * marker_field, const stk::mesh::EntityRank entity_rank)
        : my_pMesh(pMesh), my_marker_field(marker_field), my_entity_rank(entity_rank) {}
        void update_markers();
        Int get_marker(stk::mesh::Entity entity);

      protected:
        percept::PerceptMesh & my_pMesh;
        const stk::mesh::FieldBase * my_marker_field;
        const stk::mesh::EntityRank my_entity_rank;
        std::vector<stk::mesh::Entity> my_entities;
        std::vector<Int> my_markers;
      };

      class Mimic_Encr_Percept_Edge_Adapter : public percept::IEdgeAdapter
      {
      public:
        Mimic_Encr_Percept_Edge_Adapter(Mimic_Encr_Entity_Marker & element_marker, percept::PerceptMesh & pMesh, percept::UniformRefinerPatternBase & bp, stk::mesh::FieldBase * proc_rank_field=0)
        : percept::IEdgeAdapter(pMesh, bp, proc_rank_field), my_element_marker(element_marker) {}
        virtual int markEdge(const stk::mesh::Entity element, unsigned which_edge, stk::mesh::Entity node0, stk::mesh::Entity node1,
          double *coord0, double *coord1, std::vector<int>* existing_edge_marks);
      protected:
        Mimic_Encr_Entity_Marker & my_element_marker;
      };

      int
      Mimic_Encr_Percept_Edge_Adapter::markEdge(
          const stk::mesh::Entity element,
          unsigned which_edge,
          stk::mesh::Entity node0,
          stk::mesh::Entity node1,
          double *coord0,
          double *coord1,
          std::vector<int>* existing_edge_marks)
      {
        std::vector<stk::mesh::Entity> nodes;
        nodes.push_back(node0);
        nodes.push_back(node1);
        std::vector<stk::mesh::Entity> edge_elems;
        get_entities_through_relations(*m_eMesh.get_bulk_data(), nodes, stk::topology::ELEMENT_RANK, edge_elems);
        ThrowRequire(!edge_elems.empty());

        bool all_marked_for_refinement = true;
        bool all_marked_for_unrefinement = true;
        for ( UInt ie = 0 ; ie < edge_elems.size(); ++ie )
        {
          stk::mesh::Entity edge_elem = edge_elems[ie];
          ThrowRequire(m_eMesh.is_valid(edge_elem));

          Int element_marker = my_element_marker.get_marker(edge_elem);
          if (element_marker <= 0) all_marked_for_refinement = false;
          if (element_marker >= 0) all_marked_for_unrefinement = false;
        }
        ThrowRequire(!(all_marked_for_refinement && all_marked_for_unrefinement));

        int mark=percept::DO_NOTHING;
        if (all_marked_for_refinement)
        {
          mark=percept::DO_REFINE;
        }
        else if (all_marked_for_unrefinement)
        {
          mark= percept::DO_UNREFINE;
        }
        return mark;
      }

      void
      Mimic_Encr_Entity_Marker::update_markers()
      {
        my_entities.clear();
        stk::mesh::get_entities( *my_pMesh.get_bulk_data(), my_entity_rank, my_entities );
        my_markers.resize(my_entities.size());

        for (unsigned i=0; i<my_entities.size(); ++i)
        {
          //stk::mesh::Entity elem = my_entities[i];
          // if (my_pMesh.isParentElement(elem,false))
          // {
          //   my_markers[i] = percept::DO_NOTHING;
          // }
          // else
          {
            my_markers[i] = *((Int *)my_pMesh.field_data( *my_marker_field, my_entities[i] ));
          }
          //std::cout << "Storing element marker for element " << m_eMesh.identifier(elem) << " = " << my_markers[i] << std::endl;
        }
      }

      Int
      Mimic_Encr_Entity_Marker::get_marker(stk::mesh::Entity entity)
      {
        std::vector<stk::mesh::Entity>::iterator it = std::find(my_entities.begin(), my_entities.end(), entity);
        if (it == my_entities.end())
        {
          return percept::DO_NOTHING;
        }
        else
        {
          ThrowRequire(my_pMesh.is_valid(*it));
          const UInt index = std::distance(my_entities.begin(), it);
          return my_markers[index];
        }
      }

      static void mimic_encr_function_and_element_marker(percept::PerceptMesh & pMesh, double time, ScalarFieldType & function_field, stk::mesh::Field<int> & marker_field)
      {
        std::vector<stk::mesh::Entity> entities;
        stk::mesh::get_entities( *pMesh.get_bulk_data(), stk::topology::NODE_RANK, entities );

        CoordinatesFieldType* coordField = pMesh.get_coordinates_field();

        for (unsigned i=0; i<entities.size(); ++i)
        {
          stk::mesh::Entity node = entities[i];
          double *coord_data = pMesh.field_data(coordField, node);
          double *function_data = pMesh.field_data(&function_field, node);
          *function_data = (2.0*time-10.0-coord_data[0]-coord_data[1])*(2.0*time-10.0-coord_data[0]-coord_data[1]);
        }

        // Now calculate element field marker

        const double refine_lower = 0.0;
        const double refine_upper = 1.0;
        double coarsen_lower = 4.0;
        double coarsen_upper = 1000.0;

        entities.clear();
        stk::mesh::get_entities( *pMesh.get_bulk_data(), stk::topology::ELEMENT_RANK, entities );

        for (unsigned i=0; i<entities.size(); ++i)
        {
          stk::mesh::Entity elem = entities[i];

          //if (pMesh.isParentElement(elem,false)) continue;

          int elem_mark = 0;
          const percept::MyPairIterRelation elem_nodes (pMesh, elem, stk::topology::NODE_RANK );
          unsigned num_node = elem_nodes.size();
          for (unsigned inode=0; inode < num_node; inode++)
          {
            stk::mesh::Entity node = elem_nodes[ inode ].entity();
            double cur_value = *(pMesh.field_data(&function_field, node));

            //If one of the values is in the refine region, immediately mark and return
            if(refine_upper>=cur_value && cur_value>=refine_lower)
            {
              elem_mark = 1;
              break;
            }
            else if(coarsen_lower<=cur_value && cur_value<=coarsen_upper)
            {
              elem_mark = -1;
            }
          }
          Int * marker_data = stk::mesh::field_data( marker_field, elem );
          *marker_data = elem_mark;
        }

        // update ghost element values
        {
          std::vector< const stk::mesh::FieldBase *> fields(1, &marker_field);
          //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
          stk::mesh::communicate_field_data(pMesh.get_bulk_data()->aura_ghosting(), fields);
        }
      }

      TEST(regr_localRefiner, break_tri_to_tri_N_EdgeFromElementMarker_MimicEncr)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1 || p_size == 3)
        {
          PerceptMesh eMesh;
          eMesh.open(input_files_loc+"square_tri3_uns.e");

          Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
          int scalarDimension = 0; // a scalar
          stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, scalarDimension);
          ScalarFieldType * function_field =
            (ScalarFieldType *) eMesh.add_field("s_node", stk::topology::NODE_RANK, scalarDimension);
          stk::mesh::Field<int> & marker_field =  eMesh.get_fem_meta_data()->declare_field< stk::mesh::Field<int> >(stk::topology::ELEMENT_RANK, "marker_field_1");
          stk::io::set_field_role(marker_field, Ioss::Field::TRANSIENT);
          stk::mesh::put_field_on_mesh( marker_field, eMesh.get_fem_meta_data()->universal_part() , nullptr);
          eMesh.commit();
          eMesh.output_active_children_only(true);

          eMesh.save_as(output_files_loc+"mimic_encr_percept_square_adapt_tri_"+post_fix(p_size)+".e");

          Mimic_Encr_Entity_Marker element_marker(eMesh, &marker_field, stk::topology::ELEMENT_RANK);
          Mimic_Encr_Percept_Edge_Adapter breaker(element_marker, eMesh, break_tri_to_tri_N, proc_rank_field);
          breaker.setRemoveOldElements(false);
          breaker.setAlwaysInitializeNodeRegistry(false);

          for (int ipass=0; ipass < 10; ipass++)
          {
            const double time = ipass+1.0;
            mimic_encr_function_and_element_marker(eMesh, time, *function_field, marker_field);
            element_marker.update_markers();
            breaker.doBreak();

            for (int iunref_pass=0; iunref_pass < 1; iunref_pass++)
            {
              std::cout << "P[" << eMesh.get_rank() << "] iunref_pass= " << iunref_pass << std::endl;
              //mimic_encr_function_and_element_marker(eMesh, time, *function_field, marker_field);
              //element_marker.update_markers();
              //ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
              ElementUnrefineCollection elements_to_unref(*eMesh.get_bulk_data());
              breaker.buildUnrefineList(elements_to_unref);

              std::cout << "tmp srk ipass= " << ipass << " elements_to_unref.size() = " << elements_to_unref.size() << std::endl;
              breaker.unrefineTheseElements(elements_to_unref);
            }

            std::stringstream fileid_ss;
            fileid_ss << std::setfill('0') << std::setw(4) << ipass+2;
            eMesh.save_as(output_files_loc+"mimic_encr_percept_square_adapt_tri_"+post_fix(p_size)+".e-s"+fileid_ss.str(), time);
          }
        }
      }

      std::string print_field_type(const stk::mesh::DataTraits                  & arg_traits ,
                                   unsigned                            arg_rank ,
                                   const shards::ArrayDimTag * const * arg_tags )
      {
        //ThrowRequireMsg(arg_rank < 8, "Invalid field rank: " << arg_rank);

        std::ostringstream oss;
        oss << "FieldBase<" ;
        oss << arg_traits.name ;
        for ( unsigned i = 0 ; i < arg_rank ; ++i ) {
          oss << "," << arg_tags[i]->name();
        }
        oss << ">" ;
        return oss.str();
      }

      TEST(regr_localRefiner, break_tri_to_tri_N_Element_From_File)
      {
        EXCEPTWATCH;

        percept::PerceptMesh eMesh;
        eMesh.open(input_files_loc+"square_tri3_out.gold.e");

        Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);

        stk::mesh::FieldBase* proc_rank_field =
          eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK);

        //typedef stk::mesh::Field<double, stk::mesh::Cartesian3d> RFT;
        //typedef stk::mesh::Field<double> RFT;
        typedef stk::mesh::Field<double,shards::ArrayDimension>  RFT;
        //typedef stk::mesh::FieldBase RFT;
        RFT *refine_field =
          eMesh.get_fem_meta_data()->get_field<RFT>(stk::topology::ELEMENT_RANK, "marker_field");
        std::cout << "refine_field= " << refine_field << " " << typeid(*refine_field).name() << std::endl;
        std::cout << "type= " << print_field_type(refine_field->data_traits() ,
                                                  refine_field->field_array_rank() ,
                                                  refine_field->dimension_tags()) << std::endl;
        if (refine_field->type_is<double>()) std::cout << "000 double" << std::endl;
        if (refine_field->type_is<int>()) std::cout << "000 int" << std::endl;

        RefineFieldType *refine_field_int = dynamic_cast<RefineFieldType *>(eMesh.add_field_int("marker_field_int", stk::topology::ELEMENT_RANK, 0));
        (void)refine_field_int;

        eMesh.commit();

        eMesh.read_database_at_step(1); // read first time step
        eMesh.copy_element_field(refine_field_int, refine_field);
        eMesh.save_as("new.e");

        ElementRefinePredicate erp(eMesh, 0, refine_field_int, 0.0);

        PredicateBasedElementAdapter<ElementRefinePredicate>
        breaker(erp, eMesh, break_tri_to_tri_N, proc_rank_field);

        breaker.setRemoveOldElements(false);
        breaker.setAlwaysInitializeNodeRegistry(false);

        breaker.doBreak();

        breaker.deleteParentElements();

        eMesh.save_as(output_files_loc+"local_tri_N_ElementBased.e");
      }

      TEST(regr_localRefiner, break_tri_to_tri_N_Edge_From_Hessian)
      {
        EXCEPTWATCH;

        percept::PerceptMesh eMesh;
        eMesh.open(input_files_loc+"square_tri3_out.gold.e");

        Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);

        stk::mesh::FieldBase* proc_rank_field =
          eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK);

        typedef stk::mesh::Field<double> RFT;
        RFT *refine_field = dynamic_cast<RFT *>(eMesh.add_field("marker_field", stk::topology::ELEMENT_RANK));

        eMesh.commit();

        EdgeBasedInterpError ebp(eMesh, 0, refine_field);

        PredicateBasedEdgeAdapter<EdgeBasedInterpError>
          breaker(ebp, eMesh, break_tri_to_tri_N, proc_rank_field);

        breaker.setRemoveOldElements(false);
        breaker.setAlwaysInitializeNodeRegistry(false);

        breaker.doBreak();
        breaker.doBreak();
        breaker.doBreak();

        breaker.deleteParentElements();

        eMesh.save_as(output_files_loc+"local_tri_N_EdgeBased.e");
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      static void output_time_nodes_elems_error(std::ofstream &outfile,
						double time,
						percept::PerceptMesh eMesh,
						double error)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_rank = stk::parallel_machine_rank( pm );

	int num_nodes = eMesh.get_number_nodes();
	int num_elems = eMesh.get_number_elements();

	if (p_rank==0) {
	  std::cout << "time, num nodes, num elems, error = "
		    << time << " "
		    << num_nodes << " "
		    << num_elems << " "
		    << error << std::endl;

	  outfile << time << " "
		  << num_nodes << " "
		  << num_elems << " "
		  << error << std::endl;
	}
      }

      static std::pair<int,double>
      verify_transition_element_adaptivity(std::string &error_filename,
					   std::string &mesh_name,
					   std::string &function_string,
					   int num_refine,
					   double refine_fraction,
					   double unrefine_fraction = 0.0,
					   int num_steps = 0)
      {
        EXCEPTWATCH;

        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_rank = stk::parallel_machine_rank( pm );

	// I. build mesh and register fields

        percept::PerceptMesh eMesh(3);

        eMesh.open(mesh_name);
        eMesh.output_active_children_only(true);

	const int spatialDimension = eMesh.get_spatial_dim();

	stk::mesh::FieldBase* proc_rank_field =
          eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK);

	RefineFieldType * refineField = &(eMesh.get_fem_meta_data()->declare_field<percept::RefineFieldType>(stk::topology::ELEMENT_RANK, "refine_field"));
	stk::mesh::put_field_on_mesh(*refineField, eMesh.get_fem_meta_data()->universal_part(), nullptr);
	stk::io::set_field_role(*refineField, Ioss::Field::TRANSIENT);

	RefineFieldType * refineFieldOrig = &(eMesh.get_fem_meta_data()->declare_field<percept::RefineFieldType>(stk::topology::ELEMENT_RANK, "refine_field_orig"));
	stk::mesh::put_field_on_mesh(*refineFieldOrig, eMesh.get_fem_meta_data()->universal_part(), nullptr);
	stk::io::set_field_role(*refineFieldOrig, Ioss::Field::TRANSIENT);

        RefineLevelType * refine_level_field =
          &(eMesh.get_fem_meta_data()->declare_field<RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level"));
        stk::mesh::put_field_on_mesh( *refine_level_field ,
                              eMesh.get_fem_meta_data()->universal_part(), nullptr);
	stk::io::set_field_role(*refine_level_field, Ioss::Field::TRANSIENT);

	{
	  TransitionElementType& transition_element       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>((spatialDimension==2 ? stk::topology::ELEMENT_RANK : stk::topology::FACE_RANK), "transition_element");
	  stk::mesh::put_field_on_mesh( transition_element , eMesh.get_fem_meta_data()->universal_part(), nullptr);
	  stk::io::set_field_role(transition_element, Ioss::Field::TRANSIENT);
    if (spatialDimension == 3)
      eMesh.set_transition_element_field_2d(&transition_element);
    else
      eMesh.set_transition_element_field(&transition_element);
	}

	if (spatialDimension==3)
	{
	  TransitionElementType& transition_element_3       = eMesh.get_fem_meta_data()->declare_field<TransitionElementType>(stk::topology::ELEMENT_RANK, "transition_element_3");
	  stk::mesh::put_field_on_mesh( transition_element_3 , eMesh.get_fem_meta_data()->universal_part(), nullptr);
	  stk::io::set_field_role(transition_element_3, Ioss::Field::TRANSIENT);
    eMesh.set_transition_element_field(&transition_element_3);
	}

        stk::mesh::FieldBase *scalar_field =
          eMesh.add_field("scalar", stk::topology::NODE_RANK);

        stk::mesh::FieldBase *error_field =
          eMesh.add_field("error", stk::topology::ELEMENT_RANK);

#define DO_UMR 0

        // these register fields
#if DO_UMR
        //Tri3_Tri3_4 break_uniform(eMesh);
        Tet4_Tet4_8 break_uniform(eMesh);
#else
	UniformRefinerPatternBase *refinementPattern = NULL;
	if (     error_filename.find("tri3") != std::string::npos)
	  refinementPattern = new Local_Tri3_Tri3_N_HangingNode(eMesh);
	else if (error_filename.find("tet4") != std::string::npos)
	  refinementPattern = new Local_Tet4_Tet4_N_HangingNode(eMesh);
	else
	  if (p_rank==0) {
	    std::cout << "Error! could not find refinement pattern." << std::endl;
	    abort();
	  }

#endif

        eMesh.commit();

	// II. set up error norms

        stk::mesh::BulkData* bulkData = eMesh.get_bulk_data();

        StringFunction sf_scalar( function_string.c_str() , Name("sfunc"), spatialDimension, 1);
        FieldFunction ff_scalar("ffunc", scalar_field, eMesh, spatialDimension, 1);

        StringFunction error_func( "sfunc - ffunc" , Name("error"), spatialDimension, 1);

        ConstantFunction result;
        Norm<2> l2Norm(*bulkData);

        l2Norm.setNormField(error_field);

        ff_scalar.interpolateFrom(sf_scalar);

        l2Norm(error_func, result);

	std::ofstream outfile;
	if (p_rank==0) {
	  outfile.open(error_filename.c_str());
	}
	output_time_nodes_elems_error(outfile, 0.0, eMesh, result.getValue());

        //uniform case
#if DO_UMR
        UniformRefiner breaker(eMesh, break_uniform, proc_rank_field);
	num_refine = 4;

        // adaptive case
#else
	stk::mesh::Selector selector = stk::mesh::selectField(*refineField); //(eMesh.get_fem_meta_data()->universal_part());
	ElementRefinePredicate erp(eMesh, &selector, refineField, 0.0);
	TransitionElementAdapter<ElementRefinePredicate>
	  breaker(erp, eMesh, *refinementPattern, proc_rank_field);

#endif

	// III. build marker

        MarkerInfo markerInfo;

        markerInfo.errorIndicator_  = (GenericFieldType *) error_field;
        markerInfo.refineField_     = (RefineFieldType *) refineField;
        markerInfo.refineFieldOrig_ = (RefineFieldType *) refineFieldOrig; // ??
        markerInfo.refineLevelField_ = refine_level_field;
        markerInfo.maxRefinementNumberOfElementsFraction_ = -1.0;
        markerInfo.maxRefinementLevel_ = 100;
        markerInfo.refineFraction_ = refine_fraction;
        markerInfo.unrefineFraction_ = unrefine_fraction;
        markerInfo.useMarker_ = true;

        MarkerUsingErrIndFraction marker(*bulkData, markerInfo);

        breaker.setRemoveOldElements(false);
        breaker.setAlwaysInitializeNodeRegistry(false);

	// run marker once on initial mesh
	marker.mark();

        std::string fileName = output_files_loc+"local_tri_N_ElementInterpError.e";
        eMesh.save_as(fileName,0.0);

	// IV. run adaptive loop

        for (int istep=0; istep <= num_steps; istep++) {

	  // HACK! set time on string func
	  const double time = (num_steps==0) ? 0.0 : (double)istep/(double)num_steps;
	  sf_scalar.set_current_time(time);

	  for (int ipass=1; ipass <= num_refine; ipass++) {

#if DO_UMR
	    breaker.doBreak();
#else
	    breaker.refine();
	    breaker.unrefine();
#endif

	    ff_scalar.interpolateFrom(sf_scalar);

	    l2Norm(error_func, result);

	    output_time_nodes_elems_error(outfile, time, eMesh, result.getValue());

	    marker.mark();

	    std::ostringstream ostr;
	    ostr << std::setfill('0') << std::setw(5) << (istep*num_refine+ipass);

	    const double output_time = (num_steps==0) ? (double)ipass : time;

	    fileName = std::string(output_files_loc+"local_tri_N_ElementInterpError") + std::string(".e-s") + ostr.str();
	    eMesh.save_as(fileName,output_time);

	  }
	}

	breaker.deleteParentElements();
	fileName = std::string(output_files_loc+"local_tri_N_ElementInterpError_final.e");
	eMesh.save_as(fileName);

	delete refinementPattern;

	return std::pair<int,double>(eMesh.get_number_elements(), result.getValue());
      }

      TEST(regr_localRefiner, verify_transition_tri3_smooth_static)
      {
	std::string error_filename = "errors_tri3_smooth_static.dat";
	std::string mesh_name = input_files_loc+"square_tri3_uns.e";
	std::string function_string = "sin(PI*x/5)*sin(PI*y/5)";
	int num_refine = 10;
	double refine_fraction = 0.75;

	//std::pair<int,double> result =
	  verify_transition_element_adaptivity(error_filename, mesh_name, function_string, num_refine, refine_fraction);

	//1249 0.0561288
	//EXPECT_EQ(result.first, 1249);
        //EXPECT_NEAR(result.second, 0.0561288, 1e-5);
      }

      TEST(regr_localRefiner, verify_transition_tet4_smooth_static)
      {
        if (1) return;
        /*
	std::string error_filename = "errors_tet4_smooth_static.dat";
	std::string mesh_name = input_files_loc+"cube_tet4_uns.e";
	std::string function_string = "sin(PI*x/5)*sin(PI*y/5)*sin(PI*z/5)";
	int num_refine = 8;
	double refine_fraction = 0.75;

	std::pair<int,double> result =
	  verify_transition_element_adaptivity(error_filename, mesh_name, function_string, num_refine, refine_fraction);

	EXPECT_EQ(result.first, 5147);
        EXPECT_NEAR(result.second, 1.3353670331057148, 1e-8);
        */
      }

      TEST(regr_localRefiner, verify_transition_tri3_nonsmooth_static)
      {
        if (1) return;
        /*
	std::string error_filename = "errors_tri3_nonsmooth_static.dat";
	std::string mesh_name = input_files_loc+"lshape_tri3.g";
	std::string function_string = "pow(radius,2./3.)*sin((2./3.)*PI*theta)";
	int num_refine = 10;
	double refine_fraction = 0.75;

        StringFunction radius( "sqrt(x*x+y*y)" , Name("radius"), 2, 1);
        StringFunction theta ( "atan2(y,x)" ,    Name("theta"),  2, 1);

	std::pair<int,double> result =
	  verify_transition_element_adaptivity(error_filename, mesh_name, function_string, num_refine, refine_fraction);

  //0 1345 2571 0.00191135

	EXPECT_EQ(result.first, 2571);
        EXPECT_NEAR(result.second, 0.00191135, 1e-5);
        */
      }

      TEST(regr_localRefiner, verify_transition_tri3_smooth_transient)
      {
        if (1) return;
        /*
	std::string error_filename = "errors_tri3_smooth_transient.dat";
	std::string mesh_name = input_files_loc+"square_tri3_uns.e";
	std::string function_string = "sin(PI*(x-5*t)/5)*sin(PI*(y+5*sqrt(3)*t)/5)";
	int num_steps  = 4;
	int num_refine = 2;
	double refine_fraction = 0.30;
	double unrefine_fraction = 0.10;

	// for testing and parameter studies
	if (0) {
	  std::ifstream input_file("verify_transition_tri3_smooth_transient.dat");
	  input_file >> error_filename;
	  input_file >> mesh_name;
	  input_file >> function_string;
	  input_file >> num_steps;
	  input_file >> num_refine;
	  input_file >> refine_fraction;
	  input_file >> unrefine_fraction;
	}

	std::pair<int,double> result =
	  verify_transition_element_adaptivity(error_filename, mesh_name, function_string, num_refine, refine_fraction, unrefine_fraction, num_steps);

  //1 66 108 0.801269

	EXPECT_EQ(result.first, 108);
        EXPECT_NEAR(result.second, 0.801269, 1e-5);
        */
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(regr_localRefiner, test_add_children_to_parts)
      {

        stk::ParallelMachine pm = MPI_COMM_WORLD ;
        const unsigned p_size = stk::parallel_machine_size( pm );
        //const unsigned p_rank = stk::parallel_machine_rank( pm );

        if (p_size <= 3)
          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = true;
            percept::QuadFixture<double, shards::Triangle<3> > fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(-1,1, -1, 1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            Tri3_Tri3_4 break_pattern(eMesh);

            stk::mesh::Part& part = eMesh.get_fem_meta_data()->declare_part("user_defined_part", stk::topology::EDGE_RANK);
            stk::io::put_io_part_attribute(part);

            eMesh.commit();
            fixture.generate_mesh();
            eMesh.save_as(output_files_loc+"test_add_child_0.e");

            stk::mesh::Selector on_locally_owned_part =  ( eMesh.get_fem_meta_data()->locally_owned_part() );

            // add initial elements to part
            {
              eMesh.get_bulk_data()->modification_begin();
              std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part*>(0));
              std::vector<stk::mesh::Part*> remove_parts;


              add_parts[0] = &part;

              std::vector<stk::mesh::Entity> toChange;
              const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( stk::topology::EDGE_RANK );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;
                  if (!on_locally_owned_part(bucket))
                    continue;
                  const unsigned num_entity_in_bucket = bucket.size();
                  for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                    {
                      stk::mesh::Entity element = bucket[ientity];
                      toChange.push_back(element);
                    }
                }
              for (unsigned ii=0; ii < toChange.size(); ++ii)
                {
                  stk::mesh::Entity element = toChange[ii];
                  eMesh.get_bulk_data()->change_entity_parts( element, add_parts, remove_parts );
                }
              stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());
              eMesh.get_bulk_data()->modification_end();
            }
            eMesh.save_as(output_files_loc+"test_add_child_1.e");

            UniformRefiner breaker(eMesh, break_pattern, 0);
            breaker.setRemoveOldElements(false);
            breaker.setAlwaysInitializeNodeRegistry(false);
            breaker.doBreak();
            eMesh.save_as(output_files_loc+"test_add_child_2.e");

            size_t nsl = stk::mesh::count_selected_entities(stk::mesh::Selector(part) & on_locally_owned_part,
                                                            eMesh.get_bulk_data()->buckets(stk::topology::EDGE_RANK));
            size_t g_nsl = 0;
            stk::all_reduce_sum(pm, &nsl, &g_nsl, 1);

            size_t expected = n*4 + 2*n*4;
            EXPECT_EQ(g_nsl,expected);

          }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(regr_localRefiner, block_refiner_quad)
      {
        bool do_test = true;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          bool do_bootstrap_mesh = true;
          if (do_bootstrap_mesh)
          {
            const unsigned n = 5;
            const unsigned nx = n , ny = n;

            bool createEdgeSets = true;
            percept::QuadFixture<double> fixture( pm , nx , ny, createEdgeSets);
            fixture.set_bounding_box(-1,1, -1, 1);

            bool isCommitted = false;
            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data, isCommitted);

            eMesh.commit();

            fixture.generate_mesh();
            eMesh.save_as(input_files_loc+"quad_block_refine_0.e");
            //eMesh.print_info("test1",2);
            //eMesh.dump_vtk("sqquad3.vtk",false);
          }

          {
            BlockRefiner br(input_files_loc+"quad_block_refine_0.e", "quad_block_refine.e", true, true);
            br.refine("block_1",3);
          }
        }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(regr_localRefiner, block_refiner_hex)
      {
        bool do_test = true;
        //stk::ParallelMachine pm = MPI_COMM_WORLD ;

        if (do_test) {

          // if this is a fresh installation, set to true to generate the initial meshes needed for this test (once only)
          {
            BlockRefiner br(input_files_loc+"hex-7block.g", "hex_block_refine.e", false, true);
            br.refine("block_1",3);
          }
        }
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================
      TEST(regr_localRefiner, obtuse_tet)
      {
        EXCEPTWATCH;
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        //const unsigned p_rank = stk::parallel_machine_rank( pm );
        const unsigned p_size = stk::parallel_machine_size( pm );
        if (p_size == 1)
          {
              unsigned npts=4;
              static  SingleTetFixture::Point node_coord_data[  ] = {
                //{ 0.1 , 0.1 , 0.1 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 }, {2, 1, .5} };
                //{ 0.1 , 0.1 , 0.1 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 }, {2, 1, .5} };
              { 0.0 , 0.0 , 0.0 } , { 1 , 0 , 0 } , { 0.5 , .866 , 0 } , { -0.5 , 0.5 , 0.5 }, {2, 1, .5} };

              unsigned ntets=1;
              static  SingleTetFixture::TetIds tetra_node_ids[] = {
                { 1, 2, 3, 4}, {2, 3, 4, 5} };

              percept::SingleTetFixture mesh(pm, false, npts, node_coord_data, ntets, tetra_node_ids);

              stk::io::put_io_part_attribute(  mesh.m_block_tet );

              mesh.m_metaData.commit();
              mesh.populate();

              bool isCommitted = true;
              percept::PerceptMesh eMesh(&mesh.m_metaData, &mesh.m_bulkData, isCommitted);
              eMesh.save_as("obtuse_tet.e");
          }
      }


    } // namespace regression_tests
  } // namespace percept

