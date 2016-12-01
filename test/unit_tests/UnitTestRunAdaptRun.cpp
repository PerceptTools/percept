// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/Percept.hpp>

#include <stk_mesh/base/FieldBLAS.hpp>

#include <gtest/gtest.h>
#include <percept/PerceptMesh.hpp>
#include <percept/util/Loops.hpp>

#include <unit_tests/UnitTestSupport.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/TetWedgeFixture.hpp>

#include <adapt/SerializeNodeRegistry.hpp>
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

    class SetErrorField : public percept::ElementOp
    {
      percept::PerceptMesh& m_eMesh;
      stk::mesh::FieldBase *m_field;
    public:
      SetErrorField(percept::PerceptMesh& eMesh, stk::mesh::FieldBase *field) : m_eMesh(eMesh), m_field(field) {
      }

      virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)
      {
        double *f_data = stk::mesh::field_data(*dynamic_cast<ErrorFieldType *>(field), element);
        stk::mesh::FieldBase* coord_field = m_eMesh.get_coordinates_field();

        const MyPairIterRelation elem_nodes(bulkData, element, stk::topology::NODE_RANK );

        unsigned num_node = elem_nodes.size();
        double c[] = {0,0,0};
        for (unsigned inode=0; inode < num_node; inode++)
          {
            stk::mesh::Entity node = elem_nodes[ inode ].entity();
            double *c_data = m_eMesh.field_data(coord_field, node);
            c[0] += c_data[0]/double(num_node);
            c[1] += c_data[1]/double(num_node);
            //c[2] += c_data[2]/double(num_node);
          }
        double edge_len = m_eMesh.edge_length_ave(element);
        f_data[0] = 2*c[0]*c[1]*edge_len*edge_len;
        return false;  // don't terminate the loop
      }
      virtual void init_elementOp()
      {
        std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
        stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
      }
      virtual void fini_elementOp() {
        std::vector< const stk::mesh::FieldBase *> fields(1,m_field);
        stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
      }

    };

    /// test parsing, including file: option
    TEST(adapt_unit_rar, run_adapt_run)
    {
      stk::ParallelMachine pm = MPI_COMM_WORLD ;

      {
        const unsigned n = 12;
        //const unsigned nx = n , ny = n , nz = p_size*n ;
        const unsigned nx = n , ny = n;

        percept::QuadFixture<double> fixture( pm , nx , ny, true);
        fixture.set_bounding_box(0,1,0,1);
        fixture.meta_data.commit();
        fixture.generate_mesh();

        percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
        eMesh.save_as("rar_input.e");
        eMesh.close();
      }

      {
        percept::PerceptMesh eMesh(2u);
        eMesh.open("rar_input.e");
        eMesh.register_and_set_refine_fields();

        ErrorFieldType * error_field = &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
        stk::mesh::put_field( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1);
        stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);

        (void) error_field;
        eMesh.commit();

        std::string input = "{error_indicator_field: myError, "
          "marker: {type: fraction, refine_fraction: 0.2, unrefine_fraction: 0.2}, "
          "max_number_elements_fraction: 0.2, max_refinement_level: 3, extra_output: yes"
          "}";

        std::string output = "rar.yaml";

        bool debug=true;
        percept::RunAdaptRunInfo rar(eMesh, input, output, debug);
        rar.create();
        rar.create_marker();
      }

      {
        percept::PerceptMesh eMesh(2u);
        eMesh.open("rar_input.e");
        eMesh.register_and_set_refine_fields();
        ErrorFieldType * error_field = &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
        stk::mesh::put_field( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1);
        stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);
        eMesh.commit();

        bool debug=true;
        std::string input = "{file: rar_input.yaml}";
        std::string output = "rar1.yaml";
        percept::RunAdaptRunInfo rar(eMesh, input, output, debug);
        rar.create();
        rar.create_marker();

        SetErrorField set_err_field(eMesh, error_field);
        elementOpLoop(*eMesh.get_bulk_data(), set_err_field, error_field);
        eMesh.save_as("rar_err.e");

      }

    }

    /// create a sequence of meshes to adapt outside this routine
    static void generate_run_adapt_run_seq(std::string type="quad")
    {
      stk::ParallelMachine pm = MPI_COMM_WORLD ;

      {
        const unsigned n = 12;
        //const unsigned nx = n , ny = n , nz = p_size*n ;
        const unsigned nx = n , ny = n;

        if (type == "tri")
          {
            percept::QuadFixture<double, shards::Triangle<3>  > fixture( pm , nx , ny, true);
            fixture.set_bounding_box(0,1,0,1);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.save_as("rar_input.e");
            eMesh.close();
          }
        else
          {
            percept::QuadFixture<double> fixture( pm , nx , ny, true);
            fixture.set_bounding_box(0,1,0,1);
            fixture.meta_data.commit();
            fixture.generate_mesh();

            percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
            eMesh.save_as("rar_input.e");
            eMesh.close();
          }
      }

      {
        percept::PerceptMesh eMesh(2u);
        eMesh.open("rar_input.e");
        eMesh.register_and_set_refine_fields();
        //stk::mesh::FieldBase *error_field = eMesh.add_field("myError", stk::topology::ELEMENT_RANK);
        ErrorFieldType * error_field = &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
        stk::mesh::put_field( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1);
        stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);

        UniformRefinerPatternBase *localBreakPattern = 0;
        if (type == "quad")
          localBreakPattern = new Local_Quad4_Quad4_N_Transition(eMesh);
        else
          localBreakPattern = new Local_Tri3_Tri3_N_HangingNode(eMesh);

        eMesh.commit();
        eMesh.output_active_children_only(false);

        std::string input = "{family_tree_extension: ft, error_indicator_field: myError, "
          "marker: {type: fraction, refine_fraction: 0.2, unrefine_fraction: 0.2}, "
          "max_number_elements_fraction: 0.2, max_refinement_level: 3, extra_output: yes"
          // following are optional and have defaults
          ", refine_field: refine_field, refine_level_field: refine_level,"
          "transition_element_field: transition_element, parent_element_field: parent_element,"
          "new_nodes_field: new_nodes"
          "}";

        std::string output = "rar.yaml";

        bool debug = true;
        percept::RunAdaptRunInfo rar(eMesh, input, output, debug);
        rar.create();
        rar.create_marker();

        SetErrorField set_err_field(eMesh, error_field);
        elementOpLoop(*eMesh.get_bulk_data(), set_err_field, error_field);
        eMesh.save_as("rar_err.e");

        stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

        ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);
        TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, *localBreakPattern, 0);

        // special for local adaptivity
        breaker.setRemoveOldElements(false);
        breaker.setAlwaysInitializeNodeRegistry(false);


        //stk::mesh::Entity element = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), 121);
        int nref=3;
        for (int iref=0; iref < nref; ++iref)
          {

            elementOpLoop(*eMesh.get_bulk_data(), set_err_field, error_field);
            rar.m_marker->mark();

            bool enforce_what[3] = {false, true, false};
            breaker.refine( enforce_what);
            std::cout << "tet_transition number elements after pass " << iref << " = " << eMesh.get_number_elements() << std::endl;
            if (1)
              {
                char buf[1000];
                sprintf(buf, "%04d", iref);
                if (iref == 0)
                  eMesh.save_as(type+"-rar-anim.e");
                else
                  eMesh.save_as(type+"-rar-anim.e-s"+std::string(buf));
              }
          }

        delete localBreakPattern;
      }

    }

    TEST(adapt_unit_rar, run_adapt_run_seq)
    {
      generate_run_adapt_run_seq();
    }

    TEST(adapt_unit_rar, run_adapt_run_seq_1)
    {
      std::string type = "tri";
      //std::string type = "quad";
      generate_run_adapt_run_seq(type);

      EXCEPTWATCH;

      stk::ParallelMachine pm = MPI_COMM_WORLD ;
      bool debug = false;

      const unsigned p_rank = stk::parallel_machine_rank( pm );
      (void)p_rank;
      const unsigned p_size = stk::parallel_machine_size( pm );
      if (p_size <= 3)
        {
          percept::PerceptMesh eMesh(2u);
          eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);

          eMesh.open(type+"-rar-anim.e");

          eMesh.register_and_set_refine_fields();
          eMesh.add_registered_refine_fields_as_input_fields();
          ErrorFieldType * error_field = &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
          stk::mesh::put_field( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1);
          stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);
          UniformRefinerPatternBase *localBreakPattern = 0;
          if (type == "quad")
            localBreakPattern = new Local_Quad4_Quad4_N_Transition(eMesh);
          else
            localBreakPattern = new Local_Tri3_Tri3_N_HangingNode(eMesh);

          eMesh.commit();
          eMesh.output_active_children_only(false);

          stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

          ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);
          TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, *localBreakPattern, 0);

          breaker.setRemoveOldElements(false);
          breaker.setAlwaysInitializeNodeRegistry(false);

          RefinerUtil::rebuild_family_tree(eMesh, debug);
          AdaptedMeshVerifier::check_parent_element_field(eMesh, "after rebuild", debug);

          breaker.initializeRefine();
          erp.setMarkNone(true);
          breaker.initializeDB();
          erp.setMarkNone(false);

          //SetOfEntities emptySet(*eMesh.get_bulk_data());
          //breaker.replaceNodeRegistryOwnership(emptySet, eMesh.element_rank());

          AdaptedMeshVerifier adaptedMeshVerifier(debug);
          if (!adaptedMeshVerifier.isValid(eMesh, true))
            throw std::runtime_error("Invalid initial mesh");

          std::string input = "{family_tree_extension: ft, error_indicator_field: myError, "
            "marker: {type: fraction, refine_fraction: 0.2, unrefine_fraction: 0.1}, "
            "max_number_elements_fraction: 4.0, max_refinement_level: 3, extra_output: yes"
            // following are optional and have defaults
            ", refine_field: refine_field, refine_level_field: refine_level,"
            "transition_element_field: transition_element, parent_element_field: parent_element,"
            "new_nodes_field: new_nodes"
            "}";

          std::string output = "rar.yaml";

          percept::RunAdaptRunInfo rar(eMesh, input, output, debug);
          rar.create();
          rar.create_marker();

          int nref=3;
          for (int iref=0; iref < nref; ++iref)
            {
              SetErrorField set_err_field(eMesh, error_field);
              elementOpLoop(*eMesh.get_bulk_data(), set_err_field, error_field);
              rar.m_marker->mark();

              if (iref == 0)
                eMesh.save_as(type+"-rar-anim-ft.e");

              bool enforce_what[3] = {false, true, false};
              breaker.refine( enforce_what);
              std::cout << "tet_transition number elements after pass " << iref << " = " << eMesh.get_number_elements() << std::endl;
              if (1)
                {
                  char buf[1000];
                  sprintf(buf, "%04d", iref+1);
                  eMesh.save_as(type+"-rar-anim-ft.e-s"+std::string(buf));
                }
            }

          delete localBreakPattern;

        }
    }

    // this code mimics an application running and generating an error indicator field
    double compute_error_on_mesh(const std::string & input, const std::string & output)
    {
      percept::PerceptMesh eMesh(2u);
      eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);
      eMesh.open(input);

      ErrorFieldType * error_field =
        &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
      stk::mesh::put_field( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1);
      stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);

      eMesh.commit();

      SetErrorField set_err_field(eMesh, error_field);
      elementOpLoop(*eMesh.get_bulk_data(), set_err_field, error_field);

      const int num_nodes = eMesh.get_number_nodes();
      const int num_elems = eMesh.get_number_elements();
      const double error_norm = stk::mesh::field_nrm2(*error_field, eMesh.get_fem_meta_data()->universal_part(), MPI_COMM_WORLD)/sqrt((double)num_elems);

      const unsigned p_rank = stk::parallel_machine_rank( MPI_COMM_WORLD );
      if (p_rank==0) {
        std::cout << num_nodes << " "
                  << num_elems << " "
                  << error_norm
                  << std::endl;
      }

      eMesh.save_as(output);

      return error_norm;
    }

    void adapt_mesh(const std::string & input_ft, const std::string & input,
                    const std::string & output_ft, const std::string & output)
    {
      const bool has_ft_mesh = (input_ft != "");

      percept::PerceptMesh eMesh(2u);
      eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);
      eMesh.open(has_ft_mesh ? input_ft : input);
      eMesh.register_and_set_refine_fields();

      ErrorFieldType * error_field =
        &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
      stk::mesh::put_field( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1);
      stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);

      if (has_ft_mesh) {
        eMesh.add_registered_refine_fields_as_input_fields();
      }
      else {
        eMesh.add_input_field(error_field);
      }

      Teuchos::RCP<UniformRefinerPatternBase> localBreakPattern = make_local_break_pattern(eMesh);

      eMesh.commit();

      if (has_ft_mesh) {
        // second mesh with only error field
        percept::PerceptMesh eMesh_error(2u);
        eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);
        eMesh_error.open(input);

        ErrorFieldType * from_error_field =
          &eMesh_error.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
        stk::mesh::put_field( *from_error_field , eMesh_error.get_fem_meta_data()->universal_part(), 1);
        eMesh_error.add_input_field(from_error_field);

        eMesh_error.commit();

        eMesh_error.read_database_at_step(eMesh_error.get_database_time_step_count());

        copy_error_indicator(eMesh_error,eMesh,from_error_field,error_field);
      }

      stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());

      ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);
      TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, *localBreakPattern, 0);

      // special for local adaptivity
      breaker.setRemoveOldElements(false);
      breaker.setAlwaysInitializeNodeRegistry(false);

      bool debug = false;

      RefinerUtil::rebuild_family_tree(eMesh, debug);
      AdaptedMeshVerifier::check_parent_element_field(eMesh, "after rebuild", debug);

      breaker.initializeRefine();

      AdaptedMeshVerifier adaptedMeshVerifier(debug);
      if (!adaptedMeshVerifier.isValid(eMesh, true))
        throw std::runtime_error("Invalid initial mesh");

      erp.setMarkNone(true);
      breaker.initializeDB();
      erp.setMarkNone(false);

      if (1)
        {
          NodeRegistry nr(eMesh);
          RefinerUtil::rebuild_node_registry(eMesh, nr);
          YAML::Emitter yaml0, yaml1;
          SerializeNodeRegistry::serialize_write(nr, yaml1, 0);
          SerializeNodeRegistry::serialize_write(breaker.getNodeRegistry(), yaml0, 0);
          static int iter = 0;
          std::ofstream y0("y0_"+toString(iter)+".yaml");
          std::ofstream y1("y1_"+toString(iter)+".yaml");
          ++iter;
          y0 << "# " << eMesh.getProperty("in_filename") << std::endl;
          y1 << "# " << eMesh.getProperty("in_filename") << std::endl;

          y0 << yaml0.c_str() << std::endl;
          y1 << yaml1.c_str() << std::endl;

        }
      std::string rar_input = "{family_tree_extension: ft, error_indicator_field: myError, "
        "marker: {type: fraction, refine_fraction: 0.2, unrefine_fraction: 0.2}, "
        "max_number_elements_fraction: 0.2, max_refinement_level: 3, extra_output: yes"
        // following are optional and have defaults
        ", refine_field: refine_field, refine_level_field: refine_level,"
        "transition_element_field: transition_element, parent_element_field: parent_element,"
        "new_nodes_field: new_nodes"
        "}";

      std::string rar_output = "rar.yaml";

      percept::RunAdaptRunInfo rar(eMesh, rar_input, rar_output, debug);
      rar.create();
      rar.create_marker();

      rar.m_marker->mark();

      bool enforce_what[3] = {false, true, false};
      breaker.refine( enforce_what);

      eMesh.output_active_children_only(false);
      eMesh.save_as(output_ft);

      // TODO output no fields, just mesh
      eMesh.output_active_children_only(true);
      eMesh.save_as(output);
    }

    TEST(adapt_unit_rar, run_adapt_tri_seq)
    {
      stk::ParallelMachine pm = MPI_COMM_WORLD ;
      {
        const unsigned n = 3;
        const unsigned nx = n , ny = n;

        percept::QuadFixture<double, shards::Triangle<3>  > fixture( pm , nx , ny, true);
        fixture.set_bounding_box(0,1,0,1);
        fixture.meta_data.commit();
        fixture.generate_mesh();

        percept::PerceptMesh eMesh(&fixture.meta_data, &fixture.bulk_data);
        eMesh.save_as("mesh_tri3.g");
        eMesh.close();
      }

      compute_error_on_mesh("mesh_tri3.g", "error0_tri3.e");

      // first pass - only original mesh w/error
      // first arg for mesh with ft is empty
      adapt_mesh("", "error0_tri3.e",
                 "adapt1_ft_tri3.e", "adapt1_tri3.g");

      compute_error_on_mesh("adapt1_tri3.g", "error1_tri3.e");

      // additional passes - now mesh w/error + mesh w/parents and fields
      adapt_mesh("adapt1_ft_tri3.e", "error1_tri3.e",
                 "adapt2_ft_tri3.e", "adapt2_tri3.g");

      compute_error_on_mesh("adapt2_tri3.g", "error2_tri3.e");

      adapt_mesh("adapt2_ft_tri3.e", "error2_tri3.e",
                 "adapt3_ft_tri3.e", "adapt3_tri3.g");

      const double error = compute_error_on_mesh("adapt3_tri3.g", "error3_tri3.e");
      EXPECT_NEAR(error, 0.002083842387, 1e-5);
    }

    TEST(adapt_unit_rar, run_adapt_tet_wedge_seq)
    {
      stk::ParallelMachine pm = MPI_COMM_WORLD ;
      {
        percept::TetWedgeFixture fixture(pm, false, true);

        fixture.m_metaData.commit();
        fixture.populate();

        const bool isCommitted = true;
        percept::PerceptMesh eMesh(&fixture.m_metaData, &fixture.m_bulkData, isCommitted);
        eMesh.save_as("mesh_tet4_wedge6.g");
        eMesh.close();
      }

      compute_error_on_mesh("mesh_tet4_wedge6.g", "error0_tet4_wedge6.e");

      // first pass - only original mesh w/error
      // first arg for mesh with ft is empty
      adapt_mesh("", "error0_tet4_wedge6.e",
                 "adapt1_ft_tet4_wedge6.e", "adapt1_tet4_wedge6.g");

      compute_error_on_mesh("adapt1_tet4_wedge6.g", "error1_tet4_wedge6.e");

      // additional passes - now mesh w/error + mesh w/parents and fields
      adapt_mesh("adapt1_ft_tet4_wedge6.e", "error1_tet4_wedge6.e",
                 "adapt2_ft_tet4_wedge6.e", "adapt2_tet4_wedge6.g");

      compute_error_on_mesh("adapt2_tet4_wedge6.g", "error2_tet4_wedge6.e");

      adapt_mesh("adapt2_ft_tet4_wedge6.e", "error2_tet4_wedge6.e",
                 "adapt3_ft_tet4_wedge6.e", "adapt3_tet4_wedge6.g");

      const double error = compute_error_on_mesh("adapt3_tet4_wedge6.g", "error3_tet4_wedge6.e");
      EXPECT_NEAR(error, 0.0554850220463, 1e-5);
    }

    TEST(adapt_unit_rar, valgrind_error)
    {
      percept::PerceptMesh eMesh(2u);
      eMesh.set_avoid_add_all_mesh_fields_as_input_fields(true);
      eMesh.open("error0_tet4_wedge6.e");
      eMesh.register_and_set_refine_fields();

      ErrorFieldType * error_field =
        &eMesh.get_fem_meta_data()->declare_field<ErrorFieldType>(stk::topology::ELEMENT_RANK, "myError");
      stk::mesh::put_field( *error_field , eMesh.get_fem_meta_data()->universal_part(), 1);
      stk::io::set_field_role( *error_field, Ioss::Field::TRANSIENT);

      eMesh.add_input_field(error_field);

      Teuchos::RCP<UniformRefinerPatternBase> localBreakPattern = make_local_break_pattern(eMesh);

      eMesh.commit();
    }

  } // namespace unit_tests
} // namespace percept


