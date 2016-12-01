// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#include <math.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <gtest/gtest.h>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <stk_io/IossBridge.hpp>

#include <percept/Percept.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>

// Ioss
#define NO_XDMF_SUPPORT
#if defined(STK_BUILT_IN_SIERRA)
#include <main/io_info.h>
#else
//FIXME tmp
//#include <io_info.h>
#endif

//#include "RegressionTestLocalRefiner.hpp"

namespace percept
{
  namespace regression_tests
  {

    TEST(regr_PerceptCGNS, test1)
    {
      bool do_test = true;
      if (!do_test) return;
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1) return;
      eMesh.open("test1.cgns", "cgns");
      eMesh.commit();
      eMesh.save_as("test1.e");
    }

    TEST(regr_PerceptCGNS, test2)
    {
      bool do_test = true;
      if (!do_test) return;
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1) return;

      eMesh.open("test1.cgns", "cgns");

      URP_Heterogeneous_3D break_pattern(eMesh);
      eMesh.commit();

      UniformRefiner breaker(eMesh, break_pattern, 0);
      breaker.doBreak();

      eMesh.save_as("test2-ref.e");

    }

    /// Parallel
    TEST(regr_PerceptCGNS, test3)
    {
      bool do_test = true;
      if (!do_test) return;
      PerceptMesh eMesh(3u);
      //if (eMesh.get_parallel_size() > 1) return;
      eMesh.set_ioss_read_options("auto-decomp:yes");
      eMesh.open("test1.cgns", "cgns");

      URP_Heterogeneous_3D break_pattern(eMesh);
      eMesh.commit();

      UniformRefiner breaker(eMesh, break_pattern, 0);
      breaker.doBreak();

      eMesh.save_as("test3-ref.e");
    }

    // ================================================================================
    /// These functions are just demonstrations of traversing the Ioss data structures,
    ///   taken from stk_io and modified.
    ///
    void declare_ioss_field(const Ioss::Field &io_field,
                            bool use_cartesian_for_scalar)
    {
      std::string name = io_field.get_name();
      std::string field_type = io_field.transformed_storage()->name();
      size_t num_components = io_field.transformed_storage()->component_count();
      std::cout << "in-mem-test: declare_ioss_field io_field.name= " << name << " type= " << field_type << " num_components= " << num_components << std::endl;
    }

    void define_io_fields(Ioss::GroupingEntity *entity,
                          Ioss::Field::RoleType role)

    {
      Ioss::NameList names;
      entity->field_describe(role, &names);

      bool use_cartesian_for_scalar = false;
      if (role == Ioss::Field::ATTRIBUTE)
        use_cartesian_for_scalar = true;

      int ii=0;
      for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I, ++ii) {
        Ioss::Field io_field = entity->get_field(*I);
        (void)io_field;
        std::cout << "in-mem-test: define_io_fields names[" << ii << "] = " << *I << std::endl;
        declare_ioss_field(io_field, use_cartesian_for_scalar);
      }
    }


    void process_nodeblocks(Ioss::Region &region)
    {
      const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
      assert(node_blocks.size() == 1);

      Ioss::NodeBlock *nb = node_blocks[0];
      define_io_fields(nb, Ioss::Field::ATTRIBUTE);
    }

    size_t process_nodeblocks_bulk(Ioss::Region &region, int parallel_size)
    {
      const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
      assert(node_blocks.size() == 1);

      Ioss::NodeBlock *nb = node_blocks[0];

      std::vector<int> ids;
      nb->get_field_data("ids", ids);
      std::cout << "in-mem-test: process_nodeblocks_bulk(bulk): ids.size()= " << ids.size() << std::endl;

      std::cout << "in-mem-test: process_nodeblocks(bulk): ids = " << std::endl;
      size_t freq=10;
      for (size_t i=0; i < ids.size(); i++) {
        std::cout << " " << ids[i];
        if ((i+1) % freq == 0) std::cout << std::endl;
      }
      std::cout << std::endl;

      // Register node sharing information for all nodes on processor
      // boundaries.
      //
      if (parallel_size > 1)
        {
          Ioss::CommSet* io_cs = region.get_commset("commset_node");
          if (io_cs)
            {
              size_t num_sharings = io_cs->get_field("entity_processor").raw_count();
              std::cout << "in-mem-test: num_sharings= " << num_sharings << std::endl;

              size_t global_node_count = region.get_property("global_node_count").get_int();
              (void) global_node_count;

              std::vector<int> entity_proc;
              io_cs->get_field_data("entity_processor", entity_proc);

              for (size_t i = 0; i < num_sharings; ++i) {
                std::cout << "in-mem-test: process_nodeblocks(bulk): id, proc= " << entity_proc[i*2] << " " << entity_proc[i*2+1] << std::endl;
              }
            }
        }
      return ids.size();
    }

#if defined(STK_BUILT_IN_SIERRA)
    TEST(regr_PerceptCGNS, test4_in_memory_write_test)
    {
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 2) return;

      eMesh.set_ioss_read_options("auto-decomp:yes");
      eMesh.open("test1.cgns", "cgns");

      URP_Heterogeneous_3D break_pattern(eMesh);
      eMesh.commit();
      eMesh.save_as("test4.e");

      UniformRefiner breaker(eMesh, break_pattern, 0);
      breaker.doBreak();
      eMesh.save_as("test4-ref.e");

      eMesh.set_ioss_write_options("in-memory");
      eMesh.save_as("test4-ref-in-memory.e");

      Teuchos::RCP<stk::io::StkMeshIoBroker>  io_broker = eMesh.get_ioss_mesh_data_output();
      Teuchos::RCP<Ioss::Region> output_region = io_broker->get_output_io_region(eMesh.get_output_file_index());

      process_nodeblocks(*output_region.get());
      size_t ids_size = process_nodeblocks_bulk(*output_region.get(),eMesh.get_parallel_size());
      std::cout << eMesh.rank() << " ids_size= " << ids_size << std::endl;
      int64_t total_elements = 0;
      if (1)
        {
          Ioss::ElementBlockContainer ebs            = output_region->get_element_blocks();
          for (auto eb : ebs) {
            int64_t num_elem = eb->get_property("entity_count").get_int();
            total_elements += num_elem;
            std::string type       = eb->get_property("topology_type").get_string();
            std::cout << eMesh.rank() << "num elems= " << num_elem << " type= " << type << std::endl;
          }
        }
      if (eMesh.get_parallel_size() == 1)
        {
          EXPECT_EQ(ids_size, 195u);
          EXPECT_EQ(total_elements, 96);
        }
      else
        {
          EXPECT_EQ(ids_size, 105u);
          EXPECT_EQ(total_elements, 48);
        }


      Info::Interface intfc;
      std::vector<char *> argv;
      argv.push_back( const_cast<char *>("progname") );
      //argv.push_back( const_cast<char *>("-summary") );
      argv.push_back( const_cast<char *>("dummy_file.e") );
      int argc = argv.size();
      intfc.parse_options(argc, &argv[0]);
      std::cout << "intfc.summary()= " << intfc.summary() << std::endl;

      Ioss::io_info_set_db_properties(intfc, output_region->get_database());
      Ioss::io_info_file_info(intfc, *output_region);
    }
#endif

  } // namespace regression_tests
} // namespace percept

