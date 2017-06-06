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
#include <percept/PerceptUtils.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#if !STK_PERCEPT_LITE
#include <percept/structured/StructuredGridRefiner.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#endif

#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>

// Ioss
#define NO_XDMF_SUPPORT
#if defined(STK_BUILT_IN_SIERRA)
#include <main/io_info.h>
//#include <main/struc-to-unstruc.h>
#else
//FIXME tmp
//#include <io_info.h>
#endif

//#include "RegressionTestLocalRefiner.hpp"

namespace percept
{
  namespace regression_tests
  {

    TEST(regr_PerceptCGNS, test1_unstructured)
    {
      bool do_test = true;
      if (!do_test) return;
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1) return;
      eMesh.open("test1.cgns", "cgns");
      eMesh.commit();
      eMesh.save_as("test1.e");
    }

    TEST(regr_PerceptCGNS, test2_unstructured)
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
    TEST(regr_PerceptCGNS, test3_unstructured)
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
    TEST(regr_PerceptCGNS, test4_in_memory_write_test_unstructured)
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

#if 0
    TEST(regr_PerceptCGNS, test5_in_memory_write_test_unstructured)
    {
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 2) return;

      bool in_memory_write = false;
      std::pair<Ioss::DatabaseIO *, Ioss::DatabaseIO *> ret_db = Ioss::StrucToUnstruc::file_copy("sqnz_s.hdf.cgns", "a.e", "percept", "1.0", in_memory_write);
      (void)ret_db;
      //delete ret_db.first;
      //delete ret_db.second;

      if (!in_memory_write)
        ret_db.second->closeDatabase();

      Ioss::PropertyManager properties;
      if (in_memory_write)
        properties.add(Ioss::Property(std::string("MEMORY_READ"), 1));
      Ioss::DatabaseIO *    dbi = Ioss::IOFactory::create("exodus", "a.e", Ioss::READ_MODEL,
                                                    (MPI_Comm)MPI_COMM_WORLD, properties);
      if (dbi == nullptr || !dbi->ok(true)) {
        std::exit(EXIT_FAILURE);
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      Ioss::Region region(dbi, "region_1");

      eMesh.set_ioss_read_options("auto-decomp:yes");
      //eMesh.open_in_memory("a.e", ret_db.second->get_region());
      if (in_memory_write)
        eMesh.open_in_memory("a.e", &region);
      else
        eMesh.open("a.e");


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
    
    void do_snap_to_cylinder(std::shared_ptr<BlockStructuredGrid> bsg, const std::string& prefix) 
    {
      PerceptMesh eMesh(3u);

      eMesh.set_block_structured_grid(bsg);

      StructuredGrid::MTField *coord_field     = bsg->m_fields["coordinates"].get();
      StructuredGrid::MTField *coord_field_N   = bsg->m_fields["coordinates_N"].get();

      // snap some nodes at z=0 to a cylindrical boundary 
      // defined by x*x+z*z=1
      std::vector<StructuredGrid::MTNode> nodes;
      bsg->get_nodes(nodes);
      for (auto node : nodes){
        double data[3];
        get_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
        if (std::abs(data[2]) < 1e-5) // if (boundarySelector_5(node))
          {
            const double x = data[0];
            data[2] = -std::sqrt(1.0-x*x);
            set_field<StructuredGrid>(data, 3, &eMesh, coord_field, node);
          }
      }

      // save state of projected mesh
      // field, dst, src:
      bsg->copy_field(coord_field_N, coord_field);
      
      bsg->dump_vtk(prefix+"_snap");
    }

    void do_smooth_block_structured(std::shared_ptr<BlockStructuredGrid> bsg, const std::string& prefix)
    {
      PerceptMesh eMesh(3u);

      bsg->register_field("coordinates", 3);
      
      eMesh.set_block_structured_grid(bsg);
      eMesh.add_coordinate_state_fields(); // for smoothing
      
      StructuredGrid::MTField *coord_field     = bsg->m_fields["coordinates"].get();
      StructuredGrid::MTField *coord_field_NM1 = bsg->m_fields["coordinates_NM1"].get();
        
      // save state of original mesh
      // field, dst, src:
      bsg->copy_field(coord_field_NM1, coord_field);

      do_snap_to_cylinder(bsg, prefix);

      int  msq_debug     = 0; // 1,2,3 for more debug info
      bool always_smooth = true;
      int innerIter      = 10;
      
      CubeBoundarySelector boundarySelector (&eMesh);
      ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> 
        pmmpsi(&eMesh, &boundarySelector, 0, innerIter, 1.e-1, 1);
      pmmpsi.run(always_smooth, msq_debug);
      
      bsg->dump_vtk(prefix+"_smooth");
    }
    
#if !STK_PERCEPT_LITE
    TEST(regr_PerceptCGNS, test6_sgrid_refine)
    {
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1) return;

      std::array<unsigned, 3> sizes{{4,4,4}};
      std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);

      eMesh.set_block_structured_grid(bsg);

      bsg->m_sblocks[0]->dump_vtk("cube");

      StructuredGridRefiner sgr1(eMesh.get_block_structured_grid(), 0);

      sgr1.m_input->dump_vtk("cube_init");

      sgr1.do_refine();

      sgr1.m_output->dump_vtk("cube_ref");
      
      do_smooth_block_structured(sgr1.m_output, "cube_ref");
    }

    void do_refine_cgns_structured(const std::string& file, int print_level=0, bool allow_parallel=false)
    {
      stk::diag::Timer     my_timer(file, rootTimerStructured());
      stk::diag::TimeBlock my_timeblock(my_timer);

      PerceptMesh eMesh(3u);
      if (!allow_parallel && eMesh.get_parallel_size() > 1) return;

      eMesh.open(file, "cgns_structured");
      eMesh.commit();

      StructuredGridRefiner sgr(eMesh.get_block_structured_grid(), print_level);

      size_t pos = file.find(".cgns");
      VERIFY_OP_ON(pos, !=, std::string::npos, "bad filename: need .cgns extension, file= "+file);
      std::string prefix = file.substr(0,pos);
      if (print_level)
        sgr.m_input->print(std::cout, print_level-1);
      sgr.m_input->dump_vtk(prefix);
      sgr.do_refine();
      sgr.m_output->dump_vtk(prefix+"_ref");
      if (print_level)
        sgr.m_output->print(std::cout, print_level-1);
    }

    TEST(regr_PerceptCGNS, test7_sgrid_refine)
    {
      do_refine_cgns_structured("sqnz_s.hdf.cgns", 0, true);
    }

    TEST(regr_PerceptCGNS, test8_sgrid_refine)
    {
      do_refine_cgns_structured("5blocks.cgns");
      do_refine_cgns_structured("delta_vertex.cgns");
      do_refine_cgns_structured("dlr-f6.coar.cgns");
    }

    TEST(regr_PerceptCGNS, test9_sgrid_refine_rv)
    {
      do_refine_cgns_structured("rv.cgns", 2);
    }

    TEST(regr_PerceptCGNS, test10_sgrid_volume)
    {
#ifdef KOKKOS_HAVE_CUDA
      return;
#endif

      int print_level = 0;
      (void)print_level;
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1) return;

      bool in_memory_write = false;
      (void) in_memory_write;
      std::array<unsigned, 3> sizes{{2,2,2}};
      std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);
      eMesh.set_block_structured_grid(bsg);
      bsg->m_sblocks[0]->dump_vtk("cube_vol");
      bsg->dump_vtk("bsg_cube");

      StructuredCellIndex element{{0,0,0,0}};

      typename StructuredGrid::MTField *coord_field = bsg->m_fields["coordinates"].get();
      VERIFY_OP_ON(coord_field, !=, 0, "bad coord_field");

      JacobianUtilImpl<StructuredGrid> jacA;
      double jac=0.0;
      jacA(jac, eMesh, element, coord_field);
      std::cout << "jac= " << jac << std::endl;
      EXPECT_NEAR(jac, 1.0, 1.e-6);
      //std::cout << "output_xyz=\n" << sgr1.m_output->m_sblocks[0]->m_sgrid_coords << std::endl;
      if (1)
        {
          Math::Matrix rmx = Math::rotationMatrix(0, -60.0);
          Math::Matrix rmy = Math::rotationMatrix(1, 0.0);
          Math::Matrix rmz = Math::rotationMatrix(2, -30.0);

          Math::Matrix rm;
          rm = ublas::prod(rmx, rmz);
          rm = ublas::prod(rmy, rm);
          //rm = ublas::prod(scm, rm);

          using Node = typename StructuredGrid::MTNode;
          std::vector<Node> nodes;
          bsg->get_nodes(nodes);
          for (auto node : nodes)
            {
              Math::Vector v;
              get_field<StructuredGrid>(v.data(), 3, &eMesh, coord_field, node);
              v = rm * v;
              set_field<StructuredGrid>(v.data(), 3, &eMesh, coord_field, node);
            }
          bsg->dump_vtk("bsg_cube_xformed");
          jacA(jac, eMesh, element, coord_field);
          std::cout << "jac= " << jac << std::endl;
          EXPECT_NEAR(jac, 1.0, 1.e-6);
        }
    }

    void do_refine_cgns_structured_refinement_only(const std::string& file, int print_level=0, bool allow_parallel=false)
    {
      PerceptMesh eMesh(3u,  MPI_COMM_WORLD);
      if (!allow_parallel && eMesh.get_parallel_size() > 1) return;

      eMesh.open(file, "cgns_structured");
      eMesh.commit();

      StructuredGridRefiner sgr(eMesh.get_block_structured_grid(), print_level);

      size_t pos = file.find(".cgns");
      VERIFY_OP_ON(pos, !=, std::string::npos, "bad filename: need .cgns extension, file= "+file);
      std::string prefix = file.substr(0,pos);
      if (print_level)
        sgr.m_input->print(std::cout, print_level-1);
      sgr.do_refine();
      if (print_level)
        sgr.m_output->print(std::cout, print_level-1);

      eMesh.close();
    }

    TEST(DISABLED_kokkosPerformance, evaluation)
    {
      int numRepeats = 20;
      
      for (int j=0;j<numRepeats;j++){
        do_refine_cgns_structured_refinement_only("grv-aero_SI_v04a_3d-h-nw_seq16_struc.cgns");
      }
      
      for (int j=0;j<numRepeats;j++){
        do_refine_cgns_structured_refinement_only("grv-aero_SI_v04a_3d-h-nw_seq08_struc.cgns");
      }

      for (int j=0;j<numRepeats;j++){
        do_refine_cgns_structured_refinement_only("grv-aero_SI_v04a_3d-h-nw_seq04_struc.cgns");
      }
    }

#endif // STK_PERCEPT_LITE
#endif // STK_BUILT_IN_SIERRA

  } // namespace regression_tests
} // namespace percept

