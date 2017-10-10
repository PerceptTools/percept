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

#include <Ioss_StructuredBlock.h>

#include <percept/Percept.hpp>
#include <percept/PerceptUtils.hpp>
#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>

#if !STK_PERCEPT_LITE
#include <percept/structured/StructuredGridRefiner.hpp>
#if HAVE_CUBIT
#include <percept/structured/StructuredGridSnap.hpp>
#endif
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#endif

#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>

#if HAVE_CUBIT
#include <PGeom.hpp>
#endif


// Ioss
#define NO_XDMF_SUPPORT
#if defined(STK_BUILT_IN_SIERRA)
#include <main/io_info.h>
//#include <main/struc-to-unstruc.h>
#else
//FIXME tmp
//#include <io_info.h>
#endif

namespace percept
{
  namespace regression_tests
  {

    TEST(regr_PerceptCGNS, test1_unstructured)
    {
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1) return;
      eMesh.open("test1.cgns", "cgns");
      eMesh.commit();
      eMesh.save_as("test1.e");
    }

    TEST(regr_PerceptCGNS, test2_unstructured)
    {
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
      PerceptMesh eMesh(3u);
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

      std::cout << "No.Blocks: "<< bsg->m_sblocks.size() << std::endl;

    for(unsigned iBlock=0;iBlock<bsg->m_sblocks.size();iBlock++) {
        // snap some nodes at z=0 to a cylindrical boundary
        // defined by x*x+z*z=1
        Array4D coord_field_iBlock = *(coord_field->m_block_fields[iBlock]);
        Array4D::HostMirror host_coord_field = Kokkos::create_mirror_view(coord_field_iBlock);
        Kokkos::deep_copy(host_coord_field,coord_field_iBlock);

        std::vector<StructuredGrid::MTNode> block_nodes;
        bsg->get_nodes_of_sb(block_nodes,iBlock);
        for (auto node : block_nodes) {
            double data[3];
        for(unsigned iDim=0;iDim<3;iDim++)
            data[iDim]=host_coord_field(node[0],node[1],node[2],iDim);
            if (std::abs(data[2]) < 1e-5)// if (boundarySelector_5(node))
            {
                const double x = data[0];
                data[2] = -std::sqrt(1.0-x*x);
                for(unsigned iDim=0;iDim<3;iDim++)
                    host_coord_field(node[0],node[1],node[2],iDim) = data[iDim];
            }
        }
        Kokkos::deep_copy(coord_field_iBlock,host_coord_field);
      }//foreach block
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

      int innerIter      = 10;
      
      SGridBoundarySelector boundarySelector (&eMesh);
      ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid> 
        pmmpsi(&eMesh, NULL,&boundarySelector, 0, innerIter, 1.e-1, 1);
      pmmpsi.run();
      
      bsg->dump_vtk(prefix+"_smooth");
    }
    
#if !STK_PERCEPT_LITE
    TEST(regr_PerceptCGNS, test6_sgrid_refine)
    {
      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1) return;

      std::array<unsigned, 3> sizes{{4,4,4}};
      std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);

      bsg->m_sblocks[0]->dump_vtk("cube");

      StructuredGridRefiner sgr1(bsg, 0);

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
      eMesh.readBulkData();

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

    unsigned do_refine_cgns_structured_refinement_only(const std::string& file, int print_level=0, bool allow_parallel=false)
    {
      PerceptMesh eMesh(3u,  MPI_COMM_WORLD);
      if (!allow_parallel && eMesh.get_parallel_size() > 1) return 0;

      std::string dbtype("cgns_structured");
      eMesh.open(file, dbtype);
      eMesh.readBulkData();

      StructuredGridRefiner sgr(eMesh.get_block_structured_grid(), print_level);

      size_t pos = file.find(".cgns");
      VERIFY_OP_ON(pos, !=, std::string::npos, "bad filename: need .cgns extension, file= "+file);
      std::string prefix = file.substr(0,pos);
      if (print_level)
        sgr.m_input->print(std::cout, print_level-1);
      const unsigned num_new_cells = sgr.do_refine();
      if (print_level)
        sgr.m_output->print(std::cout, print_level-1);

      eMesh.close();
      return num_new_cells;
    }

    TEST(performance, refine_generic_rv_cgns)
    {
      const int numRepeats = 8;

      const int print_level=0;
      const bool allow_parallel=true;

      unsigned num_new_cells = 0;
      const unsigned expected_num_new_cells[3]={4194304, 33554432, 268435456};

      for (int j=0;j<numRepeats;j++){
          num_new_cells=do_refine_cgns_structured_refinement_only("grv-aero_SI_v04a_3d-h-nw_seq16_struc.cgns", print_level, allow_parallel);
      }
      EXPECT_EQ(num_new_cells, expected_num_new_cells[0]);

      for (int j=0;j<numRepeats;j++){
          num_new_cells=do_refine_cgns_structured_refinement_only("grv-aero_SI_v04a_3d-h-nw_seq08_struc.cgns", print_level, allow_parallel);
      }
      EXPECT_EQ(num_new_cells, expected_num_new_cells[1]);

      for (int j=0;j<numRepeats;j++){
          num_new_cells=do_refine_cgns_structured_refinement_only("grv-aero_SI_v04a_3d-h-nw_seq04_struc.cgns", print_level, allow_parallel);
      }
      EXPECT_EQ(num_new_cells, expected_num_new_cells[2]);
    }

    void do_refine_cgns_structured_refinement_only_fixture(unsigned cube_dim=4){

        PerceptMesh eMesh(3u);
        if (eMesh.get_parallel_size() > 1) return;

        std::array<unsigned, 3> sizes{{cube_dim,cube_dim,cube_dim}};
        std::shared_ptr<BlockStructuredGrid> bsg = BlockStructuredGrid::fixture_1(eMesh.parallel(), sizes);

        eMesh.set_block_structured_grid(bsg);


        StructuredGridRefiner sgr1(eMesh.get_block_structured_grid(), 0);

        std::string timer_name = "refine_fixture"+std::to_string(cube_dim*cube_dim*cube_dim);
        stk::diag::Timer     my_timer(timer_name, rootTimerStructured());
        stk::diag::TimeBlock my_timeblock(my_timer);
        sgr1.do_refine();
    }

    TEST(DISABLED_kokkosPerformance, single_block_cube){
        unsigned size=0;
        std::cout << "enter your cube dimension\n";
        std::cin>>size;
        do_refine_cgns_structured_refinement_only_fixture(size);
    }


#if HAVE_CUBIT
    void do_refine_cgns_structured_with_geom(const std::string& file, const std::string& geom_file, int print_level=0, bool allow_parallel=false)
    {
      stk::diag::Timer     my_timer(file, rootTimerStructured());
      stk::diag::TimeBlock my_timeblock(my_timer);

      PerceptMesh eMesh(3u);
      if (eMesh.get_parallel_size() > 1 && !allow_parallel) return;

      eMesh.open(file, "cgns_structured");
      eMesh.readBulkData();

      StructuredGridRefiner sgr(eMesh.get_block_structured_grid(), print_level);

      size_t pos = file.find(".cgns");
      VERIFY_OP_ON(pos, !=, std::string::npos, "bad filename: need .cgns extension, file= "+file);
      std::string prefix = file.substr(0,pos);
      if (print_level)
        sgr.m_input->print(std::cout, print_level-1);
      sgr.m_input->dump_vtk(prefix);
      sgr.do_refine();

      // snap the refined mesh to the geometry
      if (geom_file.length() != 0.0)
        {
          std::shared_ptr<PGeom> pgeom(new PGeom());

          if (pgeom->initialize(OPENNURBS_GEOMETRY_ENGINE))
            {
              pgeom->import_open_nurbs_file(geom_file.c_str());
              StructuredGridSnap str_snap(sgr.m_output);
              str_snap.snap_to_geometry(pgeom);
            }
        }

      sgr.m_output->dump_vtk(prefix+"_ref");
      if (print_level)
        sgr.m_output->print(std::cout, print_level-1);
    }
#endif

#if HAVE_CUBIT
    TEST(regr_PerceptCGNS, test_sgrid_refine_and_snap)
    {
#if !KOKKOS_HAVE_CUDA
      do_refine_cgns_structured_with_geom("half_cyl.cgns", "half_cyl.3dm", 0, true);
#endif
    }
#endif

    void do_test_sgrid_zone_connectivity(const std::string file, const int print_level = 0)
    {
        PerceptMesh eMesh(3u,  MPI_COMM_WORLD);
        if (eMesh.get_parallel_size() > 1) return;

        std::string dbtype("cgns_structured");
        eMesh.open(file, dbtype);
        eMesh.readBulkData();

        std::shared_ptr<BlockStructuredGrid> bsg = eMesh.get_block_structured_grid();

        if (print_level) bsg->print(std::cout, print_level);

        for (unsigned iblock=0; iblock < bsg->m_sblocks.size(); ++iblock)
        {
            std::cout << "iblock " << iblock << std::endl;
            unsigned noBad=0;
            std::shared_ptr<StructuredBlock> sblock = bsg->m_sblocks[iblock];
            for (unsigned izone=0; izone < sblock->m_zoneConnectivity.size(); izone++)
            {
                std::cout << "      izone " << izone << std::endl;
                Ioss::ZoneConnectivity zoneConnectivity = sblock->m_zoneConnectivity[izone];
                const unsigned dblockid = zoneConnectivity.m_donorZone - 1;
                std::shared_ptr<StructuredBlock> dblock = bsg->m_sblocks[dblockid];

                Ioss::IJK_t localBeg;
                localBeg[0] = std::min(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localBeg[1] = std::min(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localBeg[2] = std::min(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                Ioss::IJK_t localEnd;
                localEnd[0] = std::max(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localEnd[1] = std::max(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localEnd[2] = std::max(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                Array4D don_a4d = dblock->m_sgrid_coords;
                Array4D loc_a4d = sblock->m_sgrid_coords;

                for (        int I=localBeg[0]; I<=localEnd[0]; I++) {
                    for (    int J=localBeg[1]; J<=localEnd[1]; J++) {
                        for (int K= localBeg[2]; K<=localEnd[2]; K++) {
                            Ioss::IJK_t indexLocal = {{I,J,K}};

                            Kokkos::Array<int, 3> local_indices;
                            Kokkos::Array<int, 3> transform_arr;
                            Kokkos::Array<int, 3> local_range_beg;
                            Kokkos::Array<int, 3> donor_range_beg;

                            for(unsigned iii=0;iii<3;iii++)
                            {
                                local_indices[iii] = indexLocal[iii];
                                transform_arr[iii] = zoneConnectivity.m_transform[iii];
                                local_range_beg[iii] =  zoneConnectivity.m_rangeBeg[iii];
                                donor_range_beg[iii] = zoneConnectivity.m_donorRangeBeg[iii];
                            }

                            Kokkos::Array<int,3 > index_donor = device_safe_transform_block_indices(local_indices,
                                                            transform_arr,
                                                            local_range_beg,
                                                            donor_range_beg);

                            for (int d=0; d<3; d++) {
                                                            //from based 1 to based 0 indexing
                                const Double v = loc_a4d(indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,d);
                                const Double vt = don_a4d(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,d);

                                if( (device_safe_abs_flt(v-vt)>1.e-14) ){
                                    noBad++;
                                    break;
                                }
                            }
                        }
                    }
                }



                Kokkos::Array<int,3> local_sizes;
                for(unsigned ijk=0;ijk<3;ijk++)
                    local_sizes[ijk] = 1 + localEnd[ijk]-localBeg[ijk];

                Kokkos::Array<int, 3> transform_arr;
                Kokkos::Array<int, 3> local_range_beg;
                Kokkos::Array<int, 3> donor_range_beg;

                for(unsigned iii=0;iii<3;iii++)
                {
                    transform_arr[iii] = zoneConnectivity.m_transform[iii];
                    local_range_beg[iii] =  zoneConnectivity.m_rangeBeg[iii];
                    donor_range_beg[iii] = zoneConnectivity.m_donorRangeBeg[iii];
                }

                for(unsigned iNode=0;iNode<(unsigned)( localEnd[0]*localEnd[1]*localEnd[2] );iNode++)
                {
                    Kokkos::Array<unsigned, 3> indx;

                    indx[0] = (localBeg[0] -1)
                            + ( iNode % local_sizes[0] );
                    indx[1] = (localBeg[1] -1)
                            + (iNode / local_sizes[0]) % local_sizes[1];
                    indx[2] = (localBeg[2] -1)
                            +((iNode / local_sizes[0])/local_sizes[1]) % local_sizes[2];

                    Ioss::IJK_t indexLocal = {{(int)(indx[0] +1),
                                                (int)(indx[1] +1),
                                                (int)(indx[2] +1)}};

                    Kokkos::Array<int, 3> local_indices;
                    for(unsigned iii=0;iii<3;iii++)
                        local_indices[iii] = indexLocal[iii];

                    Kokkos::Array<int,3 > index_donor = device_safe_transform_block_indices(local_indices,
                            transform_arr,
                            local_range_beg,
                            donor_range_beg);

                    for (int d=0; d<3; d++) {
                        //from based 1 to based 0 indexing
                        const Double v = loc_a4d(indexLocal[0]-1,indexLocal[1]-1,indexLocal[2]-1,d);
                        const Double vt = don_a4d(index_donor[0]-1,index_donor[1]-1,index_donor[2]-1,d);

                        if( (device_safe_abs_flt(v-vt)>1.e-14) ){
                            noBad++;
                            break;
                        }
                    }
                }//foreach node
                EXPECT_EQ(noBad,unsigned(0));
            }//izone
        }//iblock
    }//do_test_sgrid_zone_connectivity

    TEST(regr_PerceptCGNS, test_sgrid_zone_connectivity)
    {
        do_test_sgrid_zone_connectivity("cube_4blocks.cgns");
        do_test_sgrid_zone_connectivity("cube_2blocks.cgns");
        do_test_sgrid_zone_connectivity("5blocks.cgns");
        do_test_sgrid_zone_connectivity("blunt_wedge.cgns");
        do_test_sgrid_zone_connectivity("bulge.cgns");
        do_test_sgrid_zone_connectivity("delta_vertex.cgns");
        //do_test_sgrid_zone_connectivity("dlr-f6.coar.cgns", 2); // test fails
        do_test_sgrid_zone_connectivity("half_cyl.cgns");
        do_test_sgrid_zone_connectivity("rv.cgns");
        do_test_sgrid_zone_connectivity("sqnz_s.hdf.cgns");

        Kokkos::View<double*,DataLayout,MemSpace> v1("v1",3);
        v1(0)=.1;
        v1(1)=.2;
        v1(2)=.3;
        auto subv1 = subview(v1,Kokkos::make_pair(3,0));
        auto subv2 = subview(v1,Kokkos::make_pair(0,3));
        std::cout << subv1.size() << " " << subv2.size()<< std::endl;
    }

    void do_test_sgrid_zone_connectivity_field_data_comm(const std::string file, const int print_level = 0)
    {
        PerceptMesh eMesh(3u,  MPI_COMM_WORLD);
        if (eMesh.get_parallel_size() > 1) return;

        std::string dbtype("cgns_structured");
        eMesh.open(file, dbtype);
        eMesh.readBulkData();

        std::shared_ptr<BlockStructuredGrid> bsg = eMesh.get_block_structured_grid();

        if (print_level) bsg->print(std::cout, print_level);

        std::cout <<file << std::endl;

        for (unsigned iblock=0; iblock < bsg->m_sblocks.size(); ++iblock)
        {   //SUM SHARED FEILDS
            std::shared_ptr<StructuredBlock> sblock = bsg->m_sblocks[iblock];
            for (unsigned izone=0; izone < sblock->m_zoneConnectivity.size(); izone++)
            {
                Ioss::ZoneConnectivity zoneConnectivity = sblock->m_zoneConnectivity[izone];
                const unsigned dblockid = zoneConnectivity.m_donorZone - 1;
                std::shared_ptr<StructuredBlock> dblock = bsg->m_sblocks[dblockid];

                Ioss::IJK_t localBeg;
                localBeg[0] = std::min(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localBeg[1] = std::min(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localBeg[2] = std::min(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                Ioss::IJK_t localEnd;
                localEnd[0] = std::max(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localEnd[1] = std::max(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localEnd[2] = std::max(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                auto loc_coords = subview(sblock->m_sgrid_coords,
                        Kokkos::make_pair(localBeg[0]-1,localEnd[0]),
                        Kokkos::make_pair(localBeg[1]-1,localEnd[1]),
                        Kokkos::make_pair(localBeg[2]-1,localEnd[2]),
                        Kokkos::make_pair(0,3));

//                Ioss::IJK_t donorBeg;// = zoneConnectivity.m_donorRangeBeg;
//                Ioss::IJK_t donorEnd;// = zoneConnectivity.m_donorRangeEnd;
                Ioss::IJK_t donorBeg;
                donorBeg[0] = std::min(zoneConnectivity.m_donorRangeBeg[0],zoneConnectivity.m_donorRangeEnd[0]);
                donorBeg[1] = std::min(zoneConnectivity.m_donorRangeBeg[1],zoneConnectivity.m_donorRangeEnd[1]);
                donorBeg[2] = std::min(zoneConnectivity.m_donorRangeBeg[2],zoneConnectivity.m_donorRangeEnd[2]);

                Ioss::IJK_t donorEnd;
                donorEnd[0] = std::max(zoneConnectivity.m_donorRangeBeg[0],zoneConnectivity.m_donorRangeEnd[0]);
                donorEnd[1] = std::max(zoneConnectivity.m_donorRangeBeg[1],zoneConnectivity.m_donorRangeEnd[1]);
                donorEnd[2] = std::max(zoneConnectivity.m_donorRangeBeg[2],zoneConnectivity.m_donorRangeEnd[2]);


                auto donor_coords = subview(dblock->m_sgrid_coords,
                        Kokkos::make_pair(donorBeg[0]-1,donorEnd[0]),
                        Kokkos::make_pair(donorBeg[1]-1,donorEnd[1]),
                        Kokkos::make_pair(donorBeg[2]-1,donorEnd[2]),
                        Kokkos::make_pair(0,3));

                if(donor_coords.size()!=loc_coords.size()){
                    std::cout << "  SUM: donor and loc sizes don't match on block " << iblock << " and zone " << izone
                              << " sizes according to subviews DONOR :" << donor_coords.size() << " LOCAL :" << loc_coords.size()
                              << " sizes according to connectivity limits DONOR: " << (donorEnd[0]-donorBeg[0])*(donorEnd[1]-donorBeg[1])*(donorEnd[2]-donorBeg[2])
                              << " LOCAL: " << (1+zoneConnectivity.m_rangeEnd[0] - zoneConnectivity.m_rangeBeg[0]) * (1+zoneConnectivity.m_rangeEnd[1] - zoneConnectivity.m_rangeBeg[1]) * (1+zoneConnectivity.m_rangeEnd[2] - zoneConnectivity.m_rangeBeg[2])
                              <<std::endl;
                    std::cout << "DONOR limits END " << donorEnd[0] << " " << donorEnd[1] << " " << donorEnd[2] << std::endl;
                    std::cout << "DONOR limits BEG " << donorBeg[0] << " " << donorBeg[1] << " " << donorBeg[2] << std::endl;

                    std::cout << "LOCAL limits END " << zoneConnectivity.m_rangeEnd[0] << " " << zoneConnectivity.m_rangeEnd[1] << " " << zoneConnectivity.m_rangeEnd[2] << std::endl;
                    std::cout << "LOCAL limits BEG " << zoneConnectivity.m_rangeBeg[0] << " " << zoneConnectivity.m_rangeBeg[1] << " " << zoneConnectivity.m_rangeBeg[2] << std::endl;
                }

                if(iblock<dblockid){
                    for (        int I=0; I<localEnd[0] - localBeg[0]; I++) {
                        for (    int J=0; J<localEnd[1] - localBeg[1]; J++) {
                            for (int K=0; K<localEnd[2] - localBeg[2]; K++) {
                                Ioss::IJK_t indexLocal = {{I,J,K}};

                                for (int d=0; d<3; d++)
                                    loc_coords(indexLocal[0],indexLocal[1],indexLocal[2],d) += donor_coords(indexLocal[0],indexLocal[1],indexLocal[2],d);
                            }//K
                        }//J
                    }//I
                }
            }//foreach zone
        }//foreach block

        for (unsigned iblock=0; iblock < bsg->m_sblocks.size(); ++iblock)
        {   //PROPOGATE FIELDS
            std::shared_ptr<StructuredBlock> sblock = bsg->m_sblocks[iblock];
            for (unsigned izone=0; izone < sblock->m_zoneConnectivity.size(); izone++)
            {
                Ioss::ZoneConnectivity zoneConnectivity = sblock->m_zoneConnectivity[izone];
                const unsigned dblockid = zoneConnectivity.m_donorZone - 1;
                std::shared_ptr<StructuredBlock> dblock = bsg->m_sblocks[dblockid];

                Ioss::IJK_t localBeg;
                localBeg[0] = std::min(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localBeg[1] = std::min(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localBeg[2] = std::min(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                Ioss::IJK_t localEnd;
                localEnd[0] = std::max(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localEnd[1] = std::max(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localEnd[2] = std::max(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                auto loc_coords = subview(sblock->m_sgrid_coords,
                        Kokkos::make_pair(localBeg[0]-1,localEnd[0]),
                        Kokkos::make_pair(localBeg[1]-1,localEnd[1]),
                        Kokkos::make_pair(localBeg[2]-1,localEnd[2]),
                        Kokkos::make_pair(0,3));

//                Ioss::IJK_t donorBeg;// = zoneConnectivity.m_donorRangeBeg;
//                Ioss::IJK_t donorEnd;// = zoneConnectivity.m_donorRangeEnd;
                Ioss::IJK_t donorBeg;
                donorBeg[0] = std::min(zoneConnectivity.m_donorRangeBeg[0],zoneConnectivity.m_donorRangeEnd[0]);
                donorBeg[1] = std::min(zoneConnectivity.m_donorRangeBeg[1],zoneConnectivity.m_donorRangeEnd[1]);
                donorBeg[2] = std::min(zoneConnectivity.m_donorRangeBeg[2],zoneConnectivity.m_donorRangeEnd[2]);

                Ioss::IJK_t donorEnd;
                donorEnd[0] = std::max(zoneConnectivity.m_donorRangeBeg[0],zoneConnectivity.m_donorRangeEnd[0]);
                donorEnd[1] = std::max(zoneConnectivity.m_donorRangeBeg[1],zoneConnectivity.m_donorRangeEnd[1]);
                donorEnd[2] = std::max(zoneConnectivity.m_donorRangeBeg[2],zoneConnectivity.m_donorRangeEnd[2]);


                auto donor_coords = subview(dblock->m_sgrid_coords,
                        Kokkos::make_pair(donorBeg[0]-1,donorEnd[0]),
                        Kokkos::make_pair(donorBeg[1]-1,donorEnd[1]),
                        Kokkos::make_pair(donorBeg[2]-1,donorEnd[2]),
                        Kokkos::make_pair(0,3));

                if(donor_coords.size()!=loc_coords.size()){
                    std::cout << "  PROP: donor and loc sizes don't match on block " << iblock << " and zone " << izone
                              << " sizes according to subviews DONOR :" << donor_coords.size() << " LOCAL :" << loc_coords.size()
                              << " sizes according to connectivity limits DONOR: " << (donorEnd[0]-donorBeg[0])*(donorEnd[1]-donorBeg[1])*(donorEnd[2]-donorBeg[2])
                              << " LOCAL: " << (1+zoneConnectivity.m_rangeEnd[0] - zoneConnectivity.m_rangeBeg[0]) * (1+zoneConnectivity.m_rangeEnd[1] - zoneConnectivity.m_rangeBeg[1]) * (1+zoneConnectivity.m_rangeEnd[2] - zoneConnectivity.m_rangeBeg[2])
                              <<std::endl;
                    std::cout << "DONOR limits END " << donorEnd[0] << " " << donorEnd[1] << " " << donorEnd[2] << std::endl;
                    std::cout << "DONOR limits BEG " << donorBeg[0] << " " << donorBeg[1] << " " << donorBeg[2] << std::endl;

                    std::cout << "LOCAL limits END " << zoneConnectivity.m_rangeEnd[0] << " " << zoneConnectivity.m_rangeEnd[1] << " " << zoneConnectivity.m_rangeEnd[2] << std::endl;
                    std::cout << "LOCAL limits BEG " << zoneConnectivity.m_rangeBeg[0] << " " << zoneConnectivity.m_rangeBeg[1] << " " << zoneConnectivity.m_rangeBeg[2] << std::endl;
                }

                if(iblock>dblockid){
                    for (        int I=0; I<localEnd[0] - localBeg[0]; I++) {
                        for (    int J=0; J<localEnd[1] - localBeg[1]; J++) {
                            for (int K=0; K<localEnd[2] - localBeg[2]; K++) {
                                Ioss::IJK_t indexLocal = {{I,J,K}};

                                for (int d=0; d<3; d++)
                                    loc_coords(indexLocal[0],indexLocal[1],indexLocal[2],d) = donor_coords(indexLocal[0],indexLocal[1],indexLocal[2],d);
                            }//K
                        }//J
                    }//I
                }
            }//foreach zone
        }//foreach block

        for (unsigned iblock=0; iblock < bsg->m_sblocks.size(); ++iblock)
        {   //VERIFY FIELDS SUMMED CORRECTLY
            unsigned noBad=0;
            std::shared_ptr<StructuredBlock> sblock = bsg->m_sblocks[iblock];
            for (unsigned izone=0; izone < sblock->m_zoneConnectivity.size(); izone++)
            {
                Ioss::ZoneConnectivity zoneConnectivity = sblock->m_zoneConnectivity[izone];
                const unsigned dblockid = zoneConnectivity.m_donorZone - 1;
                std::shared_ptr<StructuredBlock> dblock = bsg->m_sblocks[dblockid];

                Ioss::IJK_t localBeg;
                localBeg[0] = std::min(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localBeg[1] = std::min(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localBeg[2] = std::min(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                Ioss::IJK_t localEnd;
                localEnd[0] = std::max(zoneConnectivity.m_rangeBeg[0],zoneConnectivity.m_rangeEnd[0]);
                localEnd[1] = std::max(zoneConnectivity.m_rangeBeg[1],zoneConnectivity.m_rangeEnd[1]);
                localEnd[2] = std::max(zoneConnectivity.m_rangeBeg[2],zoneConnectivity.m_rangeEnd[2]);

                auto loc_coords = subview(sblock->m_sgrid_coords,
                        Kokkos::make_pair(localBeg[0]-1,localEnd[0]),
                        Kokkos::make_pair(localBeg[1]-1,localEnd[1]),
                        Kokkos::make_pair(localBeg[2]-1,localEnd[2]),
                        Kokkos::make_pair(0,3));

                //                Ioss::IJK_t donorBeg;// = zoneConnectivity.m_donorRangeBeg;
                //                Ioss::IJK_t donorEnd;// = zoneConnectivity.m_donorRangeEnd;
                Ioss::IJK_t donorBeg;
                donorBeg[0] = std::min(zoneConnectivity.m_donorRangeBeg[0],zoneConnectivity.m_donorRangeEnd[0]);
                donorBeg[1] = std::min(zoneConnectivity.m_donorRangeBeg[1],zoneConnectivity.m_donorRangeEnd[1]);
                donorBeg[2] = std::min(zoneConnectivity.m_donorRangeBeg[2],zoneConnectivity.m_donorRangeEnd[2]);

                Ioss::IJK_t donorEnd;
                donorEnd[0] = std::max(zoneConnectivity.m_donorRangeBeg[0],zoneConnectivity.m_donorRangeEnd[0]);
                donorEnd[1] = std::max(zoneConnectivity.m_donorRangeBeg[1],zoneConnectivity.m_donorRangeEnd[1]);
                donorEnd[2] = std::max(zoneConnectivity.m_donorRangeBeg[2],zoneConnectivity.m_donorRangeEnd[2]);


                auto donor_coords = subview(dblock->m_sgrid_coords,
                        Kokkos::make_pair(donorBeg[0]-1,donorEnd[0]),
                        Kokkos::make_pair(donorBeg[1]-1,donorEnd[1]),
                        Kokkos::make_pair(donorBeg[2]-1,donorEnd[2]),
                        Kokkos::make_pair(0,3));

                if(donor_coords.size()!=loc_coords.size()){
                    std::cout << "  CHECK: donor and loc sizes don't match on block " << iblock << " and zone " << izone
                            << " sizes according to subviews DONOR :" << donor_coords.size() << " LOCAL :" << loc_coords.size()
                            << " sizes according to connectivity limits DONOR: " << (donorEnd[0]-donorBeg[0])*(donorEnd[1]-donorBeg[1])*(donorEnd[2]-donorBeg[2])
                            << " LOCAL: " << (1+zoneConnectivity.m_rangeEnd[0] - zoneConnectivity.m_rangeBeg[0]) * (1+zoneConnectivity.m_rangeEnd[1] - zoneConnectivity.m_rangeBeg[1]) * (1+zoneConnectivity.m_rangeEnd[2] - zoneConnectivity.m_rangeBeg[2])
                            <<std::endl;
                    std::cout << "DONOR limits END " << donorEnd[0] << " " << donorEnd[1] << " " << donorEnd[2] << std::endl;
                    std::cout << "DONOR limits BEG " << donorBeg[0] << " " << donorBeg[1] << " " << donorBeg[2] << std::endl;

                    std::cout << "LOCAL limits END " << zoneConnectivity.m_rangeEnd[0] << " " << zoneConnectivity.m_rangeEnd[1] << " " << zoneConnectivity.m_rangeEnd[2] << std::endl;
                    std::cout << "LOCAL limits BEG " << zoneConnectivity.m_rangeBeg[0] << " " << zoneConnectivity.m_rangeBeg[1] << " " << zoneConnectivity.m_rangeBeg[2] << std::endl;
                }

                if(iblock>dblockid){
                    for (        int I=0; I<localEnd[0] - localBeg[0]; I++) {
                        for (    int J=0; J<localEnd[1] - localBeg[1]; J++) {
                            for (int K=0; K<localEnd[2] - localBeg[2]; K++) {
                                Ioss::IJK_t indexLocal = {{I,J,K}};
                                for (int d=0; d<3; d++) {
                                    const double v = loc_coords(indexLocal[0],indexLocal[1],indexLocal[2],d);
                                    const double vt = donor_coords(indexLocal[0],indexLocal[1],indexLocal[2],d);

                                    EXPECT_NEAR(v,vt,1.e-14);
                                    if(device_safe_abs_flt(v-vt)<1.e-14)
                                        noBad++;
                                }
                            }
                        }
                    }
                    EXPECT_EQ(noBad,unsigned(0));
                }
            }
        }
    }//do_test_sgrid_zone_connectivity_field_data_comm

    TEST(regr_PerceptCGNS, test_sgrid_zone_connectivity_field_data_comm)
    {
        do_test_sgrid_zone_connectivity_field_data_comm("cube_4blocks.cgns");
        do_test_sgrid_zone_connectivity_field_data_comm("cube_2blocks.cgns");
        do_test_sgrid_zone_connectivity_field_data_comm("5blocks.cgns");
        do_test_sgrid_zone_connectivity_field_data_comm("blunt_wedge.cgns");
        do_test_sgrid_zone_connectivity_field_data_comm("bulge.cgns");
        do_test_sgrid_zone_connectivity_field_data_comm("delta_vertex.cgns");
        //do_test_sgrid_zone_connectivity_field_data_comm("dlr-f6.coar.cgns", 2); // test fails
        do_test_sgrid_zone_connectivity_field_data_comm("half_cyl.cgns");
        do_test_sgrid_zone_connectivity_field_data_comm("rv.cgns");
        do_test_sgrid_zone_connectivity_field_data_comm("sqnz_s.hdf.cgns");
    }

#endif // STK_PERCEPT_LITE
#endif // STK_BUILT_IN_SIERRA

  } // namespace regression_tests
} // namespace percept

