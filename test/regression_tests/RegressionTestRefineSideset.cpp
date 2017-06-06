
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

#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/UniformRefinerPattern_Tri3_Quad4_3.hpp>
#include <adapt/UniformRefinerPattern_Tet4_Hex8_4.hpp>
#include <adapt/UniformRefinerPattern_Wedge6_Hex8_6.hpp>
#include <unit_tests/TestLocalRefiner.hpp>
#include <adapt/sierra_element/StdMeshObjTopologies.hpp>
#include <adapt/sierra_element/RefinementTopology.hpp>
#include <percept/RunEnvironment.hpp>

#include <percept/fixtures/Fixture.hpp>
#include <percept/fixtures/BeamFixture.hpp>
#include <percept/fixtures/SingleTetFixture.hpp>
#include <percept/fixtures/TetWedgeFixture.hpp>
#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/fixtures/PyramidFixture.hpp>
#include <percept/fixtures/QuadFixture.hpp>
#include <percept/fixtures/WedgeFixture.hpp>

// smoothing tests
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>

#include <adapt/UniformRefinerPattern_def.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/ReadWriteSidesetTester.hpp>


  namespace percept
  {
    namespace regression_tests
    {
    using namespace stk::unit_test_util::sideset;
    typedef stk::mesh::EntityVector PseudoEntity;

#include "RegressionTestFileLoc.hpp"

    const percept::Elem::RefinementTopology* get_refinement_topology(const CellTopologyData * const cellTopoData)
    {
        shards::CellTopology cellTopo(cellTopoData);
        percept::Elem::CellTopology elem_celltopo = percept::Elem::getBasicCellTopology(cellTopo.getName());
        const percept::Elem::RefinementTopology* ref_topo_p = Elem::getRefinementTopology(elem_celltopo);
        return ref_topo_p;
    }

    const percept::Elem::RefinementTopology* get_refinement_topology(stk::topology topology)
    {
        return get_refinement_topology(stk::mesh::get_cell_topology(topology).getCellTopologyData());
    }

    const percept::Elem::RefinementTopology* get_refinement_topology(percept::PerceptMesh &eMesh, stk::mesh::Entity element)
    {
        return get_refinement_topology(eMesh.get_cell_topology(element));
    }

    class MeshWithSideset {
    public:
        MeshWithSideset(const std::string &inputBaseName) :
            m_baseName(inputBaseName),
            m_meta(3, get_entity_rank_names()),
            m_bulk(m_meta, MPI_COMM_WORLD),
            m_eMesh(&m_meta, &m_bulk, false),
            m_isMeshLoaded(false)
        {

        }

        ~MeshWithSideset()
        {
            delete m_breakPattern;
        }

    private:
        std::vector<std::string> get_entity_rank_names()
        {
            std::vector<std::string> entityRankNames = stk::mesh::entity_rank_names();
            entityRankNames.push_back("FAMILY_TREE");
            return entityRankNames;
        }

        MeshWithSideset();
        MeshWithSideset(const MeshWithSideset&);

    private:
        std::string m_baseName;
        stk::mesh::MetaData m_meta;
        stk::unit_test_util::sideset::BulkDataTester m_bulk;
        percept::PerceptMesh m_eMesh;
        URP_Heterogeneous_3D *m_breakPattern = nullptr;
        bool m_isMeshLoaded;

    public:
        stk::unit_test_util::sideset::BulkDataTester &get_bulk() { return m_bulk; }
        stk::mesh::MetaData &get_meta() { return m_meta; }
        percept::PerceptMesh &get_mesh() { return m_eMesh; }
        UniformRefinerPatternBase &get_breaker() { return *m_breakPattern; }

        std::string get_input_filename()
        {
            std::string inputFilename(input_files_loc + m_baseName + ".e");
            return inputFilename;
        }

        std::string get_output_filename()
        {
            std::string outputFilename(output_files_loc + m_baseName + "_1.e");
            return outputFilename;
        }

        void load_mesh_and_fill_sideset_data()
        {
            if(m_isMeshLoaded) return;

            stk::unit_test_util::sideset::StkMeshIoBrokerTester inputIoBroker;

            inputIoBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
            inputIoBroker.set_rank_name_vector(get_entity_rank_names());
            inputIoBroker.set_bulk_data(m_bulk);

            inputIoBroker.add_mesh_database(get_input_filename(), stk::io::READ_MESH);
            inputIoBroker.create_input_mesh();
            inputIoBroker.add_all_mesh_fields_as_input_fields();

            m_breakPattern = new URP_Heterogeneous_3D(m_eMesh);

            m_eMesh.add_field("proc_rank", stk::topology::ELEMENT_RANK, 0);
            inputIoBroker.populate_bulk_data();

            m_eMesh.setCoordinatesField();

            m_isMeshLoaded = true;
        }

        void write_exodus_file()
        {
            ThrowRequire(m_isMeshLoaded);
            stk::unit_test_util::sideset::StkMeshIoBrokerTester outputIoBroker;
            outputIoBroker.set_bulk_data(m_bulk);
            size_t resultFileIndex = outputIoBroker.create_output_mesh(get_output_filename(), stk::io::WRITE_RESULTS);
            outputIoBroker.write_output_mesh(resultFileIndex);
        }

        void do_uniform_refinement(bool removeOldElements=true)
        {
            ThrowRequire(m_isMeshLoaded);
            UniformRefiner breaker(m_eMesh, *m_breakPattern, m_eMesh.get_field(stk::topology::ELEMENT_RANK, "proc_rank"));
            breaker.setRemoveOldElements(removeOldElements);
            //breaker.setIgnoreSideSets(true);
            breaker.doBreak();
        }
    };

      void refine_mesh(const std::string &baseName,
                       const SideSetIdAndElemIdSidesVector& expectedInputSideset,
                       const SideSetIdAndElemIdSidesVector& expectedRefinedSideset)
      {
        MeshWithSideset inputMesh(baseName);
        inputMesh.load_mesh_and_fill_sideset_data();

        stk::unit_test_util::sideset::BulkDataTester &bulk = inputMesh.get_bulk();
        stk::mesh::MetaData &meta = inputMesh.get_meta();

        compare_sidesets(inputMesh.get_input_filename(), bulk, expectedInputSideset);

        inputMesh.do_uniform_refinement();

        for(const SideSetIdAndElemIdSides& sset : expectedRefinedSideset)
	{
	    stk::mesh::SideSet sideSet = get_stk_side_set(bulk, sset.sideSet);
            bulk.create_sideset(sset.id) = sideSet;
	}

        stk::mesh::EntityVector faces;
        stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(meta.side_rank()), faces);

        for(stk::mesh::Entity face : faces)
            ASSERT_EQ(2u, bulk.num_elements(face));

        inputMesh.write_exodus_file();

        stk::mesh::MetaData meta2;
        stk::unit_test_util::sideset::BulkDataTester bulk2(meta2, MPI_COMM_WORLD);
        stk::unit_test_util::sideset::read_exo_file(bulk2, inputMesh.get_output_filename(), stk::unit_test_util::sideset::READ_SERIAL_AND_DECOMPOSE);
        EXPECT_NO_THROW(compare_sidesets(inputMesh.get_output_filename(), bulk2, expectedRefinedSideset));
      }

      struct TestCase
      {
          std::string filename;
          SideSetIdAndElemIdSidesVector originalSideset;
          SideSetIdAndElemIdSidesVector refinedSideset;
      };

      void add_test_case(const std::string &baseName,
                         const SideSetIdAndElemIdSidesVector &original,
                         const SideSetIdAndElemIdSidesVector &refined,
                         std::vector<TestCase> &testCases)
      {
          // This function is necessary to get around intel compiler issues with seagull brace initialization
          testCases.push_back({baseName, original, refined});
      }

      TEST(regr_stk_sideset_refine, DISABLED_internal_sidesets)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD;

        const unsigned p_size = stk::parallel_machine_size(pm);

        std::vector<TestCase> testCases;
        add_test_case("ADA",
                      {{1, {{ 1, 5}, { 2, 4}                                                      }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}, {27, 4}, {28, 4}, {29, 4}, {30, 4}}}},
                      testCases);
        add_test_case("ARA",
                      {{1, {{ 2, 4}                           }}},
                      {{1, {{27, 4}, {28, 4}, {29, 4}, {30, 4}}}},
                      testCases);
        add_test_case("ALA",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCases);
        add_test_case("ADe",
                      {{1, {{ 1, 5}, { 2, 1}                                                      }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}, {85, 1}, {86, 1}, {87, 1}, {88, 1}}}},
                      testCases);
        add_test_case("ARe",
                      {{1, {{ 2, 1}                           }}},
                      {{1, {{85, 1}, {86, 1}, {87, 1}, {88, 1}}}},
                      testCases);
        add_test_case("ALe",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCases);

        std::vector<TestCase> testCasesP0;
        add_test_case("ADA",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCasesP0);
        add_test_case("ARA",
                      {{1, { }}},
                      {{1, { }}},
                      testCasesP0);
        add_test_case("ALA",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCasesP0);
        add_test_case("ADe",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCasesP0);
        add_test_case("ARe",
                      {{1, { }}},
                      {{1, { }}},
                      testCasesP0);
        add_test_case("ALe",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCasesP0);

        std::vector<TestCase> testCasesP1;
        add_test_case("ADA",
                      {{1, {{ 2, 4}                               }}},
                      {{1, {{27,  4}, {28,  4}, {29,  4}, {30,  4}}}},
                      testCasesP1);
        add_test_case("ARA",
                      {{1, {{ 2, 4}                               }}},
                      {{1, {{27,  4}, {28,  4}, {29,  4}, {30,  4}}}},
                      testCasesP1);
        add_test_case("ALA",
                      {{1, { }}},
                      {{1, { }}},
                      testCasesP1);
        add_test_case("ADe",
                      {{1, {{  2, 1}                              }}},
                      {{1, {{139, 1}, {140, 1}, {141, 1}, {142, 1}}}},
                      testCasesP1);
        add_test_case("ARe",
                      {{1, {{  2, 1}                              }}},
                      {{1, {{139, 1}, {140, 1}, {141, 1}, {142, 1}}}},
                      testCasesP1);
        add_test_case("ALe",
                      {{1, { }}},
                      {{1, { }}},
                      testCasesP1);

        std::vector<TestCase> runTheseCases;
        if(p_size == 1)
            runTheseCases = testCases;
        else if(stk::parallel_machine_rank(pm) == 0)
            runTheseCases = testCasesP0;
        else
            runTheseCases = testCasesP1;

        if (p_size <= 2) {
            for(const auto& test : runTheseCases)
                refine_mesh(test.filename, test. originalSideset, test.refinedSideset);
        }
      }

      void print_parent_child_sideset_connectivity(MeshWithSideset &inputMesh,
                                                   const ElemIdSideVector &parentSideSet,
                                                   const ElemIdSideVector &childSideSet)
      {
          std::ostringstream os;
          bool checkForFamilyTree = true;

          os << "On mesh(" << inputMesh.get_bulk().parallel_rank() << "): " << inputMesh.get_input_filename() << std::endl;
          percept::PerceptMesh &eMesh = inputMesh.get_mesh();

          for(const auto &entry : parentSideSet)
          {
              stk::mesh::Entity parent = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, entry.elem_id);

              if(eMesh.get_bulk_data()->is_valid(parent))
              {
                  os << "  For element: " << eMesh.get_bulk_data()->identifier(parent) << std::endl;

                  for(unsigned i=0; i<childSideSet.size(); ++i)
                  {
                      stk::mesh::Entity child = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, childSideSet[i].elem_id);
                      stk::mesh::Entity parentOfChild = eMesh.getParent(child, checkForFamilyTree);

                      if(parent == parentOfChild)
                          os << "    (child,ordinal) : (" << childSideSet[i].elem_id << "," << childSideSet[i].side_ordinal << ")" << std::endl;
                  }
              }
          }

          std::cerr << os.str() << std::endl;
      }

      ElemIdSideVector get_refined_child_sideset(MeshWithSideset &inputMesh,
                                                 const ElemIdSideVector &parentSideSet)
      {
          percept::PerceptMesh &eMesh = inputMesh.get_mesh();

          ElemIdSideVector childSideSet;
          const bool checkForFamilyTree = true;
          const bool onlyIfElementIsParentLeaf = false;

          for(const auto &entry : parentSideSet)
          {
              stk::mesh::Entity parent = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK, entry.elem_id);

              EXPECT_TRUE(eMesh.get_bulk_data()->is_valid(parent));

              const percept::Elem::RefinementTopology* refTopo = get_refinement_topology(eMesh, parent);

              std::vector< std::pair<UInt,UInt> > refinedSideset = refTopo->get_children_on_ordinal(static_cast<UInt>(entry.side_ordinal));

              std::vector<stk::mesh::Entity> children;
              bool hasChildren = eMesh.getChildren( parent, children, checkForFamilyTree, onlyIfElementIsParentLeaf);

              EXPECT_EQ(hasChildren, (refinedSideset.size() > 0));

              if(hasChildren)
              {
                  EXPECT_GE(children.size(), refinedSideset.size());

                  for(unsigned i=0; i<refinedSideset.size(); ++i)
                  {
                      UInt childIndex = refinedSideset[i].first;
                      UInt faceIndex  = refinedSideset[i].second;

                      childSideSet.push_back({static_cast<int>(eMesh.get_bulk_data()->identifier(children[childIndex])), static_cast<int>(faceIndex)});
                  }
              }
          }

          std::sort(childSideSet.begin(), childSideSet.end(), ElemIdSideLess());

          return childSideSet;
      }

      SideSetIdAndElemIdSidesVector
      get_refined_sideset_from_parent(MeshWithSideset &inputMesh,
                                      const SideSetIdAndElemIdSidesVector& parentSidesetData)
      {
          SideSetIdAndElemIdSidesVector refinedSideset;

          refinedSideset.resize(parentSidesetData.size());

          for(unsigned i=0; i<parentSidesetData.size(); ++i)
          {
              refinedSideset[i].id = parentSidesetData[i].id;
              refinedSideset[i].sideSet = get_refined_child_sideset(inputMesh, parentSidesetData[i].sideSet);
              print_parent_child_sideset_connectivity(inputMesh, parentSidesetData[i].sideSet, refinedSideset[i].sideSet);
          }

          return refinedSideset;
      }

      void test_output_sideset(MeshWithSideset &inputMesh,const SideSetIdAndElemIdSidesVector &refinedSideset)
      {
        inputMesh.write_exodus_file();

        stk::mesh::MetaData meta2;
        stk::unit_test_util::sideset::BulkDataTester bulk2(meta2, MPI_COMM_WORLD);
        stk::unit_test_util::sideset::read_exo_file(bulk2, inputMesh.get_output_filename(), stk::unit_test_util::sideset::READ_ALREADY_DECOMPOSED);
        EXPECT_NO_THROW(compare_sidesets(inputMesh.get_output_filename(), bulk2, refinedSideset));
      }

      void test_refine_sideset(MeshWithSideset &inputMesh,
                               stk::mesh::PartVector &parts,
                               const SideSetIdAndElemIdSidesVector& originalSideset,
                               const SideSetIdAndElemIdSidesVector& expectedRefinedSideset)
      {
        stk::unit_test_util::sideset::BulkDataTester &bulk = inputMesh.get_bulk();
        stk::mesh::MetaData &meta = inputMesh.get_meta();

        bool removeOldElements = false;
        inputMesh.do_uniform_refinement(removeOldElements);

        bulk.initialize_face_adjacent_element_graph();

        stk::mesh::EntityVector faces;
        stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(meta.side_rank()), faces);

        ASSERT_EQ(0u, faces.size());
        EXPECT_TRUE(bulk.get_sideset_ids().empty());

        SideSetIdAndElemIdSidesVector refinedSideset = get_refined_sideset_from_parent(inputMesh, originalSideset);

        for(const SideSetIdAndElemIdSides& sset : refinedSideset)
        {
            stk::mesh::SideSet sideSet = get_stk_side_set(bulk, sset.sideSet);
            bulk.create_sideset(sset.id) =  sideSet;
            bulk.create_side_entities(sideSet, parts);
        }

        EXPECT_NO_THROW(compare_sidesets(inputMesh.get_input_filename(), bulk, expectedRefinedSideset));

        stk::mesh::get_selected_entities(meta.locally_owned_part(), bulk.buckets(meta.side_rank()), faces);

        for(stk::mesh::Entity face : faces)
            ASSERT_EQ(2u, bulk.num_elements(face));

        test_output_sideset(inputMesh, refinedSideset);
      }

      void refine_sideset(const std::string &base_name,
                          const SideSetIdAndElemIdSidesVector& originalSideset,
                          const SideSetIdAndElemIdSidesVector& expectedRefinedSideset)
      {
        MeshWithSideset inputMesh(base_name);

        stk::mesh::MetaData &meta = inputMesh.get_meta();

        stk::mesh::PartVector parts;
        for(const SideSetIdAndElemIdSides &sideSet : originalSideset)
        {
            std::ostringstream os;
            os << "surface_" << sideSet.id;

            stk::mesh::Part &part = meta.declare_part(os.str(), meta.side_rank());
            meta.set_part_id(part, sideSet.id);
            stk::io::put_io_part_attribute(part);
            parts.push_back(&part);
        }

        inputMesh.load_mesh_and_fill_sideset_data();

        EXPECT_TRUE(inputMesh.get_bulk().get_sideset_ids().empty());
        test_refine_sideset(inputMesh, parts, originalSideset, expectedRefinedSideset);
      }

      TEST(regr_stk_sideset_refine, get_refined_sideset)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD;

        const unsigned p_size = stk::parallel_machine_size(pm);

        std::vector<TestCase> testCases;
        add_test_case("AA",
                      {{1, {{ 1, 5}, { 2, 4}                                                      }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}, {27, 4}, {28, 4}, {29, 4}, {30, 4}}}},
                      testCases);
        add_test_case("Ae",
                      {{1, {{ 1, 5}, { 2, 1}                                                      }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}, {85, 1}, {86, 1}, {87, 1}, {88, 1}}}},
                      testCases);

        std::vector<TestCase> testCasesP0;
        add_test_case("AA",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCasesP0);
        add_test_case("Ae",
                      {{1, {{ 1, 5}                           }}},
                      {{1, {{23, 5}, {24, 5}, {25, 5}, {26, 5}}}},
                      testCasesP0);

        std::vector<TestCase> testCasesP1;
        add_test_case("AA",
                      {{1, {{ 2, 4}                               }}},
                      {{1, {{27,  4}, {28,  4}, {29,  4}, {30,  4}}}},
                      testCasesP1);
        add_test_case("Ae",
                      {{1, {{  2, 1}                              }}},
                      {{1, {{139, 1}, {140, 1}, {141, 1}, {142, 1}}}},
                      testCasesP1);

        std::vector<TestCase> runTheseCases;
        if(p_size == 1)
            runTheseCases = testCases;
        else if(stk::parallel_machine_rank(pm) == 0)
            runTheseCases = testCasesP0;
        else
            runTheseCases = testCasesP1;

        if (p_size <= 2) {
            for(const auto& test : runTheseCases)
                refine_sideset(test.filename, test.originalSideset, test.refinedSideset);
        }
      }

      stk::mesh::EntityVector get_entities(const stk::mesh::BulkData &bulk, stk::topology::rank_t rank, const stk::mesh::EntityIdVector &ids)
      {
          stk::mesh::EntityVector entities;
          for(stk::mesh::EntityId id : ids)
          {
              stk::mesh::Entity entity = bulk.get_entity(rank, id);
              EXPECT_TRUE(bulk.is_valid(entity));
              entities.push_back(entity);
          }
          return entities;
      }

      std::vector< PseudoEntity > get_pseudo_side_entities(const stk::mesh::BulkData &bulk, const stk::mesh::EntityVector &sides)
      {
          std::vector< PseudoEntity > pseudoSides;
          for(const stk::mesh::Entity side : sides)
          {
              ThrowRequire(bulk.entity_rank(side) == bulk.mesh_meta_data().side_rank());

              const stk::mesh::Entity *nodes = bulk.begin_nodes(side);
              pseudoSides.push_back( stk::mesh::EntityVector(nodes, nodes + bulk.num_nodes(side)) );
          }
          return pseudoSides;
      }

      void fill_sideset_from_elements_and_pseudo_side(const stk::mesh::BulkData &bulk,
                                                      const stk::mesh::EntityVector &elements,
                                                      const PseudoEntity &nodes,
                                                      stk::mesh::SideSet &sideset)
      {
          for(const stk::mesh::Entity element : elements)
          {
              stk::mesh::OrdinalAndPermutation ordinalAndPermutation = stk::mesh::get_ordinal_and_permutation(bulk,
                                                                                                              element,
                                                                                                              bulk.mesh_meta_data().side_rank(),
                                                                                                              nodes);

              if(stk::mesh::INVALID_CONNECTIVITY_ORDINAL != ordinalAndPermutation.first)
              {
                  sideset.push_back({element, ordinalAndPermutation.first});
                  break;
              }
          }
      }

      stk::mesh::SideSet get_sideset_from_elements_and_pseudo_sides(const stk::mesh::BulkData &bulk,
                                                                    const stk::mesh::EntityVector &elements,
                                                                    const std::vector< PseudoEntity > &pseudoSides)

      {
          stk::mesh::SideSet sideset;
          for(const PseudoEntity &pseudoSide : pseudoSides)
              fill_sideset_from_elements_and_pseudo_side(bulk, elements, pseudoSide, sideset);

          return sideset;
      }

      void test_parent_info(const stk::mesh::BulkData &bulk,
                            const stk::mesh::Entity &parentElement,
                            const stk::mesh::Entity &parentSide,
                            const stk::mesh::ConnectivityOrdinal parentSideOrdinal)
      {
          unsigned ordinalIndex = bulk.find_ordinal(parentElement, bulk.mesh_meta_data().side_rank(), parentSideOrdinal);
          unsigned numConnectivity = bulk.num_connectivity(parentElement, bulk.mesh_meta_data().side_rank());
          EXPECT_TRUE(ordinalIndex < numConnectivity);
          EXPECT_TRUE(parentSide == bulk.begin(parentElement, bulk.mesh_meta_data().side_rank())[ordinalIndex]);
      }

      stk::mesh::SideSet get_refined_sideset_from_parent_child_info(const stk::mesh::BulkData &bulk,
                                                                    const stk::mesh::Entity &parentElement,
                                                                    const stk::mesh::EntityVector &childElements,
                                                                    const stk::mesh::Entity & parentSide,
                                                                    const stk::mesh::EntityVector &childSides,
                                                                    const stk::mesh::ConnectivityOrdinal parentSideOrdinal)
      {
          test_parent_info(bulk, parentElement, parentSide, parentSideOrdinal);
          return get_sideset_from_elements_and_pseudo_sides(bulk, childElements, get_pseudo_side_entities(bulk, childSides));
      }

      void test_refined_sideset_from_parent_child_info(const std::string &inputFileName,
                                                       const stk::mesh::ConnectivityOrdinal &parentSideOrdinal,
                                                       const stk::mesh::EntityId &parentElementId,
                                                       const stk::mesh::EntityIdVector &childElementIds,
                                                       const stk::mesh::EntityId &parentSideId,
                                                       const stk::mesh::EntityIdVector &childSideIds,
                                                       const ElemIdSideVector &expectedElemIdSides)
      {
          MeshWithSideset inputMesh(inputFileName);
          inputMesh.load_mesh_and_fill_sideset_data();

          stk::unit_test_util::sideset::BulkDataTester &bulk = inputMesh.get_bulk();
          stk::mesh::MetaData &meta = inputMesh.get_meta();

          inputMesh.do_uniform_refinement(false);

          stk::mesh::Entity parentElement = bulk.get_entity(stk::topology::ELEMENT_RANK, parentElementId);
          stk::mesh::EntityVector childElements = get_entities(bulk, stk::topology::ELEMENT_RANK, childElementIds);

          stk::mesh::Entity parentSide = bulk.get_entity(meta.side_rank(), parentSideId);
          stk::mesh::EntityVector childSides = get_entities(bulk, meta.side_rank(), childSideIds);

          stk::mesh::SideSet sideSet = get_refined_sideset_from_parent_child_info(bulk, parentElement, childElements, parentSide, childSides, parentSideOrdinal);

          ASSERT_EQ(expectedElemIdSides.size(), sideSet.size()) << "for file: " << inputFileName;

          for(size_t i=0;i<sideSet.size();++i)
          {
              EXPECT_EQ(static_cast<stk::mesh::EntityId>(expectedElemIdSides[i].elem_id), bulk.identifier(sideSet[i].element)) << "for file: " << inputFileName;
              EXPECT_EQ(expectedElemIdSides[i].side_ordinal, sideSet[i].side) << "for file: " << inputFileName;
          }
      }

      TEST(regr_stk_sideset_refine, refined_child_ordinal)
      {
        stk::ParallelMachine pm = MPI_COMM_WORLD;

        const unsigned p_size = stk::parallel_machine_size(pm);

        if(1 != p_size) return;

        test_refined_sideset_from_parent_child_info("ARA",
                                                    static_cast<stk::mesh::ConnectivityOrdinal>(4),
                                                    2u,
                                                    {27, 28, 29, 30, 31, 32, 33, 34},
                                                    16u,
                                                    {17, 18, 19, 20},
                                                    {{27, 4}, {28, 4}, {29, 4}, {30, 4}});

        test_refined_sideset_from_parent_child_info("ALA",
                                                    static_cast<stk::mesh::ConnectivityOrdinal>(5),
                                                    1u,
                                                    {19, 20, 21, 22, 23, 24, 25, 26},
                                                    16u,
                                                    {17, 18, 19, 20},
                                                    {{23, 5}, {24, 5}, {25, 5}, {26, 5}});

        test_refined_sideset_from_parent_child_info("ADA",
                                                    static_cast<stk::mesh::ConnectivityOrdinal>(4),
                                                    2u,
                                                    {27, 28, 29, 30, 31, 32, 33, 34},
                                                    16u,
                                                    {17, 18, 19, 20},
                                                    {{27, 4}, {28, 4}, {29, 4}, {30, 4}});

        test_refined_sideset_from_parent_child_info("ADA",
                                                    static_cast<stk::mesh::ConnectivityOrdinal>(5),
                                                    1u,
                                                    {19, 20, 21, 22, 23, 24, 25, 26},
                                                    16u,
                                                    {17, 18, 19, 20},
                                                    {{23, 5}, {24, 5}, {25, 5}, {26, 5}});

//        ARe and ALe still fail due to percept flipping the connected ordinal on the elements (including parent)
//        test_refined_sideset_from_parent_child_info("ARe",
//                                                    static_cast<stk::mesh::ConnectivityOrdinal>(1),
//                                                    2u,
//                                                    {85, 86, 87, 88},
//                                                    16u,
//                                                    {17, 18, 19, 20},
//                                                    {{85, 1}, {86, 1}, {87, 1}, {88, 1}});
//        test_refined_sideset_from_parent_child_info("ALe",
//                                                    static_cast<stk::mesh::ConnectivityOrdinal>(5),
//                                                    1u,
//                                                    {19, 20, 21, 22, 23, 24, 25, 26},
//                                                    16u,
//                                                    {17, 18, 19, 20},
//                                                    {{23, 5}, {24, 5}, {25, 5}, {26, 5}});
      }

    }//    namespace regression_test
  }//  namespace percept

