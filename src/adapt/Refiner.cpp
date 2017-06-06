// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <exception>
#include <fstream>
#include <set>
#include <typeinfo>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <adapt/Refiner.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/FixSideSets.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/MeshUtil.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

#include <adapt/NodeRegistryDef.hpp>
#include <percept/RebalanceMesh.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <percept/mesh/geometry/kernel/GeometryKernelGregoryPatch.hpp>
#include <percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#if HAVE_CUBIT
#include <percept/mesh/geometry/kernel/GeometryKernelPGEOM.hpp>
#endif
#include <percept/mesh/geometry/kernel/MeshGeometry.hpp>
#include <percept/mesh/geometry/kernel/GeometryFactory.hpp>

#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#if ENABLE_SMOOTHER3
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherNewton.hpp>
#endif
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherAlgebraic.hpp>
#endif

#include <adapt/FixSideSetsSelector.hpp>
#include <adapt/AdaptedMeshVerifier.hpp>

#ifndef NDEBUG
#include <percept/mesh/geometry/volume/VolumeUtil.hpp>
#endif

#define DO_TIMING_REFINER DO_TIMING

#if DO_TIMING_REFINER
#define TIMING(code) code
#define TIMER(name) stk::diag::Timer timer ## name ( #name, rootTimer());  stk::diag::TimeBlock tbTimer ## name (timer ## name)
#define TIMER2(name,parentName) stk::diag::Timer timer ## name ( #name, timer ## parentName);  stk::diag::TimeBlock tbTimer ## name (timer ## name)
#else
#define TIMING(code) do {  } while(0)
#define TIMER(name) do {} while(0)
#define TIMER2(name,parentName) do {} while(0)
#endif

#define USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND 0
#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
#include "/usr/netpub/valgrind-3.8.1/include/valgrind/callgrind.h"
#endif

#include <stk_mesh/base/MeshUtils.hpp>

  namespace percept {
    extern bool s_do_transition_break;

    using namespace std;
    using namespace percept;

    NodeRegistry *s_nodeRegistry = 0;

    //#define CHK(nr) do { std::ostringstream str; str << __FILE__ << ":" << __LINE__; check(nr, str.str()); } while (0)
#define CHK(eMesh) do { } while (0)

#if 0
      static void check(NodeRegistry* nr, const std::string& where)
      {

        PerceptMesh& eMesh = nr->getMesh();
        int id = 2;
        stk::mesh::Entity entity = eMesh.get_bulk_data()->get_entity(eMesh.element_rank(), id);
        std::ostringstream str;
        //str << where;
        if(eMesh.is_valid(entity))
          {
            str << "P[" << eMesh.get_rank() << "] CHK: " << where << " elem= " << eMesh.print_entity_compact(entity);
            str << "\n";
            nr->query(str, id, (int)eMesh.side_rank(), 4);
          }
        else
          {
            str << " not valid";
          }
        std::cout << str.str() << std::endl;
      }
#endif


    stk::diag::Timer& Refiner::rootTimer()
    {
      if (m_alternateRootTimer)
        {
          return *m_alternateRootTimer;
        }
      else
        {
          return m_timer;
        }
    }

    Refiner::Refiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase &  bp, stk::mesh::FieldBase *proc_rank_field) :
      m_eMesh(eMesh), m_breakPattern(),
      m_nodeRegistry(0),
      m_proc_rank_field(proc_rank_field), m_doRemove(true), m_ranks(), m_ignoreSideSets(false),
      m_geomFile(""), m_geomSnap(false),
      m_refinementInfo(this),
      m_doQueryOnly(false),
      m_progress_meter_frequency(20),
      m_doProgress(false),
      m_alwaysInitNodeRegistry(true),
      m_doSmoothGeometry(false),
      m_removeGeometryBlocks(false),
      m_fixAllBlockBoundaries(false),
      m_needsRemesh(true),
      m_doLevelBasedUnrefinement(false)
      ,m_alternateRootTimer(0)
      ,m_modBegEndRootTimer(0)
      ,m_refinerSelector(0)
      ,m_doAddChildrenToParts(true)
      ,m_avoidFixSideSets(false)
      ,m_avoidFixSideSetChecks(false)
      ,m_avoidClearDanglingNodes(false)
      ,m_onlyOneLevelUnrefine(false)
      ,m_doRebalance(false)
      ,m_rebalThreshold(1.0)
      ,m_removeFromNewNodesPart(true)
      ,m_do_new_elements(true)
      ,m_timerSet(sierra::Diag::TIMER_ALL)
      ,m_timer(stk::diag::createRootTimer("Refiner", m_timerSet))

    {
      bp.setSubPatterns(m_breakPattern, eMesh);
      m_nodeRegistry = new NodeRegistry (m_eMesh, this);
      s_nodeRegistry = m_nodeRegistry;
      m_nodeRegistry->initialize();
      m_nodeRegistry->init_comm_all();
      m_allocated = false;
      m_timer.start();
    }

    //NLM new constructor for Python users
    //takes in a Pattern refine_type to determine which UniformRefinerPatternBase will be used

    Refiner::Refiner(percept::PerceptMesh& eMesh, Pattern refine_type, stk::mesh::FieldBase *proc_rank_field) :
      m_eMesh(eMesh), m_breakPattern(),
      m_nodeRegistry(0),
      m_proc_rank_field(proc_rank_field), m_doRemove(true), m_ranks(), m_ignoreSideSets(false),
      m_geomFile(""), m_geomSnap(false),
      m_refinementInfo(this),
      m_doQueryOnly(false),
      m_progress_meter_frequency(20),
      m_doProgress(false),
      m_alwaysInitNodeRegistry(true),
      m_doSmoothGeometry(false),
      m_removeGeometryBlocks(false),
      m_fixAllBlockBoundaries(false),
      m_needsRemesh(true),
      m_doLevelBasedUnrefinement(false)
      ,m_alternateRootTimer(0)
      ,m_modBegEndRootTimer(0)
      ,m_refinerSelector(0)
      ,m_doAddChildrenToParts(true)
      ,m_avoidFixSideSets(false)
      ,m_avoidFixSideSetChecks(false)
      ,m_avoidClearDanglingNodes(false)
      ,m_onlyOneLevelUnrefine(false)
      ,m_doRebalance(false)
      ,m_rebalThreshold(1.0)
      ,m_removeFromNewNodesPart(true)
      ,m_do_new_elements(true)
      ,m_timerSet(sierra::Diag::TIMER_ALL)
      ,m_timer(stk::diag::createRootTimer("Refiner", m_timerSet))
    {
      m_nodeRegistry = new NodeRegistry (m_eMesh, this);
      s_nodeRegistry = m_nodeRegistry;
      m_nodeRegistry->initialize();
      m_nodeRegistry->init_comm_all();
      m_allocated = true;
      m_timer.start();
      UniformRefinerPatternBase* bp;
      switch(refine_type)
        {
        case LINE2_LINE2_2:
          bp = (UniformRefinerPatternBase*) new Line2_Line2_2(eMesh);
          break;
        case BEAM2_BEAM2_2:
          bp = (UniformRefinerPatternBase*) new Beam2_Beam2_2(eMesh);
          break;
        case SHELLLINE2_SHELLLINE2_2:
          bp = (UniformRefinerPatternBase*) new ShellLine2_ShellLine2_2(eMesh);
          break;
        case SHELLLINE3_SHELLLINE3_2:
          bp = (UniformRefinerPatternBase*) new ShellLine3_ShellLine3_2(eMesh);
          break;
        case QUAD4_QUAD4_4_OLD:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad4_4_Old(eMesh);
          break;
        case QUAD4_QUAD4_4:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad4_4(eMesh);
          break;
        case QUAD4_QUAD4_4_SIERRA:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad4_4_Sierra(eMesh);
          break;
        case TRI3_TRI3_4:
          bp = (UniformRefinerPatternBase*) new Tri3_Tri3_4(eMesh);
          break;
        case SHELLTRI3_SHELLTRI3_4:
          bp = (UniformRefinerPatternBase*) new ShellTri3_ShellTri3_4(eMesh);
          break;
        case SHELLTRI6_SHELLTRI6_4:
          bp = (UniformRefinerPatternBase*) new ShellTri6_ShellTri6_4(eMesh);
          break;
        case SHELLQUAD4_SHELLQUAD4_4:
          bp = (UniformRefinerPatternBase*) new ShellQuad4_ShellQuad4_4(eMesh);
          break;
        case SHELLQUAD8_SHELLQUAD8_4:
          bp = (UniformRefinerPatternBase*) new ShellQuad8_ShellQuad8_4(eMesh);
          break;
        case TET4_TET4_8:
          bp = (UniformRefinerPatternBase*) new Tet4_Tet4_8(eMesh);
          break;
        case HEX8_HEX8_8:
          bp = (UniformRefinerPatternBase*) new Hex8_Hex8_8(eMesh);
          break;
        case WEDGE6_WEDGE6_8:
          bp = (UniformRefinerPatternBase*) new Wedge6_Wedge6_8(eMesh);
          break;
        case PYRAMID5_PYRAMID5_10:
          bp = (UniformRefinerPatternBase*) new Pyramid5_Pyramid5_10(eMesh);
          break;
        case LINE3_LINE3_2:
          bp = (UniformRefinerPatternBase*) new Line3_Line3_2(eMesh);
          break;
        case BEAM3_BEAM3_2:
          bp = (UniformRefinerPatternBase*) new Beam3_Beam3_2(eMesh);
          break;
        case TRI6_TRI6_4:
          bp = (UniformRefinerPatternBase*) new Tri6_Tri6_4(eMesh);
          break;
        case QUAD9_QUAD9_4:
          bp = (UniformRefinerPatternBase*) new Quad9_Quad9_4(eMesh);
          break;
        case QUAD8_QUAD8_4:
          bp = (UniformRefinerPatternBase*) new Quad8_Quad8_4(eMesh);
          break;
        case HEX27_HEX27_8:
          bp = (UniformRefinerPatternBase*) new Hex27_Hex27_8(eMesh);
          break;
        case HEX20_HEX20_8:
          bp = (UniformRefinerPatternBase*) new Hex20_Hex20_8(eMesh);
          break;
        case TET10_TET10_8:
          bp = (UniformRefinerPatternBase*) new Tet10_Tet10_8(eMesh);
          break;
        case WEDGE15_WEDGE15_8:
          bp = (UniformRefinerPatternBase*) new Wedge15_Wedge15_8(eMesh);
          break;
        case WEDGE18_WEDGE18_8:
          bp = (UniformRefinerPatternBase*) new Wedge18_Wedge18_8(eMesh);
          break;
        case PYRAMID13_PYRAMID13_10:
          bp = (UniformRefinerPatternBase*) new Pyramid13_Pyramid13_10(eMesh);
          break;
        case QUAD4_QUAD9_1:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad9_1(eMesh);
          break;
        case QUAD4_QUAD8_1:
          bp = (UniformRefinerPatternBase*) new Quad4_Quad8_1(eMesh);
          break;
        case BEAM2_BEAM3_1:
          bp = (UniformRefinerPatternBase*) new Beam2_Beam3_1(eMesh);
          break;
        case SHELLQUAD4_SHELLQUAD8_1:
          bp = (UniformRefinerPatternBase*) new ShellQuad4_ShellQuad8_1(eMesh);
          break;
        case TRI3_TRI6_1:
          bp = (UniformRefinerPatternBase*) new Tri3_Tri6_1(eMesh);
          break;
        case TET4_TET10_1:
          bp = (UniformRefinerPatternBase*) new Tet4_Tet10_1(eMesh);
          break;
        case HEX8_HEX27_1:
          bp = (UniformRefinerPatternBase*) new Hex8_Hex27_1(eMesh);
          break;
        case HEX8_HEX20_1:
          bp = (UniformRefinerPatternBase*) new Hex8_Hex20_1(eMesh);
          break;
        case WEDGE6_WEDGE15_1:
          bp = (UniformRefinerPatternBase*) new Wedge6_Wedge15_1(eMesh);
          break;
        case WEDGE6_WEDGE18_1:
          bp = (UniformRefinerPatternBase*) new Wedge6_Wedge18_1(eMesh);
          break;
        case PYRAMID5_PYRAMID13_1:
          bp = (UniformRefinerPatternBase*) new Pyramid5_Pyramid13_1(eMesh);
          break;
        case QUAD4_TRI3_2:
          bp = (UniformRefinerPatternBase*) new Quad4_Tri3_2(eMesh);
          break;
        case QUAD4_TRI3_4:
          bp = (UniformRefinerPatternBase*) new Quad4_Tri3_4(eMesh);
          break;
        case QUAD4_TRI3_6:
          bp = (UniformRefinerPatternBase*) new Quad4_Tri3_6(eMesh);
          break;
        case HEX8_TET4_24:
          bp = (UniformRefinerPatternBase*) new Hex8_Tet4_24(eMesh);
          break;
        case HEX8_TET4_6_12:
          bp = (UniformRefinerPatternBase*) new Hex8_Tet4_6_12(eMesh);
          break;
        case TET4_HEX8_4:
          bp = (UniformRefinerPatternBase*) new Tet4_Hex8_4(eMesh);
          break;
        case WEDGE6_HEX8_6:
          bp = (UniformRefinerPatternBase*) new Wedge6_Hex8_6(eMesh);
          break;
        case TRI3_QUAD4_3:
          bp = (UniformRefinerPatternBase*) new Tri3_Quad4_3(eMesh);
          break;
        default:
          throw std::runtime_error("Refiner::Refiner Unrecognized refinement pattern");
          break;
        }

      bp->setSubPatterns(m_breakPattern, eMesh);

      m_proc_rank_field = proc_rank_field;
    }

    Refiner::~Refiner()
    {
      if (m_nodeRegistry)
        delete m_nodeRegistry;
      if (m_allocated)
      {
        //for (unsigned int i = 0; i < m_breakPattern.size(); i++)
        // SRK only delete the first one - others are handled by the parent 0'th break pattern
        for (unsigned int i = 0; i < 1; i++)
        {
          delete m_breakPattern[i];
        }
      }
    }

    void Refiner::doProgressPrint(const std::string& msg)
    {
      if (m_doProgress && m_eMesh.get_rank() == 0)
        {
          size_t now=0, hwm=0;
          stk::get_memory_usage(now, hwm);
          std::cout << msg << " cpu: " << m_eMesh.cpu_time() << " [sec] mem= " << stk::human_bytes(now) << " [now_proc0] " << stk::human_bytes(hwm) << " [hwm_proc0]" << std::endl;
        }
    }

#define EXTRA_PRINT_UR_GETBLOCKS 0

    void Refiner::
    setGeometryFile(std::string file_name) { m_geomFile = file_name;
      m_geomSnap = true; }


    void Refiner::
    setRemoveOldElements(bool do_remove) { m_doRemove = do_remove; }

    bool Refiner::
    getRemoveOldElements() { return m_doRemove; }

    void Refiner::
    setDoProgressMeter(bool do_progress) { m_doProgress = do_progress; }

    bool Refiner::
    getDoProgressMeter() { return m_doProgress; }

    void Refiner::
    setIgnoreSideSets(bool ignore_ss)
    {
      m_ignoreSideSets= ignore_ss;
    }

    bool Refiner::
    getIgnoreSideSets()
    {
      return m_ignoreSideSets;
    }

    void Refiner::
    removeFromOldPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;

      std::string oldPartName = breakPattern->getOldElementsPartName()+toString(rank);
      stk::mesh::Part *oldPart = m_eMesh.get_fem_meta_data()->get_part(oldPartName);
      //std::cout << "tmp removeFromOldPart:: oldPartName= " << oldPartName << std::endl;
      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null");
        }

      stk::mesh::PartVector remove_parts(1, oldPart);
      stk::mesh::PartVector add_parts;
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

      std::vector<stk::mesh::Entity> elems;
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k)) // && fromPartsSelector(**k) )
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                {
                  stk::mesh::Entity element = bucket[i_element];
                  elems.push_back(element);
                }
            }
        }

      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          //std::cout << "tmp removeFromOldPart:: element = " << *elems[ielem] << std::endl;
          m_eMesh.get_bulk_data()->change_entity_parts( elems[ielem], add_parts, remove_parts );
        }
    }

    void Refiner::
    addOldElementsToPart(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern, unsigned *elementType)
    {
      EXCEPTWATCH;
      std::string oldPartName = breakPattern->getOldElementsPartName()+toString(rank);
      stk::mesh::Part *oldPart = m_eMesh.get_fem_meta_data()->get_part(oldPartName);

#define DEBUG_REMOVE_OLD_PARTS 0

      if (DEBUG_REMOVE_OLD_PARTS) std::cout << "tmp addOldElementsToPart:: oldPartName= " << oldPartName << std::endl;

      if (!oldPart)
        {
          std::cout << "oldPartName= " << oldPartName << std::endl;
          throw std::runtime_error("oldpart is null");
        }

      stk::mesh::PartVector add_parts(1, oldPart);
      stk::mesh::PartVector remove_parts;
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

      // The list of Parts that this break pattern will refine.  Only remove elements belonging to these parts.
      stk::mesh::Selector fromPartsSelector = stk::mesh::selectUnion( breakPattern->getFromParts() );

      std::vector<stk::mesh::Entity> elems;
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      unsigned nele=0;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k) && fromPartsSelector(**k) )
            {
              stk::mesh::Bucket & bucket = **k ;
              const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(bucket);
              shards::CellTopology topo(bucket_cell_topo_data);

              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                {
                  EXCEPTWATCH;
                  stk::mesh::Entity element = bucket[i_element];
                  if (!m_eMesh.is_valid(element))
                    {
                      std::cout << "element = 0" << std::endl;
                      throw std::runtime_error("element = 0");
                    }

                  if (elementType && (topo.getKey() != *elementType))
                    {
                    }
                  else
                    {
                      elems.push_back(element);
                      ++nele;
                      if (DEBUG_REMOVE_OLD_PARTS) {
                        std::cout << "tmp adding to oldParts = "; m_eMesh.print(element); std::cout << std::endl;
                      }
                    }
                }
            }
        }


      for (unsigned ielem=0; ielem < elems.size(); ielem++)
        {
          if (DEBUG_REMOVE_OLD_PARTS) {
            std::cout << "tmp addOldElementsToPart element = "; m_eMesh.print(elems[ielem]); std::cout << std::endl;
          }
          m_eMesh.get_bulk_data()->change_entity_parts( elems[ielem], add_parts, remove_parts );
        }

    }

    void Refiner::
    trace_print(std::string msg)
    {
      size_t heap_in_Mb = 0;
      size_t memory_in_Mb = 0;
      double cpu = 0.0;

      if (TRACE_STAGE_PRINT)
        {
          memory_in_Mb = Util::memory(heap_in_Mb);
          memory_in_Mb = memory_in_Mb / (1024*1024);
          heap_in_Mb = heap_in_Mb / (1024*1024);
        }

      cpu = Util::cpu_time();
      std::cout
        << msg
        << " cpu_time= " << cpu/(60.) << " [min] "
        << " mem= " << memory_in_Mb << " [Mb] "
        //<< " heap= " << heap_in_Mb << " [Mb] "
        << std::endl;

    }


    struct myVec
    {
      double *data;
      int len;
      int res;
    };

    static void doPrintSizes()
    {
      if (0)
        {
          std::cout
            << "sizeof(myVec) = " << sizeof(myVec) << " "
            << "sizeof(Relation) = " << sizeof(stk::mesh::Relation) << " "
            << "sizeof(Entity) = " << sizeof(stk::mesh::Entity) << " "
            << "\nsizeof(EntityKey) = " << sizeof(stk::mesh::EntityKey) << " "
            << "\nsizeof(RelationVector) = " << sizeof(stk::mesh::RelationVector) << " "
            << "\nsizeof(EntityCommInfoVector) = " << sizeof(stk::mesh::EntityCommInfoVector) << " "
            << "\nsizeof(Bucket *) = " << sizeof(stk::mesh::Bucket *) << " "
            << "\nsizeof(unsigned) = " << sizeof(unsigned) << " "
            << "\nsizeof(size_t) = " << sizeof(size_t) << " "
            << "\nsizeof(EntityModificationLog) = " << sizeof(stk::mesh::EntityState) << std::endl;

        }
    }

    void Refiner::
    checkBreakPatternValidityAndBuildRanks(std::vector<stk::mesh::EntityRank>& ranks)
    {
      ranks.resize(0);

      if (0) doPrintPatterns();

      if (m_doRemove)
        mod_begin();

      for (unsigned ibp = 0; ibp < m_breakPattern.size(); ibp++)
        {
          if (m_breakPattern[ibp])
            {
              stk::mesh::EntityRank irank = m_breakPattern[ibp]->getPrimaryEntityRank();
              stk::mesh::EntityRank irank_prev = static_cast<stk::mesh::EntityRank>(percept::EntityRankEnd);

              if (ibp > 0) irank_prev = m_breakPattern[ibp-1]->getPrimaryEntityRank();
              if (irank > irank_prev)
                {
                  for (unsigned jbp = 0; jbp < m_breakPattern.size(); jbp++)
                    {
                      if (m_breakPattern[jbp])
                        {
                          stk::mesh::EntityRank jrank = m_breakPattern[jbp]->getPrimaryEntityRank();
                          std::cout << "tmp jbp= " << jbp << " jrank= " << jrank
                                    << " fromTopo= " << m_breakPattern[jbp]->getFromTopology()->name
                                    << " toTopo= " << m_breakPattern[jbp]->getToTopology()->name
                                    << std::endl;
                        }
                    }

                  std::cout << "tmp irank= " << irank << " irank_prev= " << irank_prev << std::endl;
                  throw std::logic_error("m_breakPattern: must be in decreasing order of rank");
                }
              ranks.push_back(irank);
              if (m_doRemove)
                {
                  unsigned elementType = m_breakPattern[ibp]->getFromTypeKey();
                  addOldElementsToPart(irank, m_breakPattern[ibp], &elementType);
                }
            }
          else
            {
              std::cout << "ibp = " << ibp << std::endl;
              throw std::logic_error("m_breakPattern is null");
            }
        }
      if (m_doRemove)
        mod_end(0,"RefinerCBR");

    }

    void Refiner::
    doPrintPatterns()
    {
          for (unsigned jbp = 0; jbp < m_breakPattern.size(); jbp++)
            {
              if (m_breakPattern[jbp])
                {
                  stk::mesh::EntityRank jrank = m_breakPattern[jbp]->getPrimaryEntityRank();
                  std::cout << "tmp jbp= " << jbp << " jrank= " << jrank
#ifndef __clang__
                            << " class= " << m_eMesh.demangle(typeid(*m_breakPattern[jbp]).name())
#endif
                            << " fromTopo= " << m_breakPattern[jbp]->getFromTopology()->name
                            << " toTopo= " << m_breakPattern[jbp]->getToTopology()->name
                            << std::endl;
                }
            }
    }

    void Refiner::
    refineMethodApply(NodeRegistry::ElementFunctionPrototype function, const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks, const CellTopologyData * const bucket_topo_data)
    {
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(element);

      shards::CellTopology cell_topo(cell_topo_data);
      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;
          stk::mesh::EntityRank needed_entity_rank = needed_entity_ranks[ineed_ent].first;

          if (needed_entity_rank == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_rank == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_rank == stk::topology::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              /// note: at this level of granularity we can do single edge refinement, hanging nodes, etc.
              //SubDimCell_SDCEntityType subDimEntity;
              //getSubDimEntity(subDimEntity, element, needed_entity_rank, iSubDimOrd);
              bool doMark = true;
              if (needed_entity_ranks[ineed_ent].third.size())
                {
                  VERIFY_OP_ON(needed_entity_ranks[ineed_ent].third.size(), ==, numSubDimNeededEntities, "bad size");
                  if (!needed_entity_ranks[ineed_ent].third[iSubDimOrd])
                    doMark = false;
                }
              if (doMark)
                (m_nodeRegistry ->* function)(element, needed_entity_ranks[ineed_ent], iSubDimOrd, true, bucket_topo_data);

            } // iSubDimOrd
        } // ineed_ent
    }

    bool isShell(PerceptMesh& eMesh, stk::mesh::Entity element)
    {
      // shards::CellTopology element_topo(eMesh.get_cell_topology(element));
      // int topoDim = UniformRefinerPatternBase::getTopoDim(element_topo);

      // bool isShell = false;
      // if (topoDim < (int)eMesh.entity_rank(element))
      //   {
      //     isShell = true;
      //   }
      // return isShell;
      stk::topology::topology_t topo = eMesh.topology(element);
      switch(topo)
        {
        case stk::topology::SHELL_LINE_2:
        case stk::topology::SHELL_LINE_3:
        case stk::topology::SHELL_TRI_3:
        case stk::topology::SHELL_TRI_4:
        case stk::topology::SHELL_TRI_6:
        case stk::topology::SHELL_QUAD_4:
        case stk::topology::SHELL_QUAD_8:
        case stk::topology::SHELL_QUAD_9:
          return true;
        default:
          return false;
        }
      return false;
    }



    // called by doMark
    void Refiner::
    preMark(int iter, int num_registration_loops) {}

    //void Refiner::
    //mark(int iter, int num_registration_loops) {}

    bool Refiner::
    postMark(int iter, int num_registration_loops) { return false; }

    void Refiner::
    doBreak(int num_registration_loops)
    {
      EXCEPTWATCH;

      stk::diag::setEnabledTimerMetricsMask(stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME);
      stk::diag::Timer timerAdapt_("percept::DoBreak", rootTimer());
      stk::diag::TimeBlock tbTimerAdapt_(timerAdapt_);
      stk::diag::TimeBlock tbRoot_(rootTimer());

      m_eMesh.get_bulk_data()->delete_face_adjacent_element_graph();

      initializeRefine();

      if (m_doQueryOnly) return;

      doMark(num_registration_loops);

      doRefine();
    }

    void Refiner::
    initializeRefine()
    {
      stk::diag::Timer timerAdapt_("percept::Refine", rootTimer());
      stk::diag::TimeBlock timerAdaptBlock_(timerAdapt_);
      stk::diag::Timer timerInitRefine_("percept::InitRefine", timerAdapt_);
      stk::diag::TimeBlock timerInitRefineBlock_(timerInitRefine_);

      REF_LTRACE("initializeRefine: start");

      {
        TIMER2(IR_dumpDB,InitRefine_);
        m_nodeRegistry->dumpDB("start of doBreak");
      }

      {
        TIMER2(IR_req_same_proc,InitRefine_);
        if (m_eMesh.getProperty("percept_Refiner_require_sides") != "false")
          require_sides_on_same_proc_as_pos_perm_element();
      }

      if (0) doPrintSizes();

      {
        TIMER2(IR_side_part,InitRefine_);

        if (m_eMesh.getProperty("use_side_map") == "true")
          get_side_part_relations(m_eMesh, false, m_side_part_map);
      }

      {
        TIMER2(IR_checkPolarity,InitRefine_);
        m_eMesh.setProperty("AdaptedMeshVerifier::checkPolarity","initializeRefine");
        AdaptedMeshVerifier::checkPolarity(m_eMesh);
      }

      REF_LTRACE("initializeRefine: after checkPolarity");

      {
        TIMER2(IR_checkBPVali,InitRefine_);
        std::vector<stk::mesh::EntityRank>& ranks = m_ranks;
        // check logic of break pattern setup and also build ranks used vector
        checkBreakPatternValidityAndBuildRanks(ranks);
      }

      {
        TIMER2(IR_getRInf, InitRefine_);
        getRefinementInfo(m_ranks);
      }
    }

    void Refiner::
    getRefinementInfo(std::vector<stk::mesh::EntityRank>& ranks)
    {
      ///////////////////////////////////////////////////////////
      ///// Get info on refinements that will be done
      ///////////////////////////////////////////////////////////
      if (1)
        {

          unsigned num_orig_nodes=0;
          {
            std::vector<size_t> count1 ;

            stk::mesh::count_entities(stk::mesh::Selector(m_eMesh.get_fem_meta_data()->universal_part()) , *m_eMesh.get_bulk_data(), count1 );
            if (count1.size() < 3)
              {
                throw std::logic_error("logic error in Refiner m_refinementInfoaByType");
              }
            num_orig_nodes = count1[0];
          }

          m_refinementInfo.m_refinementInfoByType.resize(ranks.size());

          stk::mesh::PartVector fromPartsAll;

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              stk::mesh::PartVector * fromParts = &(m_breakPattern[irank]->getFromParts());
              if (fromParts)
                {
                  for (unsigned ipart = 0; ipart < fromParts->size(); ipart++)
                    {
                      // print message about beams && shells
                      if (1 && m_eMesh.get_rank()==0)
                        {
                          stk::mesh::Part *part = (*fromParts)[ipart];
                          const CellTopologyData * part_cell_topo_data = m_eMesh.get_cell_topology(*part);
                          if (part_cell_topo_data)
                            {
                              shards::CellTopology part_cell_topo(part_cell_topo_data);
                              std::string ptopo_name = part_cell_topo.getName();
                              bool printAll = false;
                              if (printAll)
                                {
                                  std::cout << "P[0] INFO: will refine block: " << part->name() << " with topology= " << part_cell_topo.getName() << std::endl;
                                }
                            }
                        }

                      fromPartsAll.push_back((*fromParts)[ipart]);
                    }
                }
            }

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              EXCEPTWATCH;
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              shards::CellTopology cell_topo(m_breakPattern[irank]->getFromTopology());
              //std::cout << "a1 toposrk cell_topo= " << cell_topo.getName() << std::endl;

              stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->locally_owned_part());
              if (fromPartsAll.size())
                {
                  selector = stk::mesh::Selector();
                  for (unsigned ipart = 0; ipart < fromPartsAll.size(); ipart++)
                    {
                      stk::mesh::Part *part = fromPartsAll[ipart];
                      const CellTopologyData * part_cell_topo_data = m_eMesh.get_cell_topology(*part);
                      if (part_cell_topo_data)
                        {
                          shards::CellTopology part_cell_topo(part_cell_topo_data);
                          if (part_cell_topo.getKey() == elementType)
                            {
                              selector = selector | *part;
                            }
                        }
                    }
                  selector = selector & (stk::mesh::Selector(m_eMesh.get_fem_meta_data()->locally_owned_part()));
                }
              std::vector<size_t> count ;
              stk::mesh::count_entities( selector, *m_eMesh.get_bulk_data(), count );
              if (count.size() < 3)
                {
                  throw std::logic_error("logic error in Refiner m_refinementInfoByType");
                }
              unsigned n_ele = count[ ranks[irank] ];

              m_refinementInfo.m_refinementInfoByType[irank].m_rank = ranks[irank];
              m_refinementInfo.m_refinementInfoByType[irank].m_numOrigElems = n_ele;

              m_refinementInfo.m_refinementInfoByType[irank].m_numNewElems = n_ele * m_breakPattern[irank]->getNumNewElemPerElem();
              m_refinementInfo.m_refinementInfoByType[irank].m_topology = cell_topo;
              m_refinementInfo.m_refinementInfoByType[irank].m_numOrigNodes = num_orig_nodes;
              m_refinementInfo.m_refinementInfoByType[irank].m_numNewNodes = 0; // can't predict this
            }

          // sum info from all procs
          {
            stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();

            for (unsigned irank = 0; irank < ranks.size(); irank++)
              {
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfo.m_refinementInfoByType[irank].m_numOrigElems ) );
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfo.m_refinementInfoByType[irank].m_numNewElems ) );
                stk::all_reduce( pm, stk::ReduceSum<1>( &m_refinementInfo.m_refinementInfoByType[irank].m_numOrigNodes ) );
              }
          }

        }
      m_elementRankTypeInfo.resize(ranks.size());

      fillElementRankTypeInfo(ranks);
    }

    void Refiner::
    doMark(int num_registration_loops)
    {
#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_START_INSTRUMENTATION;
  CALLGRIND_TOGGLE_COLLECT;
#endif

      stk::diag::Timer timerAdapt_("percept::Refine", rootTimer());
      stk::diag::Timer timerDoMark_("percept::DoMark", timerAdapt_);
      stk::diag::TimeBlock timerAdaptBlock_(timerAdapt_);
      stk::diag::TimeBlock timerDoMarkBlock_(timerDoMark_);

      getRefinementInfo().full_stats_before_mark();

      REF_LTRACE("doMark: start");

      static SubDimCellData empty_SubDimCellData;
      std::vector<stk::mesh::EntityRank>& ranks = m_ranks;

      // do elements first, then any faces or edge elements

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // do the top level, all elements of this rank operation
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      fillElementRankTypeInfo(ranks);

      // FIXME warn if a topology shows up without a break pattern

      ///////////////////////////////////////////////////////////
      /////  // start top-level ranks
      ///////////////////////////////////////////////////////////

      for (int ireg=0; ireg < num_registration_loops; ++ireg) {
        // start top-level ranks

        {
          TIMER2(preMark,DoMark_);
          preMark(ireg, num_registration_loops);
        }

        ///////////////////////////////////////////////////////////
        /////  // node registration step
        ///////////////////////////////////////////////////////////

        {
          // node registration step
          EXCEPTWATCH;
          TIMER2(nodeReg,DoMark_);
          if (0) std::cout << "tmp srk ireg= " << ireg << std::endl;

          /**/  TRACE_PRINT("Refiner: beginRegistration (top-level rank)... ");
          m_nodeRegistry->dumpDB("before init");
          if (m_alwaysInitNodeRegistry)
            {
              m_nodeRegistry->initialize();
            }
          else
            {
              //m_nodeRegistry->clear_element_owner_data();
              m_nodeRegistry->init_entity_repo();
            }

          m_nodeRegistry->init_comm_all();

          m_eMesh.adapt_parent_to_child_relations().clear();

          // register non-ghosted elements needs for new nodes, parallel create new nodes
          m_nodeRegistry->beginRegistration(ireg,num_registration_loops);
          {
            TIMER2(LB_NodeReg,DoMark_);

          //std::cout << "tmp srk INFO: ranks.size()= " << ranks.size() << " m_breakPattern.size()= " << m_breakPattern.size() << std::endl;
          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              {
                EXCEPTWATCH;

                ElementRankTypeInfo& e_info = m_elementRankTypeInfo[irank];
                VERIFY_OP_ON(ranks[irank], ==, e_info.first,"er1");
                VERIFY_OP_ON(elementType, ==, e_info.second,"er2");

                vector<NeededEntityType> needed_entity_ranks;
                m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                bool count_only = false;
                bool doAllElements = true;

                //CHK(m_nodeRegistry);

                unsigned num_elem_not_ghost_0_incr = doForAllElements(irank, "[0/16] Register New Nodes",
                                                                      ranks[irank], &NodeRegistry::registerNeedNewNode,
                                                                      elementType, needed_entity_ranks,
                                                                      count_only, doAllElements);

                //CHK(m_nodeRegistry);

                if (0) std::cout << "tmp srk1 irank= " << irank << " ranks[irank]= " << ranks[irank]
                                 << " nodeRegistry size= " << m_nodeRegistry->getMap().size()
                                 << " num_elem_not_ghost_0_incr= " << num_elem_not_ghost_0_incr
                                 << std::endl;
                //m_nodeRegistry->checkDB(std::string("after doForAllElements rank= " + toString(ranks[irank]) ) );
              }
            }

          }

          if (DO_MEMORY && m_eMesh.get_do_print_memory()) {
            std::string hwm = m_eMesh.print_memory_both();
            if (!m_eMesh.get_rank()) std::cout << "MEM: " << hwm << " before createNewNodesInParallel= " << std::endl;
          }

          {
            TIMER2(NR_EndReg,DoMark_);
            m_nodeRegistry->endRegistration(ireg,num_registration_loops);
          }

          if (DO_MEMORY && m_eMesh.get_do_print_memory()) {
            std::string hwm = m_eMesh.print_memory_both();
            if (!m_eMesh.get_rank()) std::cout << "MEM: " << hwm << " after  createNewNodesInParallel= " << std::endl;
          }

        }
        m_nodeRegistry->dumpDB("after registration");

#define CHECK_DEBUG 0
        if (CHECK_DEBUG)
          {
            MPI_Barrier( MPI_COMM_WORLD );
            std::cout << "P["<< m_eMesh.get_rank()
                      <<"] ========================================================================================================================" << std::endl;
            m_nodeRegistry->checkDB("after registerNeedNewNode");
            check_db("after registerNeedNewNode");
            MPI_Barrier( MPI_COMM_WORLD );
            std::cout << "P["<< m_eMesh.get_rank()
                      <<"] ========================================================================================================================" << std::endl;
          }

        ///////////////////////////////////////////////////////////
        /////  Check for remote
        ///////////////////////////////////////////////////////////

        {   // beginCheckForRemote()
          EXCEPTWATCH;
          TIMER2(checkForRemote,DoMark_);

          /**/                                                        TRACE_PRINT("Refiner: beginCheckForRemote (top-level rank)... ");

          // now register ghosted elements needs for new nodes (this does a pack operation)
          m_nodeRegistry->beginCheckForRemote(ireg,num_registration_loops);
          {
            TIMER2(LB_CheckForRem,DoMark_);
            //unsigned num_elem = 0;
          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
              //if (ranks[irank] >= m_eMesh.face_rank())
              {
                EXCEPTWATCH;

                vector<NeededEntityType> needed_entity_ranks;
                m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                bool count_only = false;
                bool doAllElements = false;  // only do ghost elements if false
                if (!m_nodeRegistry->s_use_new_ownership_check)
                  {
                    doAllElements = false;
                  }
                else
                  {
                    doAllElements = true;
                  }
                doForAllElements(irank, "[1/16] Check For Comm", ranks[irank], &NodeRegistry::checkForRemote, elementType, needed_entity_ranks, count_only, doAllElements);
              }
            }
          }
          m_nodeRegistry->endCheckForRemote(ireg,num_registration_loops);                /**/   TRACE_PRINT("Refiner: endCheckForRemote (top-level rank)... ");

          if (CHECK_DEBUG)
            {
              std::cout << "P["<< m_eMesh.get_rank()
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB("after checkForRemote");
              check_db("after checkForRemote");
              std::cout << "P["<< m_eMesh.get_rank()
                        <<"] ========================================================================================================================" << std::endl;
            }

        }

        ///////////////////////////////////////////////////////////
        /////  Get from remote
        ///////////////////////////////////////////////////////////
        /// communicate all-to-all the new node creation information which also updates the node registry so it can
        /// be queried locally now for any ghost or non-ghost element

        { // get from remote

          {
            TIMER2(beginAndGetFromRemote,DoMark_);

            m_nodeRegistry->beginGetFromRemote(ireg,num_registration_loops);
            {
              TIMER2(LB_GetFromRem,DoMark_);
              //unsigned num_elem = 0;
              for (unsigned irank = 0; irank < ranks.size(); irank++)
                {
                  unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
                  {
                    EXCEPTWATCH;

                    vector<NeededEntityType> needed_entity_ranks;
                    m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

                    bool count_only = false;

                    bool doAllElements = false;   // if false, ghost elements only
                    if (!m_nodeRegistry->s_use_new_ownership_check)
                      {
                        doAllElements = false;
                      }
                    else
                      {
                        doAllElements = true;   // if false, ghost elements only
                      }
                    doForAllElements(irank, "[2/16] Get From Remote", ranks[irank], &NodeRegistry::getFromRemote,  elementType, needed_entity_ranks,  count_only, doAllElements);
                  }
                }
            }
          }

          {
            TIMER2(endGetFromRemote, DoMark_);
            m_nodeRegistry->endGetFromRemote(ireg,num_registration_loops);
          }

          m_nodeRegistry->dumpDB("after endGetFromRemote");

          //stk::diag::printTimersTable(std::cout, perceptTimer(), stk::diag::METRICS_ALL, false);

          if (CHECK_DEBUG)
            {
              MPI_Barrier( MPI_COMM_WORLD );
              std::cout << "P["<< m_eMesh.get_rank()
                        <<"] ========================================================================================================================" << std::endl;
              m_nodeRegistry->checkDB("end getFromRemote");
              check_db("end getFromRemote");
              MPI_Barrier( MPI_COMM_WORLD );

              std::cout << "P["<< m_eMesh.get_rank()
                        <<"] ========================================================================================================================" << std::endl;
            }
        }  // get from remote

        bool break_loop = false;
        {
          TIMER2(postMark, DoMark_);
          break_loop = postMark(ireg, num_registration_loops);
        }

        {
          mod_begin(&timerDoMark_);
          mod_end(&timerDoMark_,"Mark");
        }

        if (break_loop)
          break;

      } // start top-level ranks - num_registration_loops end loop

      getRefinementInfo().full_stats_after_mark();

      REF_LTRACE("doMark: end");
    } // doMark

    void Refiner::
    fillElementRankTypeInfo(std::vector<stk::mesh::EntityRank>& ranks)
    {
      m_elementRankTypeInfo.resize(ranks.size());

      for (unsigned irank = 0; irank < ranks.size(); irank++)
        {
          EXCEPTWATCH;
          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
          shards::CellTopology cell_topo(m_breakPattern[irank]->getFromTopology());

          if (TRACE_STAGE_PRINT) std::cout << "tmp Refiner:: irank = " << irank << " ranks[irank] = " << ranks[irank]
                                           << " elementType= " << elementType
                                           << " cell_topo= " << cell_topo.getName()
                                           << std::endl;

          m_elementRankTypeInfo[irank] = ElementRankTypeInfo(ranks[irank], elementType);
        }
    }

    void Refiner::
    doRefine()
    {
#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_START_INSTRUMENTATION;
  CALLGRIND_TOGGLE_COLLECT;
#endif
      stk::diag::Timer timerAdapt_("percept::Refine", rootTimer());
      stk::diag::TimeBlock timerAdaptBlock_(timerAdapt_);
      stk::diag::Timer timerDoRefine_("percept::DoRefine", timerAdapt_);
      stk::diag::TimeBlock timerDoRefineBlock_(timerDoRefine_);

      getRefinementInfo().full_stats_before_refine();

      m_eMesh.m_nodeRegistry = (void *) &getNodeRegistry();
      REF_LTRACE("doRefine: start");

      stk::mesh::BulkData& bulkData = *m_eMesh.get_bulk_data();
      static SubDimCellData empty_SubDimCellData;
      std::vector<stk::mesh::EntityRank>& ranks = m_ranks;

      m_eMesh.initializeIdServer();

      FixSideSetsSelectorRefine fss_ref(m_eMesh);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // for each element type, in top-down rank order, do the rest of the refinement operations
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
        mod_begin(&timerDoRefine_);
      }

      if (CHECK_DEBUG)
        {
          m_nodeRegistry->checkDB("after mod begin() after end getFromRemote");
        }

      size_t num_new_elements[2] = {0,0};
      for (unsigned irank = 0; irank < ranks.size(); irank++)
        {
          EXCEPTWATCH;

          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
          if (TRACE_STAGE_PRINT) {
            auto breakPattern = m_breakPattern[irank];
            std::cout << "tmp Refiner:: irank = " << irank
                      << " ranks[irank] = " << ranks[irank] << " "
                      << PerceptMesh::demangle(typeid(*breakPattern).name())
                      << " elementType= " << elementType
                      << " fromTopo= " << m_breakPattern[irank]->getFromTopology()->name
                      << " toTopo= " << m_breakPattern[irank]->getToTopology()->name
                      << std::endl;
          }

          std::vector<stk::mesh::EntityRank> ranks_one(1, ranks[irank]);

          ElementRankTypeInfo& e_info = m_elementRankTypeInfo[irank];
          VERIFY_OP_ON(elementType, ==, e_info.second, "mismatch");

          // loop over elements, build faces, edges in threaded mode (guaranteed no mem conflicts)
          // (note: invoke UniformRefinerPattern: what entities are needed)
          vector<NeededEntityType> needed_entity_ranks;
          m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);
          vector<stk::mesh::Entity> new_elements, ft_new_elements;

          {
            EXCEPTWATCH;

            // count num new elements needed on this proc (served by UniformRefinerPattern)
            bool count_only = true;
            bool doAllElements = true;
            unsigned num_elem_not_ghost = 0;
            {
              TIMER2(LB_RegNewNodes,DoRefine_);
              TIMER2(RegNewNodes,DoRefine_);
              /**/                                                TRACE_PRINT("Refiner: registerNeedNewNode count_only(true) ranks[irank]==ranks[0]... ");
              num_elem_not_ghost = doForAllElements(irank, "[3/16] Register New Nodes [Count only]", ranks[irank], &NodeRegistry::registerNeedNewNode,  elementType, needed_entity_ranks, count_only, doAllElements);
              /**/                                                TRACE_PRINT("Refiner: registerNeedNewNode count_only(true) ranks[irank]==ranks[0]... done ");
            }

            //CHK(m_nodeRegistry);

            //unsigned num_elem_needed_old = num_elem_not_ghost * m_breakPattern[irank]->getNumNewElemPerElem();
            unsigned num_elem_needed = 0;
            {
              TIMER2(EstimateNumNewElems,DoRefine_);
              getRefinementInfo().full_stats_before_estimate(ranks[irank], num_elem_not_ghost);
              num_elem_needed = m_breakPattern[irank]->estimateNumberOfNewElements(m_eMesh, ranks[irank], *m_nodeRegistry, num_elem_not_ghost);
              unsigned nelNeedMarked = num_elem_needed / m_breakPattern[irank]->getNumNewElemPerElem();
              getRefinementInfo().full_stats_after_estimate(ranks[irank], nelNeedMarked);
              num_new_elements[0] += num_elem_needed;
            }

            //if (0)
            //  {
            //    std::cout << "P[" << m_eMesh.get_rank() << "] num_new_elements_old, new= " << num_elem_needed_old << " " << num_elem_needed
            //              << " rank= " << ranks[irank]
            //              << " m_breakPattern= " << PerceptMesh::demangle(typeid(*m_breakPattern[irank]).name())
            //              << std::endl;
            //  }

            // FIXME TMP
#define DEBUG_1 0
            if (DEBUG_1 | CHECK_DEBUG)
              {
                std::cout << "tmp Refiner::doBreak: irank= " << irank << " ranks[irank]= " << ranks[irank] << " bp= ["
                          << m_breakPattern[irank]->getFromTopoPartName() << " ==> "
                          << m_breakPattern[irank]->getToTopoPartName() << "] num_elem_needed= " << num_elem_needed
                          << " num_elem_not_ghost= " << num_elem_not_ghost
                          << std::endl;
              }

            // create new entities on this proc
            {
              TIMER2(CreateEntities,DoRefine_);
              m_nodeRegistry->beginLocalMeshMods();
              new_elements.resize(0);                                                /**/ TRACE_PRINT("Refiner: createEntities... ranks[irank]==ranks[0] ");
              ft_new_elements.resize(0);
              auto breakPattern = m_breakPattern[irank];
              std::string bpName = PerceptMesh::demangle(typeid(*breakPattern).name());

              size_t hwm_max=0, hwm_min=0, hwm_avg=0, hwm_sum=0;
              size_t dhwm_max=0, dhwm_min=0, dhwm_avg=0, dhwm_sum=0;

              // if (DO_MEMORY && m_eMesh.get_do_print_memory()) {
              //   m_eMesh.get_memory_high_water_mark_across_processors(m_eMesh.parallel(), hwm_max, hwm_min, hwm_avg, hwm_sum);
              //   std::string hwm = m_eMesh.print_memory_both();
              //   if (!m_eMesh.get_rank()) std::cout << "MEM: " << hwm << " before createEntities= " << bpName << std::endl;
              // }

              if (UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE && ranks[irank] == m_eMesh.side_rank())
                {
                  new_elements.resize(0);
                }
              else
                {
#if USE_CREATE_ENTITIES
                  m_eMesh.createEntities( ranks[irank], num_elem_needed, new_elements);
#else
                  m_eMesh.getEntitiesUsingIdServer( ranks[irank], num_elem_needed, new_elements);
#endif
                }

              if (DO_MEMORY && m_eMesh.get_do_print_memory()) {
                m_eMesh.get_memory_high_water_mark_across_processors(m_eMesh.parallel(), dhwm_max, dhwm_min, dhwm_avg, dhwm_sum);
                dhwm_max -= hwm_max;
                dhwm_min -= hwm_min;
                dhwm_avg -= hwm_avg;
                dhwm_sum -= hwm_sum;
                size_t dhwm_tot = dhwm_sum;
                std::string hwm = m_eMesh.print_memory_both();
                size_t num_elem_needed_tot = num_elem_needed;
                stk::all_reduce( m_eMesh.parallel(), stk::ReduceSum<1>( &num_elem_needed_tot ) );

                if (!m_eMesh.get_rank() && num_elem_needed>0) {
                  std::cout << "MEM: " << hwm << " after  createEntities= " << bpName
                            << "\nMEM: " << double(dhwm_tot)/double(std::max(num_elem_needed_tot,size_t(1))) << " = memory per Entity for "
                            << num_elem_needed_tot << " " << ranks[irank] << " Topo: " << m_breakPattern[irank]->getFromTopoPartName()
                            << std::endl;
                }
              }

              if (0) std::cout << "Refiner::createEntities new_elements.size= " << new_elements.size() << " num_elem_needed= " << num_elem_needed << std::endl;
              if (1)
                {
                  const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
#if USE_CREATE_ENTITIES
                  m_eMesh.createEntities( FAMILY_TREE_RANK, num_elem_needed, ft_new_elements);  /**/ TRACE_PRINT("Refiner: createEntities... ranks[irank]==ranks[0] done ");
#else
                  m_eMesh.getEntitiesUsingIdServer( FAMILY_TREE_RANK, num_elem_needed, ft_new_elements);
#endif
                }
              // if (DO_MEMORY && m_eMesh.get_do_print_memory()) {
              //   std::string hwm = m_eMesh.print_memory_both();
              //   if (!m_eMesh.get_rank()) std::cout << "MEM: " << hwm << " after  FT createEntities= " << bpName << std::endl;
              // }
              m_nodeRegistry->endLocalMeshMods();
            }

          }
          m_nodeRegistry->dumpDB("after endLocalMeshMods");
          if (CHECK_DEBUG)
            {
              m_nodeRegistry->checkDB("after endLocalMeshMods");
            }

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global element ops: here's where we e.g. connect the new elements by declaring new relations
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal... ");
          /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL);

          //CHK(m_nodeRegistry);

          size_t num_actual_new_elems = 0;
          {
            TIMER2(LB_CreateElements,DoRefine_);
            vector<stk::mesh::Entity>::iterator new_elements_pool_end_iter;
            vector<stk::mesh::Entity>::iterator ft_new_elements_pool_end_iter;

            num_actual_new_elems =
              createElementsAndNodesAndConnectLocal(irank, ranks[irank], m_breakPattern[irank], e_info.second, needed_entity_ranks, new_elements, ft_new_elements,
                                                    &new_elements_pool_end_iter,
                                                    &ft_new_elements_pool_end_iter);

            if (ranks[irank] == m_eMesh.element_rank())
              {
                fss_ref.add_elements(new_elements.begin(), new_elements_pool_end_iter);
              }
            num_new_elements[1] += num_actual_new_elems;
          }

          //CHK(m_nodeRegistry);
          /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL);
          /**/                                                TRACE_PRINT("Refiner: createElementsAndNodesAndConnectLocal...done ");

          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          ///  Global node loop operations:  this is where we perform ops like adding new nodes to the right parts, interpolating fields, etc.
          /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

          if (TRACE_STAGE_PRINT && !m_eMesh.get_rank()) {
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL, "CONNECT_LOCAL");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_createNewNeededNodes, "CONNECT_LOCAL_createNewNeededNodes");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_createNewElements, "CONNECT_LOCAL_createNewElements");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_URP_createOrGetNode, "CONNECT_LOCAL_URP_createOrGetNode");
            Util::trace_cpu_time_and_mem_print(CONNECT_LOCAL_URP_declare_relation, "CONNECT_LOCAL_URP_declare_relation");
          }

        } // irank

      if (0)
        {
          stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();
          stk::all_reduce( pm, stk::ReduceSum<2>( &num_new_elements[0] ) );
          std::cout << "P[" << m_eMesh.get_rank() << "] num_new_elements_created, needed= "
                    << num_new_elements[0] << " " << num_new_elements[1] << std::endl;
        }

      /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.]... ");
#if !STK_ADAPT_URP_LOCAL_NODE_COMPS
      if (1)
        {
          EXCEPTWATCH;
          // only need to do this once: the map is fully built and we loop over the map's faces/edges, which are fixed after the getFromRemote step

          doProgressPrint("Stage: [5/16] Add to existing parts...");

          {
            TIMER2(AddToParts,DoRefine_);
            m_nodeRegistry->addToExistingPartsNew();
          }
          //std::cout << "tmp prolongate... " << std::endl;
#if CHECK_DEBUG
          //check_db("after doBreak");
          m_nodeRegistry->checkDB("before prolongate");
#endif

          doProgressPrint("Stage: [6/16] Prolongate fields...");

          {
            TIMER2(ProlongCoord,DoRefine_);
            m_nodeRegistry->prolongate(m_eMesh.get_coordinates_field());
          }
          //std::cout << "tmp prolongate...done " << std::endl;
          //std::cout << "tmp prolongateFields... " << std::endl;

          {
            TIMER2(ProlongFields,DoRefine_);
            m_nodeRegistry->prolongateFields();
          }

          //std::cout << "tmp prolongateFields...done " << std::endl;
#if defined(STK_BUILT_IN_SIERRA)
          if (m_rbar_names.size())
            m_nodeRegistry->add_rbars(m_rbar_names);
#endif
        }
      //std::cout << "tmp dump_elements 1" << std::endl;
      // m_eMesh.dump_elements();
#endif

      if (0)
        {
          stk::mesh::fixup_ghosted_to_shared_nodes(bulkData);
          bulkData.modification_end();

          //CHK(m_nodeRegistry);
          m_eMesh.setProperty("dump_vtk:sides_only", "true");
          m_eMesh.dump_vtk("after-ref-b4-fix-sides.vtk", false, 0, true);
          m_eMesh.setProperty("dump_vtk:sides_only", "false");
          m_eMesh.dump_vtk("after-ref-b4-fix.vtk", false, 0, true);
          //m_eMesh.save_as("after-ref-b4-fix.e");
          bulkData.modification_begin();
        }
      /**/                                                TRACE_PRINT("Refiner: addToExistingParts [etc.] ...done ");

      REF_LTRACE("doRefine: fix_side_sets...");

      doProgressPrint("Stage: [7/16] Remove empty elements...");

      if (1)
      {
        TIMER2(RemEmptyElem,DoRefine_);
        removeEmptyElements();
      }

      bool reduced_mod_end = true;
      if (m_eMesh.getProperty("percept_reduced_mod_end") == "false")
        reduced_mod_end = false;

      /***********************/                           TRACE_PRINT("Refiner: fix_side_sets_2 ");
      {
        doProgressPrint("Stage: [8/16] Fix side sets...");

        mod_end(&timerDoRefine_,"Refine0");
        mod_begin(&timerDoRefine_);

        if (!getIgnoreSideSets())
          {
            TIMER2(FSS_refine,DoRefine_);
            m_eMesh.setProperty("FixSideSets::fix_side_sets_2","refiner");

            //fix_side_sets_2(false,0,0, 0, "Ref1");
            fix_side_sets_2(false,0,0, &fss_ref, "Ref1");
          }
        m_eMesh.adapt_parent_to_child_relations().clear();
      }
      /***********************/                           TRACE_PRINT("Refiner: fix_side_sets_2...done ");

      REF_LTRACE("doRefine: fix_side_sets...done");

      //std::cout << "tmp dump_elements 2" << std::endl;
      //m_eMesh.dump_elements();

#if CHECK_DEBUG
      std::cout << "m_doRemove= " << m_doRemove << std::endl;
      check_db("b4 remove");
#endif

      if (m_doRemove)
        {
          doProgressPrint("Stage: [9/16] Remove old elements...");

          EXCEPTWATCH;

          TIMER2(DoRemove,DoRefine_);

#if CHECK_DEBUG
          //check_db("after doBreak");
          m_nodeRegistry->checkDB("before removeOldElements");
#endif

          for (unsigned irank = 0; irank < ranks.size(); irank++)
            {
              //std::cout << "removeOldElements irank= " << irank << " ranks[]= " << ranks[irank] << std::endl;
#if PERCEPT_USE_FAMILY_TREE
              if (irank == 0)
                removeFamilyTrees();
#endif
              removeOldElements(irank, ranks[irank], m_breakPattern[irank]);
              renameNewParts(ranks[irank], m_breakPattern[irank]);
              fixSurfaceAndEdgeSetNames(ranks[irank], m_breakPattern[irank]);
            }
        }


      // remove any elements that are empty (these can exist when doing local refinement)
#if CHECK_DEBUG
      //check_db("after doBreak");
      m_nodeRegistry->checkDB("before removeEmptyElements");
#endif
      {
        doProgressPrint("Stage: [10/16] Remove empty elements part 2...");

        TIMER2(RemEmptyElem,DoRefine_);
        removeEmptyElements();
      }

#if CHECK_DEBUG
      //check_db("after doBreak");
      m_nodeRegistry->checkDB("after removeEmptyElements");
#endif

      if (m_removeFromNewNodesPart)
        {
          TIMER2(RemoveFromNewNodesPart,DoRefine_);
          removeFromNewNodesPart();
        }

      doProgressPrint("Stage: [11/16] Modification_end...");
      {
        mod_end(&timerDoRefine_,"Refine1");
        mod_begin(&timerDoRefine_);
      }

      doProgressPrint("Stage: [12/16] Remove unattached nodes...");

      // remove nodes not referred to by elements
      {
        TIMER2(RemDangling,DoRefine_);
        removeDanglingNodes();
      }

      REF_LTRACE("doRefine: removeDanglingNodes...done");

      doProgressPrint("Stage: [13/16] Remove empty family trees...");

      {
        TIMER2(RemEmptyFT,DoRefine_);
        removeEmptyFamilyTrees();
      }

#if CHECK_DEBUG
      //check_db("after doBreak");
      m_nodeRegistry->checkDB("after removeDanglingNodes");
#endif

      set_active_part();
      if (m_doAddChildrenToParts)
        add_children_to_parts();

      reset_family_tree_to_node_relations();

      doProgressPrint("Stage: [14/16] Modification end # 2...");

      /**/                                                TRACE_PRINT("Refiner: mod_end...start... ");
      if (!reduced_mod_end)
        {
          // force a flush of all pending deletes, etc
          mod_end(&timerDoRefine_,"Flush1");
          mod_begin(&timerDoRefine_);
          mod_end(&timerDoRefine_,"Flush2");
        }
      else
        {
          mod_end(&timerDoRefine_,"Refine3");
        }

      /**/                                                TRACE_PRINT("Refiner: mod_end...done ");

      if (m_eMesh.get_create_edges())
        {
          if (DEBUG_1) std::cout << "Adapt: create_missing_edges... start" << std::endl;
          RefinerUtil::create_missing_edges(m_eMesh);
          if (DEBUG_1) std::cout << "Adapt: create_missing_edges... end" << std::endl;
        }
      //std::cout << "tmp dump_elements 3" << std::endl;
      //m_eMesh.dump_elements();

      doProgressPrint("Stage: [15/16] Checks on ownership...");

      if (1)
      {
        TIMER2(ref_check_own,DoRefine_);
        check_parent_ownership();
        check_sides_on_same_proc_as_owned_element("doRefine:b4setPEF", true);
        //require_sides_on_same_proc_as_pos_perm_element(true);
      }

      m_eMesh.set_parent_element_field();
      if (1)
        {
          RefinerUtil::save_node_registry(m_eMesh, *m_nodeRegistry, "Refiner: end");
          if (m_eMesh.getProperty("NodeRegistry_rebuild_test") == "true")
            {
              static int cnt = 0;
              ++cnt;
              std::string file1 = "node_reg.1."+toString(cnt)+".e";
              m_eMesh.save_as(file1);
              PerceptMesh eMeshNR;
              eMeshNR.set_ioss_read_options("large");
              eMeshNR.open_read_only(file1);
              NodeRegistry nr1(eMeshNR);
              RefinerUtil::rebuild_node_registry(eMeshNR, nr1, true, &m_eMesh, m_nodeRegistry, true);
            }
        }

      m_eMesh.setProperty("AdaptedMeshVerifier::checkPolarity", "doRefine end");
      AdaptedMeshVerifier::checkPolarity(m_eMesh);
      //check_sidesets_2("end of doBreak");

#if  defined(STK_PERCEPT_HAS_GEOMETRY)
      bool use_ref_mesh = true;
      if (m_eMesh.getProperty("smooth_use_reference_mesh") == "0")
        use_ref_mesh = false;
      snapAndSmooth(m_geomSnap, m_geomFile, use_ref_mesh);
#endif

      REF_LTRACE("doRefine: doRebalance...");

      doProgressPrint("Stage: [16/16] Rebalance...");

      /**/                                                TRACE_PRINT( "Refiner:doBreak ... done");
      {
        TIMER2(Rebalance,DoRefine_);
        doRebalance();
      }
      //std::cout << "tmp m_nodeRegistry.m_gee_cnt= " << m_nodeRegistry->m_gee_cnt << std::endl;
      //std::cout << "tmp m_nodeRegistry.m_gen_cnt= " << m_nodeRegistry->m_gen_cnt << std::endl;
      getRefinementInfo().countCurrentNodes(m_eMesh);
      getRefinementInfo().full_stats_after_refine();

      getNodeRegistry().init_entity_repo();
      //getNodeRegistry().clear_dangling_elements();

      m_nodeRegistry->dumpDB("after doBreak");
#if CHECK_DEBUG
      //check_db("after doBreak");
      m_nodeRegistry->checkDB("after doBreak");
#endif

      doProgressPrint("Stage: [16/16] Rebalance...done, Refinement done.");

#if USE_PERCEPT_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

    } // doRefine

    void Refiner::doRebalance()
    {
      if (m_doRebalance)
        {
          if (m_eMesh.get_rank() == 0)
            {
              std::cout << "Refiner:: rebalancing... weights field= " << m_eMesh.m_weights_field << std::endl;
            }
          if (1)
            {
              RefinerUtil::save_node_registry(m_eMesh, *m_nodeRegistry, "doRebalance");

              if (m_eMesh.getProperty("NodeRegistry_rebuild_test_simple") == "true")
                {
                  PerceptMesh& eMeshNR = m_eMesh;
                  NodeRegistry nr1(eMeshNR);
                  RefinerUtil::rebuild_node_registry(eMeshNR, nr1, true, &m_eMesh, m_nodeRegistry, true);
                  std::cout << "success" << std::endl;
                }

              if (m_eMesh.getProperty("NodeRegistry_rebuild_test") == "true")
                {
                  static int cnt = 0;
                  ++cnt;
                  std::string file1 = "node_reg.3."+toString(cnt)+".e";
                  m_eMesh.save_as(file1);
                  PerceptMesh eMeshNR;
                  eMeshNR.set_ioss_read_options("large");
                  eMeshNR.open_read_only(file1);
                  NodeRegistry nr1(eMeshNR);
                  RefinerUtil::rebuild_node_registry(eMeshNR, nr1, true, &m_eMesh, m_nodeRegistry, true);
                }
            }

          RebalanceMesh rb(m_eMesh, m_eMesh.m_weights_field, false);
          const double imb_before = rb.compute_imbalance();
          if (imb_before > m_rebalThreshold)
            {
              const double imb_after = rb.rebalance();
              if (m_eMesh.get_rank() == 0)
                {
                  std::cout << "Refiner:: imbalance before= " << imb_before << " imbalance after= " << imb_after << std::endl;
                }
              m_eMesh.m_markNone = true;
              m_nodeRegistry->setCheckForGhostedNodes(true);
              initializeDB();
              m_nodeRegistry->setCheckForGhostedNodes(false);
              m_eMesh.m_markNone = false;

              if (1)
                {
                  std::vector< const stk::mesh::FieldBase *> fields;
                  fields.push_back(m_eMesh.m_refine_field);
                  fields.push_back(m_eMesh.m_transition_element_field);
                  stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
                }
            }
          else
            {
              if (m_eMesh.get_rank() == 0)
                {
                  std::cout << "Refiner:: no rebalance done, imbalance before= " << imb_before << " which is < threshold =  " << m_rebalThreshold << std::endl;
                }
            }
        }
    }

    void Refiner::
    finalizeRefine()
    {
    }

    /// Delete all elements that aren't child elements
    void Refiner::deleteParentElements()
    {
      doProgressPrint("Stage: Deleting parent elements...");

      //check_sidesets_2(" deleteParentElements:: start");
      //check_sidesets(" deleteParentElements:: start");
      //check_sidesets_1(" deleteParentElements:: start");

      std::vector<stk::mesh::EntityRank> ranks_to_be_deleted;
      ranks_to_be_deleted.push_back(stk::topology::ELEMENT_RANK);
      ranks_to_be_deleted.push_back(m_eMesh.side_rank());
      if (m_eMesh.get_spatial_dim() == 3)
        ranks_to_be_deleted.push_back(m_eMesh.edge_rank());

      //std::cout << "tmp srk ranks_to_be_deleted= " << ranks_to_be_deleted << std::endl;

      elements_to_be_destroyed_type parents(*m_eMesh.get_bulk_data());
      for (unsigned irank=0; irank < ranks_to_be_deleted.size(); irank++)
        {

          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( ranks_to_be_deleted[irank] );
          int npar=0;
          int nchild=0;
          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              // only do "old" elements
              //if (!oldPartSelector(bucket))
              //  continue;

              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (!m_eMesh.isParentElement(element, false))
                  //if (!m_eMesh.hasFamilyTree(element) || m_eMesh.isChildElement(element, true))
                  //if (!m_eMesh.hasFamilyTree(element) || m_eMesh.isChildElementLeaf(element, true))
                    {
                      // it has no family tree, so it's a leaf, or it has no children
                      ++nchild;
                    }
                  else
                    {
                      ++npar;
                      parents.insert(element);
                    }
                }
            }
          //std::cout << "tmp removeElements(parents) irank, size= " << ranks_to_be_deleted[irank] << " " << npar << " nchild= " << nchild << std::endl;

        }

      mod_begin();

      removeFamilyTrees();

      removeElements(parents);

      fix_side_sets_2(false,0,0,0,"Ref2");

      mod_end(0, "RefinerDelPar");

      doProgressPrint("Stage: Deleting parent elements, number of parent elements = ");

    }

    void Refiner::removeEmptyElements()
    {
      if (CHECK_DEBUG) std::cout << "removeEmptyElements start.... " << std::endl;

      for (unsigned irank=0; irank < m_ranks.size(); irank++)
        {
          elements_to_be_destroyed_type list(*m_eMesh.get_bulk_data());

          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_ranks[irank] );

          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (0 == m_eMesh.get_bulk_data()->num_connectivity(element, stk::topology::NODE_RANK))
                    {
                      list.insert(element);
                    }
                }
            }

          removeElements(list);
        }
      if (CHECK_DEBUG) std::cout << "removeEmptyElements ....end " << std::endl;

    }

    void Refiner::removeEmptyFamilyTrees()
    {
      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);

      elements_to_be_destroyed_type list(*m_eMesh.get_bulk_data());

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( FAMILY_TREE_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity element = bucket[iElement];
              if (m_eMesh.get_bulk_data()->has_no_relations(element))
                {
                  list.insert(element);
                }
            }
        }

      for (elements_to_be_destroyed_type::iterator it = list.begin(); it != list.end(); ++it)
        {
          if ( ! m_eMesh.get_bulk_data()->destroy_entity( *it ) )
            {
              throw std::runtime_error("couldn't destroy family tree entity");
            }
        }
    }

    void Refiner::removeDanglingNodes()
    {
      if (m_avoidClearDanglingNodes)
        return;

      SetOfEntities node_list(*m_eMesh.get_bulk_data());
      SetOfEntities pseudos(*m_eMesh.get_bulk_data());

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          stk::mesh::BulkData &mesh = bucket.mesh();

          const unsigned num_nodes_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_nodes_in_bucket; iElement++)
            {
              stk::mesh::Entity node = bucket[iElement];
              size_t num_rels = mesh.count_relations(node);
              if (0) std::cout << "node= " << m_eMesh.identifier(node) << " delete= " << (0==num_rels) << std::endl;
              if (0 == num_rels)
                {
                  node_list.insert(node);
                }
            }
        }

      //if (1 && !m_eMesh.get_rank()) std::cout << "P[" << m_eMesh.get_rank() << "] tmp number of dangling nodes = " << node_list.size() << " and pseudos= " << pseudos.size() << std::endl;
      //!srk
      getNodeRegistry().clear_dangling_nodes(&node_list);

      if (0)
        {
          for (SetOfEntities::iterator itbd = pseudos.begin(); itbd != pseudos.end();  ++itbd)
            {
              stk::mesh::Entity pseudo_p = *itbd;

              if ( ! m_eMesh.get_bulk_data()->destroy_entity( pseudo_p ) )
                {
                  throw std::logic_error("Refiner::removeDanglingNodes couldn't remove pseudo");

                }
            }
        }

      for (SetOfEntities::iterator itbd = node_list.begin(); itbd != node_list.end();  ++itbd)
        {
          stk::mesh::Entity node_p = *itbd;

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( node_p ) )
            {
              throw std::logic_error("Refiner::removeDanglingNodes couldn't remove node");

            }
        }

      // check for any null entities
      //std::cout << "check for any null entities..." << std::endl;
      const stk::mesh::BucketVector & elem_buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = elem_buckets.begin() ; k != elem_buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_nodes_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_nodes_in_bucket; iElement++)
            {
              stk::mesh::Entity elem = bucket[iElement];
              const percept::MyPairIterRelation rels (m_eMesh, elem, m_eMesh.node_rank());
              for (unsigned j=0; j < rels.size(); j++)
                {
                  if (!m_eMesh.is_valid(rels[j].entity())) throw std::runtime_error("bad node in an element");
                }
            }
        }


    }

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#if 0
    static stk::mesh::Selector getNodeWasSnappedSelector(MeshGeometry& mesh_geometry)
    {
      stk::mesh::Selector selector;
      const std::vector<GeometryEvaluator*>& geomEvaluators = mesh_geometry.getGeomEvaluators();

      for (unsigned s=0; s < geomEvaluators.size(); s++)
        {
          selector |= geomEvaluators[s]->mMesh;
        }
      return selector;
    }
#endif
#endif


#if defined( STK_PERCEPT_HAS_GEOMETRY )
    void Refiner::snapAndSmooth(bool geomSnap, std::string geomFile, bool use_ref_mesh)
    {
      //SMOOTHING_OPTIONS option = SNAP_PLUS_SMOOTH;
      SMOOTHING_OPTIONS option = USE_LINE_SEARCH_WITH_MULTIPLE_STATES;

      GeometryKernel *geomKernel = 0;
      if (geomFile.length() != 0.0)
        {
        if (geomFile.find(".3dm") != std::string::npos)
          {
            std::string m2gFile = geomFile.substr(0,geomFile.length()-3) + "m2g";

            struct stat s;
            if (0 == stat(m2gFile.c_str(), &s))
              {
#if HAVE_CUBIT
                geomKernel = new GeometryKernelPGEOM();
#else
                throw std::runtime_error("CUBIT not supported on this platform");
#endif
              }
            else
              {
                geomKernel = new GeometryKernelOpenNURBS();
              }
          }
        else if (geomFile.find(".e") != std::string::npos || geomFile.find(".g") != std::string::npos || geomFile.find(".exo") != std::string::npos)
          {
            geomKernel = new GeometryKernelGregoryPatch(m_eMesh, false);
          }
        else if(geomFile.find(".sat") != std::string::npos)
          {
#ifdef HAVE_ACIS
    	    geomKernel = new GeometryKernelPGEOM();
#else
    	    throw std::runtime_error("ACIS not supported on this platform");
#endif
          }
        else
          {
            VERIFY_MSG("invalid file extension on --input_geometry file \n   "
                       "-- valid extensions are .3dm (OpenNURBS) or .e,.g,.exo \n"
                       "for GregoryPatch Exodus files or .sat for ACIS (assumes \n"
                       "there is also a file with the same name ending in  .m2g) - file= " + geomFile);
          }
        }

      // set to 0.0 for no checks, > 0.0 for a fixed check delta, < 0.0 (e.g. -0.5) to check against local edge length average times this |value|
      double doCheckMovement = 0.0;
      //double doCheckMovement = -1.0;

      // anything exceeding a value > 0.0 will be printed
      double doCheckCPUTime = 0.0;
      //double doCheckCPUTime = 0.1;

      MeshGeometry mesh_geometry(m_eMesh, geomKernel, doCheckMovement, doCheckCPUTime);
      GeometryFactory factory(geomKernel, &mesh_geometry);
      if (geomFile != "") {
        factory.read_file(geomFile, &m_eMesh);
      }

      switch(option) {
      case SNAP_PLUS_SMOOTH:
        {
          mesh_geometry.snap_points_to_geometry(&m_eMesh);
          if (doCheckMovement != 0.0)
            mesh_geometry.print_node_movement_summary();

          if (m_doSmoothGeometry)
            {
              smoothGeometry(&mesh_geometry, 0, option, use_ref_mesh);
              mesh_geometry.snap_points_to_geometry(&m_eMesh);
            }
        }
        break;
      case USE_LINE_SEARCH_WITH_MULTIPLE_STATES:
        {
          //VERIFY_OP_ON(m_eMesh.get_coordinates_field()->number_of_states(), ==, 3, "Must use PerceptMesh::set_num_coordinate_field_states(3) to use new smoothing.");
          stk::mesh::FieldBase *nm1_field = m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1");

          if (m_doSmoothGeometry && nm1_field)
            {
              // make a copy of current non-snapped state (dst,src)
              m_eMesh.copy_field(m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1"), m_eMesh.get_coordinates_field() );
            }

          // do the snap
          if (geomSnap)
            {
              mesh_geometry.snap_points_to_geometry(&m_eMesh);
            }
          if (doCheckMovement != 0.0)
            mesh_geometry.print_node_movement_summary();

          if (m_doSmoothGeometry && nm1_field)
            {
              if (!geomSnap)
                mesh_geometry.pre_process(&m_eMesh);

              // make a copy of current snapped state
              m_eMesh.copy_field(m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_N"), m_eMesh.get_coordinates_field() );

              // reset current state to non-snapped state
              m_eMesh.copy_field(m_eMesh.get_coordinates_field(), m_eMesh.get_field(stk::topology::NODE_RANK, "coordinates_NM1") );

              // option to smooth without geometry
              if (geomFile == "")
                {
                  bool option_fix_all_internal_and_outer_boundary_nodes=getFixAllBlockBoundaries();
                  stk::mesh::Selector boundarySelector;
                  if (option_fix_all_internal_and_outer_boundary_nodes)
                    {
                      stk::mesh::Part *skin_part = m_eMesh.get_skin_part("inner_skin_part", true);
                      boundarySelector = boundarySelector | *skin_part;
                    }
                  else
                    {
                      // build a selector from all surface parts
                      const stk::mesh::PartVector parts = m_eMesh.get_fem_meta_data()->get_parts();
                      for (unsigned ip=0; ip < parts.size(); ip++)
                        {
                          bool stk_auto= stk::mesh::is_auto_declared_part(*parts[ip]);

                          if (stk_auto) continue;
                          unsigned per = parts[ip]->primary_entity_rank();

                          if (per == m_eMesh.side_rank())
                            {
                              std::cout << "INFO::smoothing: freezing points on boundary: " << parts[ip]->name() << std::endl;
                              boundarySelector = boundarySelector | *parts[ip];
                            }
                        }
                    }
		  stk::mesh::Selector owned = stk::mesh::MetaData::get(*m_eMesh.get_bulk_data()).locally_owned_part();
                  boundarySelector = boundarySelector & owned;
                  smoothGeometry(0, &boundarySelector, option, use_ref_mesh);
                }
              else
                smoothGeometry(&mesh_geometry, 0, option, use_ref_mesh);

            }

        }
        break;
      }

      if (geomKernel) delete geomKernel;
    }
#endif

#if  defined(STK_PERCEPT_HAS_GEOMETRY)
    void Refiner::smoothGeometry(MeshGeometry* mesh_geometry, stk::mesh::Selector* selector, SMOOTHING_OPTIONS option, bool use_ref_mesh)
    {
      bool do_smoothing = true;
      if (do_smoothing)
        {
          int  msq_debug             = 2; // 1,2,3 for more debug info
          bool always_smooth         = true;

          switch(option) {
          case SNAP_PLUS_SMOOTH:
            {
              /// deprecated
            }
            break;
          case USE_LINE_SEARCH_WITH_MULTIPLE_STATES:
            {
              // geometry used for classification of fixed/non-fixed nodes
              int niter = 1001;
              double tol = 1.e-4;
              std::string sni = m_eMesh.getProperty("smoother_niter");
              std::string snt = m_eMesh.getProperty("smoother_tol");
              if (sni.length()) niter = boost::lexical_cast<int>(sni);
              if (snt.length()) tol = boost::lexical_cast<double>(snt);

              if (m_eMesh.getProperty("smoother_type") == "Newton")
                {
#if ENABLE_SMOOTHER3
                  percept::ReferenceMeshSmootherNewton pmmpsi(&m_eMesh, selector, mesh_geometry, niter, tol);
                  pmmpsi.m_do_animation = 0;
                  pmmpsi.m_use_ref_mesh = use_ref_mesh;
                  pmmpsi.run( always_smooth, msq_debug);
#endif
                }
              else if (m_eMesh.getProperty("smoother_type") == "algebraic")
                {
                  double drop_off_coeffs[3] = {1,1,1};  // FIXME
                  int nlayers_drop_off = 40;
                  int niter = 1;
                  percept::ReferenceMeshSmootherAlgebraic pmmpsi(&m_eMesh, selector, mesh_geometry, niter, 1.e-4, drop_off_coeffs, nlayers_drop_off);
                  pmmpsi.m_do_animation = 0;
                  pmmpsi.m_use_ref_mesh = use_ref_mesh;
                  pmmpsi.run( always_smooth, msq_debug);
                }
              else
                {
                  percept::ReferenceMeshSmootherConjugateGradientImpl<STKMesh> pmmpsi(&m_eMesh, selector, mesh_geometry, niter, tol);
                  //pmmpsi.m_do_animation = 0;
                  pmmpsi.m_use_ref_mesh = use_ref_mesh;
                  pmmpsi.run( always_smooth, msq_debug);
                }
            }
            break;
          }
          return;
        }
    }
#endif

    unsigned Refiner::
    doForAllElements(unsigned irank, std::string function_info,
                     stk::mesh::EntityRank rank, NodeRegistry::ElementFunctionPrototype function,
                     unsigned elementType,
                     vector<NeededEntityType>& needed_entity_ranks,
                     bool only_count, bool doAllElements)
    //bool only_count=false, bool doAllElements=true)
    {
      EXCEPTWATCH;
      unsigned num_elem = 0;

      int progress_meter_num_total = 0;
      if (m_doProgress && m_eMesh.get_rank() == 0)
        {
          m_doProgress = false;
          progress_meter_num_total = doForAllElements(irank, function_info, rank, function, elementType, needed_entity_ranks,  true, doAllElements);
          m_doProgress = true;
          if (progress_meter_num_total)
            {
              std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
              ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
              notifyObservers(&pd);
            }
        }
      int progress_meter_when_to_post = progress_meter_num_total / m_progress_meter_frequency;
      if (0 == progress_meter_when_to_post)
        progress_meter_when_to_post = 1;
      double d_progress_meter_num_total = progress_meter_num_total;

      stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->universal_part());
      stk::mesh::PartVector * fromParts = &(m_breakPattern[irank]->getFromParts());
      if (fromParts)
        {
          if (DEBUG_1) std::cout << "tmp srk1 rank= " << rank << " fromParts->size()= " << fromParts->size() << std::endl;
          if (0)
            {
              std::string ostr;
              for (unsigned ii=0; ii < fromParts->size(); ++ii)
                {
                  ostr += " " + (*fromParts)[ii]->name();
                }
              std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk1 rank= " << rank << " fromParts->size()= " << fromParts->size() << "\n" << ostr << std::endl;
            }

          selector = stk::mesh::selectUnion(*fromParts);
        }
      // if (m_excludeParts.size())
      //   {
      //     selector = selector & !stk::mesh::selectUnion(m_excludeParts);
      //   }

      size_t nn_ele=0;
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );
      if (0) std::cout << "tmp srk1 rank= " << rank << " buckets.size= " << buckets.size() << std::endl;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          if (selector(bucket))
            {
              const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(bucket);
              shards::CellTopology topo(bucket_cell_topo_data);
              if (topo.getKey() == elementType)
                {
                  unsigned num_elements_in_bucket = bucket.size();
                  nn_ele += num_elements_in_bucket;
                }
            }
        }
      std::vector<stk::mesh::Entity> elements(nn_ele);
      nn_ele = 0;
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          if (selector(bucket))
            {
              const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(bucket);
              shards::CellTopology topo(bucket_cell_topo_data);
              if (topo.getKey() == elementType)
                {
                  unsigned num_elements_in_bucket = bucket.size();
                  for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                    {
                      stk::mesh::Entity element = bucket[iElement];
                      elements[nn_ele++] = element;
                    }
                }
            }
        }

      if (m_refinerSelector)
        {
          const bool stats = false;
          if (stats)
            {
              getRefinementInfo().full_stats_before_filter(rank, elements.size());
            }
          filterUsingRefinerSelector(rank, elements);
          if (stats)
            {
              getRefinementInfo().full_stats_after_filter(rank, elements.size());
            }
        }

      for (size_t iElement=0; iElement < elements.size(); ++iElement)
        {
          stk::mesh::Entity element = elements[iElement];
          const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(element);
          shards::CellTopology topo(bucket_cell_topo_data);

          bool elementIsGhost = m_eMesh.isGhostElement(element);
          if (!elementIsGhost)
            ++num_elem;

          VERIFY_OP_ON(m_eMesh.is_valid(element), ==, true, "doForAllElements bad element");
          if (!only_count && (doAllElements || elementIsGhost))
            {
              refineMethodApply(function, element, needed_entity_ranks, bucket_cell_topo_data);
            }

          if (m_doProgress && m_eMesh.get_rank() == 0 && iElement && (iElement % progress_meter_when_to_post == 0) )
            {
              double progress_meter_percent = 100.0*((double)num_elem)/std::max(d_progress_meter_num_total,1.0);
              std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
              ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, oss.str());
              notifyObservers(&pd);
              if (0) std::cout << "progress_meter_percent = " << progress_meter_percent << std::endl;
            }
        } // elements

      if (m_doProgress && m_eMesh.get_rank() == 0 && progress_meter_num_total)
        {
          std::ostringstream oss; oss << function_info << " [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
          notifyObservers(&pd);
        }

      return num_elem;
    }

    void Refiner::filterUsingRefinerSelector(stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& elements)
    {
      if (m_refinerSelector)
        {
          if (m_refinerSelector->use_batch_filter())
            {
              m_refinerSelector->batch_filter(rank, elements);
            }
          else
            {
              size_t nele=0;
              for (size_t ii=0; ii < elements.size(); ++ii)
                {
                  if ((*m_refinerSelector)(elements[ii]) || m_eMesh.isGhostElement(elements[ii]))
                    {
                      ++nele;
                    }
                }
              //std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk filterUsingRefinerSelector:: orig size= " << elements.size()
              //          << " filtered size= " << nele << std::endl;
              std::vector<stk::mesh::Entity> elements_new(nele);
              nele=0;
              for (size_t ii=0; ii < elements.size(); ++ii)
                {
                  if ((*m_refinerSelector)(elements[ii]) || m_eMesh.isGhostElement(elements[ii]))
                    {
                      elements_new[nele++] = elements[ii];
                    }
                }
              //std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk 1 filterUsingRefinerSelector:: orig size= " << elements.size()
              //          << " filtered size= " << nele << std::endl;
              elements = elements_new;
            }
        }
    }

    size_t Refiner::
    createElementsAndNodesAndConnectLocal(unsigned irank, stk::mesh::EntityRank rank, UniformRefinerPatternBase *breakPattern,
                                          unsigned elementType,   vector<NeededEntityType>& needed_entity_ranks,
                                          vector<stk::mesh::Entity>& new_elements_pool,
                                          vector<stk::mesh::Entity>& ft_element_pool,
                                          vector<stk::mesh::Entity>::iterator * new_elements_pool_end_iter,
                                          vector<stk::mesh::Entity>::iterator * ft_new_elements_pool_end_iter
                                          )
    {
      EXCEPTWATCH;
      static NewSubEntityNodesType s_new_sub_entity_nodes(percept::EntityRankEnd);

      NewSubEntityNodesType& new_sub_entity_nodes = s_new_sub_entity_nodes;

      vector<stk::mesh::Entity>::iterator element_pool_it = new_elements_pool.begin();
      vector<stk::mesh::Entity>::iterator ft_element_pool_it = ft_element_pool.begin();

      int jele = 0;
      int numPrints = 20;

      stk::mesh::Selector selector(m_eMesh.get_fem_meta_data()->universal_part());
      stk::mesh::PartVector * fromParts = &(m_breakPattern[irank]->getFromParts());
      if (fromParts)
        {
          selector = stk::mesh::selectUnion(*fromParts);
        }
      if (m_excludeParts.size())
        {
          selector = selector & !stk::mesh::selectUnion(m_excludeParts);
          if (!m_eMesh.get_rank() && m_eMesh.getProperty("MeshAdapt.debug") == "true")
            std::cout << "Refiner:createElementsAndNodesAndConnectLocal with excluded parts, selector= " << selector << std::endl;
        }

      // create new elements and connect them up

      std::vector<stk::mesh::Entity> elems;

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          if (selector(bucket))
            {
              const CellTopologyData * const bucket_cell_topo_data = m_eMesh.get_cell_topology(bucket);
              shards::CellTopology topo(bucket_cell_topo_data);
              if (topo.getKey() == elementType)
                {
                  unsigned num_elements_in_bucket = bucket.size();
                  jele += num_elements_in_bucket;
                  for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                    {
                      stk::mesh::Entity element = bucket[iElement];
                      elems.push_back(element);
                    }
                }
            }
        }

      int nele = jele;
      jele = 0;
      int printEvery = nele/numPrints;
      if (printEvery == 0) printEvery = 1;
      if (0)
        {
          std::cout << "Refiner::createElementsAndNodesAndConnectLocal: rank= " << rank
                    << " num elements = " << nele
                    << " printEvery= " << printEvery
                    << std::endl;
        }
      if (m_doProgress && m_eMesh.get_rank() == 0 && elems.size())
        {
          std::ostringstream oss; oss << "[4/16] Create Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
          notifyObservers(&pd);
        }

      if (nele == 0)
        {
          *new_elements_pool_end_iter = element_pool_it;
          *ft_new_elements_pool_end_iter = ft_element_pool_it;
          return 0;
        }

      stk::mesh::Entity first_element = elems[0];
      const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(first_element);
      shards::CellTopology cell_topo(cell_topo_data);

      for (int iElement = 0; iElement < nele; iElement++)
        {
          stk::mesh::Entity element = elems[iElement];
          const stk::mesh::Entity element_p = element;

          if (!m_eMesh.is_valid(element_p))
            {
              throw std::runtime_error("Refiner::createElementsAndNodesAndConnectLocal");
            }

          if (0 && (jele % printEvery == 0))
            {
              std::cout << "Refiner::createElementsAndNodesAndConnectLocal: element # = " << jele << " ["
                        << (((double)jele)/((double)nele)*100.0) << " %]" << std::endl;
            }
          if (m_doProgress && m_eMesh.get_rank() == 0 && (jele % printEvery == 0))
            {
              std::ostringstream oss; oss << "[4/16] Create Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
              ProgressMeterData pd(ProgressMeterData::RUNNING, 100.0*((double)jele)/((double)std::max(nele,1)), oss.str());
              notifyObservers(&pd);
            }

          if (m_proc_rank_field && rank == stk::topology::ELEMENT_RANK)
            {
              double *fdata = stk::mesh::field_data( *static_cast<const ScalarFieldType *>(m_proc_rank_field) , element );
              fdata[0] = double(m_eMesh.owner_rank(element));
              //if (1 || eMesh.owner_rank(element) == 3)
              //  std::cout << "tmp eMesh.owner_rank(element) = " << eMesh.owner_rank(element) << std::endl;
            }
          // FIXME

          // skip elements that are already a parent (if there's no family tree yet, it's not a parent, so avoid throwing an error if isParentElement)
          const bool check_for_family_tree = false;
          bool isParent = m_eMesh.isParentElement(element, check_for_family_tree);
          if (0)
            {
              const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
              percept::MyPairIterRelation element_to_family_tree_relations (m_eMesh, element, FAMILY_TREE_RANK);
              if (element_to_family_tree_relations.size() == 1)
                {
                  std::cout << "tmp isParent = " << isParent << " isChild = " << m_eMesh.isChildElement(element) << " element_to_family_tree_relations.size() = " << element_to_family_tree_relations.size() << std::endl;
                }
            }

          if (isParent)
            continue;


          if (!m_eMesh.isGhostElement(element))
            {
              //std::cout << "P["<< m_eMesh.get_rank() << "] eMesh.owner_rank(element) = " << eMesh.owner_rank(element) << std::endl;
              /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL_createNewNeededNodes);

              if (createNewNeededNodeIds(cell_topo_data, element, needed_entity_ranks, new_sub_entity_nodes, breakPattern))
                {
                  std::cout << "typeid= " << typeid(*breakPattern).name() << std::endl;
                  throw std::logic_error("needed_entity_ranks[ineed_ent].second");
                }

              /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL_createNewNeededNodes);

              /**/                                                TRACE_CPU_TIME_AND_MEM_0(CONNECT_LOCAL_createNewElements);

              //if (0)
              //  std::cout << "tmp Refiner: element = " << m_eMesh.identifier(element) << " topo= " << m_eMesh.bucket(element).topology()
              //            << " bp = " << PerceptMesh::demangle(typeid(*breakPattern).name())
              //            << " new_elements_pool.size= " << new_elements_pool.size() << std::endl;

              //breakPattern->m_ep_begin = new_elements_pool.begin();
              //breakPattern->m_ep_end = new_elements_pool.end();
              vector<stk::mesh::Entity>::iterator current_it = element_pool_it;
              breakPattern->createNewElements(m_eMesh, *m_nodeRegistry, element, new_sub_entity_nodes, element_pool_it, ft_element_pool_it, m_proc_rank_field);

              /**/                                                TRACE_CPU_TIME_AND_MEM_1(CONNECT_LOCAL_createNewElements);
              vector<stk::mesh::Entity>::iterator new_it = element_pool_it;
              size_t num_new_elemements = new_it - current_it;
              if (0) std::cout << " tmp Refiner: num_new_elemements= " << num_new_elemements << std::endl;

#ifndef NDEBUG
              if (0)
                {
                  m_nodeRegistry->prolongateCoordsAllSubDims(element);

                  vector<stk::mesh::Entity> newElements(current_it, new_it);
                  newElements.insert(newElements.begin(), element);

                  for (unsigned ii=0; ii < newElements.size(); ++ii)
                    {
                      stk::mesh::Entity newElement = newElements[ii];
                      //if (ii == 0) m_nodeRegistry->prolongateCoords(newElement, stk::topology::ELEMENT_RANK, 0u);
                      VolumeUtil jacA;
                      double jacobian = 0.0;
                      jacA(jacobian, m_eMesh, newElement, m_eMesh.get_coordinates_field());
                      if (jacobian < 0)
                        {
                          if (1) std::cout << "tmppyr Refiner " << m_eMesh.demangle(typeid(*breakPattern).name()) << " createNewElements: id= " << m_eMesh.identifier(newElement) << " jac= " << jacobian << std::endl;
                          m_eMesh.print_entity(newElement);
                          std::cout << " tmp Refiner: num_new_elemements= " << num_new_elemements << std::endl;
                          //std::cout << "tmppyr Refiner negative volume, stacktrace=\n" << eMesh.demangled_stacktrace() << std::endl;
                          std::cout << "tmppyr Refiner parts= " << m_eMesh.print_entity_parts_string(newElement, "\n");
                          std::cout << "tmppyr Refiner s_do_transition_break= " << s_do_transition_break << std::endl;

                          throw std::runtime_error("neg vol Refine");
                        }
                    }
                }
#endif
            }

          ++jele;
        }

      if (!UniformRefinerPatternBase::USE_DECLARE_ELEMENT_SIDE)
        {
          ptrdiff_t ndiff = new_elements_pool.end() - element_pool_it;
          VERIFY_OP_ON(ndiff, >=, 0,  "ERROR: bad element pool, element_pool.size= "+toString(new_elements_pool.size()));
        }
      size_t num_new_elems = element_pool_it - new_elements_pool.begin();

      *new_elements_pool_end_iter = element_pool_it;
      *ft_new_elements_pool_end_iter = ft_element_pool_it;

      if (m_doProgress && m_eMesh.get_rank() == 0 && elems.size())
        {
          std::ostringstream oss; oss << "[4/16] Create Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
          notifyObservers(&pd);
        }
      return num_new_elems;
    }

    /// create a list of nodes from the new nodes that can be easily deciphered by the UniformRefinerPattern
    /// Returns the 3D array new_sub_entity_nodes[entity_rank][ordinal_of_sub_dim_entity][ordinal_of_node_on_sub_dim_entity]

    bool Refiner::
    createNewNeededNodeIds(const CellTopologyData * const cell_topo_data,
                           const stk::mesh::Entity element, vector<NeededEntityType>& needed_entity_ranks, NewSubEntityNodesType& new_sub_entity_nodes, UniformRefinerPatternBase *breakPattern)
    {
      EXCEPTWATCH;

      NodeRegistry& nodeRegistry = *m_nodeRegistry;

      const percept::MyPairIterRelation elem_nodes (m_eMesh, element, stk::topology::NODE_RANK);

      // CHECK - cache this
      for (unsigned ineed_ent=0; ineed_ent < needed_entity_ranks.size(); ineed_ent++)
        {
          unsigned numSubDimNeededEntities = 0;

          // special case of face in 3d or edge in 2d
          if (needed_entity_ranks[ineed_ent].first == m_eMesh.entity_rank(element))
            {
              numSubDimNeededEntities = 1;
            }
          else if (needed_entity_ranks[ineed_ent].first == m_eMesh.edge_rank())
            {
              numSubDimNeededEntities = cell_topo_data->edge_count;
            }
          else if (needed_entity_ranks[ineed_ent].first == m_eMesh.face_rank())
            {
              numSubDimNeededEntities = cell_topo_data->side_count;
            }
          else if (needed_entity_ranks[ineed_ent].first == stk::topology::ELEMENT_RANK)
            {
              numSubDimNeededEntities = 1;
            }

          if (needed_entity_ranks[ineed_ent].first >= new_sub_entity_nodes.size())
            {
              throw std::logic_error("Refiner::createNewNeededNodeIds logic err #1");
            }
          new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first].resize(numSubDimNeededEntities);

          if (0)
            {
              std::cout << "P[" << m_eMesh.get_rank() << "]  needed_entity_ranks[ineed_ent]= " << needed_entity_ranks[ineed_ent].first
                        << " , " << needed_entity_ranks[ineed_ent].second << " numSubDimNeededEntities= " << numSubDimNeededEntities
                        << std::endl;
            }

          // ensure nodes don't get inadvertently reused
          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
            }

          for (unsigned iSubDimOrd = 0; iSubDimOrd < numSubDimNeededEntities; iSubDimOrd++)
            {
              // CHECK
              NodeIdsOnSubDimEntityType* nodeIds_onSE_ptr = nodeRegistry.getNewNodesOnSubDimEntity(element, needed_entity_ranks[ineed_ent].first, iSubDimOrd);
              if (nodeIds_onSE_ptr == 0)
                {
                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
                  continue;
                }
              NodeIdsOnSubDimEntityType& nodeIds_onSE = *nodeIds_onSE_ptr;

              if (nodeIds_onSE.size() == 0)
                {
                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
                  continue;
                }

              if (needed_entity_ranks[ineed_ent].third.size())
                {
                  VERIFY_OP_ON(needed_entity_ranks[ineed_ent].third.size(), ==, numSubDimNeededEntities, "bad size");
                  if (!needed_entity_ranks[ineed_ent].third[iSubDimOrd])
                    {
                      new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(0);
                      continue;
                    }
                }

              if (!m_eMesh.is_valid(nodeIds_onSE[0])) {

                if (nodeIds_onSE.m_entity_id_vector[0] == 0)
                  {
                    // for debugging cases that may have inconsistent edge marking schemes
#define DEBUG_ALLOW_0_ENTITY_ID_VECTOR 0
                    if (DEBUG_ALLOW_0_ENTITY_ID_VECTOR)
                      {
                        continue;
                      }
                    else
                      {
                        SubDimCell_SDCEntityType subDimEntity(m_eMesh);
                        m_nodeRegistry->getSubDimEntity(subDimEntity, element, needed_entity_ranks[ineed_ent].first, iSubDimOrd);
                        std::cout << "P[" << m_eMesh.get_rank() << "] nodeId ## = 0 << "
                                  << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0]
                                  << " element= " << m_eMesh.id(element)
                                  << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                                  << " needed_entity_ranks.third= " << needed_entity_ranks[ineed_ent].third
                                  << " breakPattern= " << m_eMesh.demangle(typeid(*breakPattern).name())
                                  << " iSubDimOrd = " << iSubDimOrd
                                  << " node0 = " << m_eMesh.identifier(subDimEntity[0])
                                  << " node1 = " << (subDimEntity.size() > 1 ? m_eMesh.identifier(subDimEntity[1]) : 0)
                                  << " node2 = " << (subDimEntity.size() > 2 ? m_eMesh.identifier(subDimEntity[2]) : 0)
                                  <<  std::endl;
                        m_eMesh.dump_vtk("bad-5.0.vtk");
                        stk::mesh::Entity neigh = m_eMesh.get_face_neighbor(element, iSubDimOrd);
                        std::cout << " element= " << m_eMesh.print_entity_compact(element) << "\nneigh= " << m_eMesh.print_entity_compact(neigh) << std::endl;
                        std::set<stk::mesh::Entity> ll_debug;
                        ll_debug.insert(element);
                        ll_debug.insert(neigh);
                        std::string ff = "bad-"+toString(m_eMesh.get_rank())+".vtk";
                        m_eMesh.dump_vtk(ff, false, &ll_debug);
                        VERIFY_MSG("Refiner::createNewNeededNodeIds logic err #5.0, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
                      }
                  }

                stk::mesh::Entity node1 = m_eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, nodeIds_onSE.m_entity_id_vector[0]);
                nodeIds_onSE[0] = node1;

                if (!m_eMesh.is_valid(node1))
                  {
                    if (!m_nodeRegistry->getUseCustomGhosting())
                    {
                      static stk::mesh::PartVector empty_parts;
                      node1 = m_eMesh.get_bulk_data()->declare_node(nodeIds_onSE.m_entity_id_vector[0], empty_parts);
                    }

                    if (!m_eMesh.is_valid(node1))
                    {
                      std::cout << "P[" << m_eMesh.get_rank() << "] nodeId ## = 0 << "
                              << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0] << " node1= " << node1
                              << " element= " << element
                              << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                              << " iSubDimOrd = " << iSubDimOrd
                              <<  std::endl;
                      throw std::logic_error("Refiner::createNewNeededNodeIds logic error #0");
                    }
                  }

              }

              unsigned num_new_nodes_needed = needed_entity_ranks[ineed_ent].second;

              if (num_new_nodes_needed < 1)
                {
                  //std::cout << "needed_entity_ranks[ineed_ent].second = " << num_new_nodes_needed << std::endl;
                  //throw std::logic_error("needed_entity_ranks[ineed_ent].second");
                  return true;
                }

              if (iSubDimOrd >= new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first].size())
                {
                  throw std::logic_error("Refiner::createNewNeededNodeIds logic err #2");
                }
              //std::cout << "tmp elementid, iSubDimOrd, num_new_nodes_needed = " << element) << " " << iSubDimOrd << " " << num_new_nodes_needed << std::endl;
              new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd].resize(num_new_nodes_needed);
              if (num_new_nodes_needed > nodeIds_onSE.size())
                {

                  std::cout << "Refiner::createNewNeededNodeIds logic err #3:  num_new_nodes_needed= " << num_new_nodes_needed
                            << " nodeIds_onSE.size() = " << nodeIds_onSE.size()
                            << " needed_entity_ranks.first= " << needed_entity_ranks[ineed_ent].first
                            << " needed_entity_ranks.second= " << needed_entity_ranks[ineed_ent].second
                            << " numSubDimNeededEntities= " << numSubDimNeededEntities
                            << " iSubDimOrd = " << iSubDimOrd
                            << std::endl;
                  std::cout << "P[" << m_eMesh.get_rank() << "] element= ";
                  m_eMesh.print(element);
                  for (unsigned kn=0; kn < elem_nodes.size(); kn++)
                    {
                      m_eMesh.dump_vtk(elem_nodes[kn].entity(), std::string("node_dump_")+toString(kn)+".vtk");
                    }

                  throw std::logic_error("Refiner::createNewNeededNodeIds logic err #3");
                }

              for (unsigned i_new_node = 0; i_new_node < num_new_nodes_needed; i_new_node++)
                {
                  if (!m_eMesh.is_valid(nodeIds_onSE[i_new_node]))
                    {
                      if (nodeIds_onSE.m_entity_id_vector[i_new_node] == 0)
                        {
                          if (DEBUG_ALLOW_0_ENTITY_ID_VECTOR)
                            {
                              continue;
                            }
                          else
                            {
                              SubDimCell_SDCEntityType subDimEntity(m_eMesh);
                              m_nodeRegistry->getSubDimEntity(subDimEntity, element, needed_entity_ranks[ineed_ent].first, iSubDimOrd);
                              std::cout << "P[" << m_eMesh.get_rank() << "] nodeId ## = 0 << "
                                        << " nodeIds_onSE.m_entity_id_vector[0] = " << nodeIds_onSE.m_entity_id_vector[0]
                                        << " element= " << m_eMesh.identifier(element)
                                        << " needed_entity_ranks= " << needed_entity_ranks[ineed_ent].first
                                        << " iSubDimOrd = " << iSubDimOrd
                                //<< " breakPattern= " << typeid(*breakPattern).name()
                                        << " node0 = " << m_eMesh.identifier(subDimEntity[0])
                                        << " node1 = " << m_eMesh.identifier(subDimEntity[1])
                                        <<  std::endl;
                              //m_eMesh.print_entity(std::cout, element, 0);
                              m_eMesh.dump_vtk("bad-5.1.vtk");
                              stk::mesh::Entity neigh = m_eMesh.get_face_neighbor(element, iSubDimOrd);
                              std::cout << " element= " << m_eMesh.print_entity_compact(element) << "\nneigh= " << m_eMesh.print_entity_compact(neigh) << std::endl;
                              std::set<stk::mesh::Entity> ll_debug;
                              ll_debug.insert(element);
                              ll_debug.insert(neigh);
                              std::string ff = "bad-"+toString(m_eMesh.get_rank())+".vtk";
                              m_eMesh.dump_vtk(ff, false, &ll_debug);
                              VERIFY_MSG("Refiner::createNewNeededNodeIds logic err #5.1, nodeIds_onSE.m_entity_id_vector[i_new_node] == 0");
                            }
                        }
                      stk::mesh::Entity node1 = m_eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, nodeIds_onSE.m_entity_id_vector[i_new_node]);

                      if (!m_eMesh.is_valid(node1))
                        {
                          if (!m_nodeRegistry->getUseCustomGhosting())
                          {
                            static stk::mesh::PartVector empty_parts;
                            node1 = m_eMesh.get_bulk_data()->declare_node(nodeIds_onSE.m_entity_id_vector[i_new_node], empty_parts);
                          }
                          if (!m_eMesh.is_valid(node1))
                          {
                            throw std::logic_error("Refiner::createNewNeededNodeIds logic err #4");
                          }
                        }
                      nodeIds_onSE[i_new_node] = node1;
                      VERIFY_OP_ON(m_eMesh.identifier(node1), ==, nodeIds_onSE.m_entity_id_vector[i_new_node], "Refiner::createNewNeededNodeIds logic err #4.1");
                    }
                  new_sub_entity_nodes[needed_entity_ranks[ineed_ent].first][iSubDimOrd][i_new_node] = m_eMesh.identifier(nodeIds_onSE[i_new_node]);

                }
            }
        }
      return false;
    }


    /** Creates a map of element sides to their higher-dimensional base elements
     */

#define EXTRA_PRINT_UR_BESDB 0

    void Refiner::
    buildElementSideDB(SubDimCellToDataMap& cell_2_data_map)
    {

    }


    /** Sets orientations and associativity of elements to sub-dimensional faces/edges after refinement.
     */
#define EXTRA_PRINT_UR_FES 0

    struct FindPair {
      const std::string& m_side_part_name;
      const std::string& m_elem_part_name;
      FindPair(const std::string& side_part_name, const std::string& elem_part_name) : m_side_part_name(side_part_name), m_elem_part_name(elem_part_name) {}

      bool operator()(std::pair<stk::mesh::Part* const, stk::mesh::PartVector >& iter)
      {
        bool side_found =  iter.first->name() == m_side_part_name;
        bool elem_found = false;
        for (unsigned i=0; i < iter.second.size(); i++)
          {
            if (m_elem_part_name == iter.second[i]->name())
              {
                elem_found = true;
                break;
              }
          }
        return side_found && elem_found;
      }
    };

#define DEBUG_GSPR 0

    static void add_if_not_present(stk::mesh::Part *side_part, stk::mesh::Part *elem_part, SidePartMap& side_part_map)
    {
      SidePartMap::iterator found = std::find_if(side_part_map.begin(), side_part_map.end(), FindPair(side_part->name(), elem_part->name()));
      if (found != side_part_map.end())
        return;

      stk::mesh::PartVector& epv_add = side_part_map[side_part];
      epv_add.push_back(elem_part);
      if (DEBUG_GSPR) std::cout << "Refiner::get_side_part_relations: add_to_map = " << side_part->name() << " elem_part= " << elem_part->name() << std::endl;
    }


    // determine side part to elem part relations
    //static
    void Refiner::
    get_side_part_relations(PerceptMesh& eMesh, bool checkParentChild, SidePartMap& side_part_map, bool debug)
    {
      EXCEPTWATCH;

      //stk::mesh::EntityRank node_rank = eMesh.node_rank();
      stk::mesh::EntityRank edge_rank = eMesh.edge_rank();
      stk::mesh::EntityRank side_rank = eMesh.side_rank();
      stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;

      int spatialDim = eMesh.get_spatial_dim();

      stk::mesh::EntityRank side_rank_iter_begin = side_rank;
      stk::mesh::EntityRank side_rank_iter_end = side_rank;
      if (spatialDim == 3)
        {
          side_rank_iter_begin = edge_rank;
        }

      if (debug || DEBUG_GSPR)
        {
          bool doAll = true;
          const stk::mesh::PartVector parts = eMesh.get_fem_meta_data()->get_parts();
          for (unsigned ip=0; ip < parts.size(); ip++)
            {
              bool stk_auto= stk::mesh::is_auto_declared_part(*parts[ip]);
              //const CellTopologyData *const topology = eMesh.get_cell_topology(*parts[ip]);
              if (stk_auto && !doAll)
                continue;
              if (eMesh.get_rank() == 0) std::cout << "tmp srk get_side_part_relations:: parts[ip]-> == " << parts[ip]->name() << std::endl;
            }
        }

      // get super-relations (side_part.name() --> elem_part.name())
      for (stk::mesh::EntityRank side_rank_iter = side_rank_iter_begin; side_rank_iter <= side_rank_iter_end; side_rank_iter++)
        {
          const stk::mesh::BucketVector & side_buckets = eMesh.get_bulk_data()->buckets( side_rank_iter );
          for ( stk::mesh::BucketVector::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
            {
              stk::mesh::Bucket & side_bucket = **it_side_bucket ;
              stk::mesh::PartVector const& side_parts = side_bucket.supersets();

              if (DEBUG_GSPR)
                {
                  std::cout << "\nside_bucket.supersets() =  " << std::endl;
                  for (unsigned isp=0; isp < side_parts.size(); isp++)
                    {
                      bool stk_auto= stk::mesh::is_auto_declared_part(*side_parts[isp]);
                      const CellTopologyData *const topology = eMesh.get_cell_topology(*side_parts[isp]);
                      unsigned per = side_parts[isp]->primary_entity_rank();
                      std::cout << "superset= " << side_parts[isp]->name() << " stk_auto= " << stk_auto << " topology= " << topology << " primary_entity_rank= " << per << std::endl;
                    }
                }

              stk::mesh::PartVector elem_parts = side_parts; // since high-rank parts are in side_parts already
              for (unsigned isp=0; isp < side_parts.size(); isp++)
                {
                  if ( stk::mesh::is_auto_declared_part(*side_parts[isp]) )
                    continue;

                  const AutoPart *side_auto_part = side_parts[isp]->attribute<AutoPart>();
                  if (side_auto_part)
                    continue;

                  unsigned per = side_parts[isp]->primary_entity_rank();
                  if (per != side_rank_iter)
                    continue;

                  for (unsigned iep=0; iep < elem_parts.size(); iep++)
                    {
                      if ( stk::mesh::is_auto_declared_part(*elem_parts[iep]) )
                        continue;

                      const AutoPart *auto_part = elem_parts[iep]->attribute<AutoPart>();
                      if (elem_parts[iep]->name().find(UniformRefinerPatternBase::getOldElementsPartName()) != std::string::npos)
                        {
                          if (!auto_part) throw std::runtime_error("Refiner::get_side_part_relations: bad old part attribute for auto");
                        }

                      if (auto_part)
                        {
                          continue;
                        }

                      if (elem_parts[iep]->primary_entity_rank() == element_rank)
                        {
                          if (elem_parts[iep] != side_parts[isp])
                            {
                              SidePartMap::iterator found = side_part_map.find(side_parts[isp]);
                              if (found == side_part_map.end())
                                {
                                  side_part_map[side_parts[isp]] = stk::mesh::PartVector(1, elem_parts[iep]);
                                }
                              else
                                {
                                  stk::mesh::PartVector& epv = found->second;
                                  stk::mesh::PartVector::iterator epv_found = std::find(epv.begin(), epv.end(), elem_parts[iep]);
                                  if (epv_found == epv.end())
                                    {
                                      epv.push_back(elem_parts[iep]);
                                    }
                                }
                            }
                        }
                    }
                }
            }

          // add parts created by UnrefinerPattern
          const stk::mesh::PartVector parts = eMesh.get_fem_meta_data()->get_parts();
          SidePartMap side_part_map_copy = side_part_map;

          SidePartMap::iterator iter;
          for (iter = side_part_map_copy.begin(); iter != side_part_map_copy.end(); iter++)
            {
              stk::mesh::Part *side_part = iter->first;
              std::string side_part_name = side_part->name();
              const stk::mesh::PartVector *side_pv  = side_part->attribute<stk::mesh::PartVector>();
              const stk::mesh::PartVector& epv = iter->second;
              for (unsigned iepv=0; iepv < epv.size(); iepv++)
                {
                  stk::mesh::Part *elem_part = epv[iepv];
                  std::string elem_part_name = elem_part->name();
                  if (DEBUG_GSPR) std::cout << "looking for other parts: side_part = " << side_part->name() << " elem_part= " << elem_part_name << std::endl;

                  const stk::mesh::PartVector *elem_pv  = elem_part->attribute<stk::mesh::PartVector>();

                  stk::mesh::PartVector side_pv1;
                  if (side_pv) side_pv1 = *side_pv;
                  stk::mesh::PartVector elem_pv1;
                  if (elem_pv) elem_pv1 = *elem_pv;
                  side_pv1.push_back(side_part);
                  elem_pv1.push_back(elem_part);
                  for (unsigned iside_pv = 0; iside_pv < side_pv1.size(); iside_pv++)
                    {
                      for (unsigned ielem_pv = 0; ielem_pv < elem_pv1.size(); ielem_pv++)
                        {
                          stk::mesh::Part *sidep = side_pv1[iside_pv];
                          stk::mesh::Part *elemp = elem_pv1[ielem_pv];
                          add_if_not_present(sidep, elemp, side_part_map);
                        }
                    }
                }
            }
        }

      if ((debug || DEBUG_GSPR) && !eMesh.get_rank())
        {
          SidePartMap::iterator iter;
          for (iter = side_part_map.begin(); iter != side_part_map.end(); iter++)
            {
              stk::mesh::PartVector& epv = iter->second;
              for (unsigned iepv=0; iepv < epv.size(); iepv++)
                {
                  std::cout << "Refiner::get_side_part_relations: side_part = " << std::setw(50) << iter->first->name() << " elem_part= " << std::setw(50) << epv[iepv]->name() << std::endl;
                }
            }
        }
    }

#undef EXTRA_PRINT_UR_FES

    void Refiner::removeFamilyTrees()
    {
      EXCEPTWATCH;

      doProgressPrint("Stage: Removing family_trees... ");

      elements_to_be_destroyed_type elements_to_be_destroyed(*m_eMesh.get_bulk_data());

      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( FAMILY_TREE_RANK );

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element   = bucket[iElement];
                  stk::mesh::Entity element_p = element;

                  elements_to_be_destroyed.insert(element_p);
                }
            }
        }
      removeElements(elements_to_be_destroyed);

      doProgressPrint("Stage: Removing family_trees, size() = ");
    }

    void Refiner::
    removeOldElements(unsigned irank, stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;
      if (CHECK_DEBUG) std::cout << "removeOldElements start ... " << std::endl;

      const stk::mesh::Part *oldPart = m_eMesh.getPart(breakPattern->getOldElementsPartName()+toString(rank));

      if (1 && oldPart)
        {
          const CellTopologyData * const cell_topo_data = m_eMesh.get_cell_topology(*oldPart);
          std::string ct_name = (cell_topo_data ? cell_topo_data->name : "");
          //std::cout << "tmp removeOldElements::name= " << oldPart->name() << " for rank= " << rank << " topology= " << ct_name << std::endl;
        }

      if (!oldPart)
        {
          std::cout << "name= " << breakPattern->getOldElementsPartName()+toString(rank) << std::endl;
          throw std::runtime_error("oldPart is null");
        }

      stk::mesh::Selector removePartSelector (*oldPart);

      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( rank );

      elements_to_be_destroyed_type elements_to_be_destroyed(*m_eMesh.get_bulk_data());

      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (removePartSelector(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();

              if (0)
                {
                  std::string str;
                  stk::mesh::PartVector const& pv = bucket.supersets();
                  for (unsigned ip = 0; ip < pv.size(); ip++)
                    {
                      str += " "+pv[ip]->name();
                    }
                  std::cout << "P[" << m_eMesh.get_rank() << "] removing elements in bucket of parts: " << str << std::endl;
                }

              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  stk::mesh::Entity element_p = element;

                  if (!m_eMesh.isGhostElement(element))
                    {
                      elements_to_be_destroyed.insert(element_p);
                      //std::cout << "tmp removing elem = "; m_eMesh.print(element_p);
                    }
                }
            }
        }

      //getNodeRegistry().clear_elements_to_be_deleted(&elements_to_be_destroyed);
      removeElements(elements_to_be_destroyed, irank);
      if (CHECK_DEBUG) std::cout << "removeOldElements  ... end " << std::endl;

    }

    void Refiner::removeElements(elements_to_be_destroyed_type& elements_to_be_destroyed, unsigned irank)
    {
      elements_to_be_destroyed_type elements_to_be_destroyed_pass2(*m_eMesh.get_bulk_data());

      int progress_meter_num_total = elements_to_be_destroyed.size();
      if (m_doProgress && m_eMesh.get_rank() == 0 && progress_meter_num_total)
        {
          std::ostringstream oss; oss << "Delete Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::INIT, 0.0, oss.str());
          notifyObservers(&pd);
        }
      unsigned num_elem = 0;
      int progress_meter_when_to_post = progress_meter_num_total / m_progress_meter_frequency;
      if (0 == progress_meter_when_to_post)
        progress_meter_when_to_post = 1;
      double d_progress_meter_num_total = progress_meter_num_total;

      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed.begin(); itbd != elements_to_be_destroyed.end();  ++itbd)
        {
          stk::mesh::Entity element_p = *itbd;

          if (m_doProgress && m_eMesh.get_rank() == 0 && (num_elem % progress_meter_when_to_post == 0) )
            {
              double progress_meter_percent = 100.0*((double)num_elem)/std::max(d_progress_meter_num_total,1.0);
              std::ostringstream oss; oss << "Delete Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
              ProgressMeterData pd(ProgressMeterData::RUNNING, progress_meter_percent, oss.str());
              notifyObservers(&pd);
            }

          ++num_elem;

          if (0)
            {
              std::cout << "tmp removeElements removing element_p = " << element_p << std::endl;
              if (m_eMesh.is_valid(element_p)) std::cout << "tmp removeElements removing id= " << m_eMesh.identifier(element_p) << std::endl;
            }

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( element_p ) )
            {
              elements_to_be_destroyed_pass2.insert(element_p);
              //throw std::logic_error("Refiner::removeElements couldn't remove element");

            }
        }

      if (m_doProgress && m_eMesh.get_rank() == 0 && progress_meter_num_total)
        {
          std::ostringstream oss; oss << "Delete Elements " << " pass [" << 100.0*((double)irank)/((double)m_ranks.size()) << " %]" << " cpu: " << m_eMesh.cpu_time() << " [sec]";
          ProgressMeterData pd(ProgressMeterData::FINI, 0.0, oss.str());
          notifyObservers(&pd);
        }

      //std::cout << "tmp Refiner::removeElements pass2 size = " << elements_to_be_destroyed_pass2.size() << std::endl;
      for (elements_to_be_destroyed_type::iterator itbd = elements_to_be_destroyed_pass2.begin();
           itbd != elements_to_be_destroyed_pass2.end();  ++itbd)
        {
          stk::mesh::Entity element_p = *itbd;

          if (1)
            {
              const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
              stk::mesh::EntityRank rank= m_eMesh.entity_rank(element_p);
              for (stk::mesh::EntityRank higher_rank = static_cast<stk::mesh::EntityRank>(rank+1); higher_rank <= FAMILY_TREE_RANK; ++higher_rank)
                {
                  while (true)
                    {
                      percept::MyPairIterRelation rels (m_eMesh, element_p, higher_rank);
                      if (!rels.size())
                        break;
                      stk::mesh::Entity to_rel = rels[0].entity();
                      stk::mesh::RelationIdentifier to_id = rels[0].relation_ordinal();

                      bool del = m_eMesh.get_bulk_data()->destroy_relation( to_rel, element_p, to_id);
                      if (!del)
                        throw std::runtime_error("removeOldElements:: destroy_relation failed 1");
                    }
                }
            }

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( element_p ) )
            {
              shards::CellTopology cell_topo(m_eMesh.get_cell_topology(element_p));
              std::cout << "tmp Refiner::removeElements couldn't remove element in pass2,...\n tmp destroy_entity returned false: cell= " << cell_topo.getName() << std::endl;
              const percept::MyPairIterRelation elem_relations (m_eMesh, element_p, static_cast<stk::mesh::EntityRank>(m_eMesh.entity_rank(element_p)+1));
              std::cout << "tmp elem_relations[rank+1].size() = " << elem_relations.size() << std::endl;
              const percept::MyPairIterRelation elem_relations_elem (m_eMesh, element_p, m_eMesh.element_rank());
              std::cout << "tmp elem_relations[3].size() = " << elem_relations_elem.size() << std::endl;

              throw std::logic_error("Refiner::removeElements couldn't remove element, destroy_entity returned false.");
            }
        }
    }

    /// fix names of surfaces (changing for example surface_hex8_quad4 to surface_tet4_tri3)
    void Refiner::
    fixSurfaceAndEdgeSetNames(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;
      stk::mesh::PartVector toParts = breakPattern->getToParts();

      //std::cout << "toParts.size()= " << toParts.size() << " typeid= " << typeid(*breakPattern).name()  << std::endl;

      for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
        {
          //const std::string & partName = toParts[i_part]->name();
          std::string * toPartName_p = const_cast<std::string *> (&toParts[i_part]->name());

          std::string toPartName = toParts[i_part]->name();
          if ( toPartName.find("surface_", 0) == std::string::npos)
            {
              if (0) std::cout << "tmp fixSurfaceAndEdgeSetNames:: skipping toPartName= " << toPartName << " typeid= " << typeid(*breakPattern).name()  << std::endl;
              continue;
            }

          std::string newToPartName = toPartName;

          StringStringMap::iterator map_it;
          StringStringMap str_map =  breakPattern->fixSurfaceAndEdgeSetNamesMap();
          if (0) std::cout << "tmp fixSurfaceAndEdgeSetNamesMap:: str_map.size()= " << str_map.size()
            //<< " " << breakPattern->getFromTopoPartName() << "__" << breakPattern->getToTopoPartName()
                           << " typeid= " << typeid(*breakPattern).name()
                           << std::endl;

          for (map_it = str_map.begin(); map_it != str_map.end(); map_it++)
            {
              std::string from_str = map_it->first;
              std::string to_str = map_it->second;
              Util::replace(newToPartName, from_str, to_str);
              if (0)
                std::cout << "tmp fixSurfaceAndEdgeSetNamesMap: old= " << toPartName << " new= " << newToPartName << std::endl;
            }

          *toPartName_p = newToPartName;

          if (0)
            std::cout << "tmp fixSurfaceAndEdgeSetNamesMap:: P[" << m_eMesh.get_rank() << "] new part name= " << toParts[i_part]->name()
                      << " old part name = " << toPartName
                      << std::endl;
        }
    }

    // FIXME this is a hack to rename parts
    /// Renames as follows:
    ///   originalPartName -> originalPartName_uo_1000    The original part holds the elements to be converted, and is renamed to be the "old" part
    ///   originalPartName_urpconv -> originalPartName    The new part has the same name as the original part with urpconv appended, which
    ///                                                      is then changed back to the original part name
    ///   fromPartName -> fromPartName+"_uo_1000"
    ///   toPartName_urpconv -> toPartName
    ///
    /// So, after the renaming, the original part name holds the new elements, and the original elements are
    ///   in the part with the original name appended with _uo_1000.  These parts are ignored on subsequent input.
    ///
#define DEBUG_RENAME_NEW_PARTS 0
#define DEBUG_RENAME_NEW_PARTS_1 0
    static std::string strip_hashes(std::string in, stk::mesh::PartVector& parts, const std::string& sep, bool check=true)
    {
      std::string out=in;
      size_t pos=0;
      pos = in.find(sep);
      if (pos != std::string::npos)
        {
          std::string o1 = in.substr(pos+1);
          size_t pos2 = o1.find(sep);
          out = in.substr(0,pos)+o1.substr(pos2+1);
        }

      if (DEBUG_RENAME_NEW_PARTS_1) std::cout << "tmp srk in= " << std::setw(50) << in << " out= " << std::setw(50) << out << std::endl;
      if (check && out != in)
        {
          for (unsigned ii=0; ii < parts.size(); ++ii)
            {
              if (parts[ii]->name() == out)
                {
                  std::cout << "bad part[" << ii << "] = " << parts[ii]->name() << std::endl;
                  for (unsigned jj=0; jj < parts.size(); ++jj)
                    {
                      std::cout << "part[" << jj << "] = " << parts[jj]->name() << std::endl;
                    }
                  throw std::runtime_error("bad name change");
                }
            }
        }
      return out;
    }

    void Refiner::
    renameNewParts(stk::mesh::EntityRank rank, UniformRefinerPatternBase* breakPattern)
    {
      EXCEPTWATCH;

      stk::mesh::PartVector toParts = breakPattern->getToParts();
      stk::mesh::PartVector fromParts = breakPattern->getFromParts();
      bool do_strip_hashes = breakPattern->m_do_strip_hashes;
      bool do_strip_hashes_from = false;

      stk::mesh::PartVector all_parts = m_eMesh.get_fem_meta_data()->get_parts();

      if (DEBUG_RENAME_NEW_PARTS)
        {
          std::cout << "\n\ntmp srk getFromTopoPartName()= " << breakPattern->getFromTopoPartName() << " getToTopoPartName()= " << breakPattern->getToTopoPartName() << std::endl;
          std::cout << "fromParts.size() = " << fromParts.size() << " toParts.size() = " << toParts.size() << std::endl;
          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              std::cout << "tmp toParts[i_part]->name() = " << toParts[i_part]->name() << std::endl;
            }
          for (unsigned i_part = 0; i_part < fromParts.size(); i_part++)
            {
              std::cout << " fromParts[i_part]->name() = " << fromParts[i_part]->name()  << std::endl;
            }
        }

      VERIFY_OP_ON(fromParts.size(), ==, toParts.size(), "sizes wrong");

      if (fromParts.size() == toParts.size())
        {
          for (unsigned i_part = 0; i_part < toParts.size(); i_part++)
            {
              if (DEBUG_RENAME_NEW_PARTS) std::cout << "tmp before: fromPartName= " << fromParts[i_part]->name()
                                                    << " toPartName= " << toParts[i_part]->name() << std::endl;

              std::string toPartName = toParts[i_part]->name();
              if (toParts[i_part]->name() == fromParts[i_part]->name())
                {
                  continue;
                }
              std::string fromPartName = toPartName;
              int len = fromPartName.length();
              int clen = breakPattern->getAppendConvertString().length();
              fromPartName.erase(len - clen, clen);
              //mesh::Part *fromPart = m_eMesh.get_non_const_part(fromPartName);
              stk::mesh::Part *fromPart = fromParts[i_part];
              VERIFY_OP_ON(fromPart, !=, 0, std::string("Refiner::renameNewParts null fromPart found, fromPart= ")+fromPartName);

              if (1)
                {
                  std::string newToPartName = fromPartName;
                  std::string newFromPartName = fromPartName + breakPattern->getAppendOriginalString();
                  if (do_strip_hashes) newToPartName = strip_hashes(newToPartName, all_parts, breakPattern->getConvertSeparatorString(), false);
                  if (do_strip_hashes_from) newFromPartName = strip_hashes(newFromPartName, all_parts, breakPattern->getConvertSeparatorString(), false);
                  stk::io::set_alternate_part_name(*toParts[i_part], newToPartName);
                  stk::io::set_alternate_part_name(*fromParts[i_part], newFromPartName);
                }

              if (DEBUG_RENAME_NEW_PARTS) {
                std::cout << "tmp  after: fromPartName= " << fromParts[i_part]->name() << " toPartName= " << toParts[i_part]->name() << std::endl;
                std::cout << "tmp P[" << m_eMesh.get_rank() << "] fromPartName: " << fromPartName << " part= " << toParts[i_part]->name()
                          << " old part name = " << fromPart->name()
                          << std::endl;
              }
            }
        }
    }


    RefinementInfo&
    Refiner::
    getRefinementInfo()
    {
      return m_refinementInfo;
    }

    void
    Refiner::
    setQueryPassOnly(bool doQueryOnly)
    {
      m_doQueryOnly = doQueryOnly;
    }


    void Refiner::
    add_children_to_parts()
    {
      static std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part*>(0));
      static std::vector<stk::mesh::Part*> remove_parts;

      std::vector<stk::mesh::Entity> children;
      const std::vector< stk::mesh::Part * > & parts = m_eMesh.get_fem_meta_data()->get_parts();
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

      std::vector<std::pair<stk::mesh::Entity,stk::mesh::Part*> > toChange;

      unsigned nparts = parts.size();
      for (unsigned ipart=0; ipart < nparts; ipart++)
        {
          stk::mesh::Part& part = *parts[ipart];

          if (stk::mesh::is_auto_declared_part(part))
            continue;

          bool percept_auto_part = part.attribute<percept::AutoPart>() != 0;
          if (percept_auto_part)
            continue;

          const stk::mesh::EntityRank part_rank = part.primary_entity_rank();

          if (part_rank == stk::topology::NODE_RANK || part_rank == stk::topology::INVALID_RANK)
            continue;

          stk::mesh::Selector selector(part);
          //std::cout << "part= " << part.name() << std::endl;
          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( part_rank );

          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              if (selector(bucket))
                {
                  const unsigned num_entity_in_bucket = bucket.size();
                  for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                    {
                      stk::mesh::Entity element = bucket[ientity];
                      if (m_eMesh.hasFamilyTree(element))
                        {
                          m_eMesh.getChildren(element, children, false);
                          for (unsigned ich = 0; ich < children.size(); ich++)
                            {
                              if (on_locally_owned_part(m_eMesh.bucket(children[ich])) && !selector(m_eMesh.bucket(children[ich]))
                                  && m_eMesh.bucket(children[ich]).topology() == m_eMesh.bucket(element).topology())
                                {
                                  toChange.push_back(std::make_pair(children[ich],&part));
                                }
                            }
                        }
                    }
                }
            }
        }
      for (unsigned ii=0; ii < toChange.size(); ++ii)
        {
          stk::mesh::Entity element = toChange[ii].first;
          add_parts[0] = toChange[ii].second;
          m_eMesh.get_bulk_data()->change_entity_parts( element, add_parts, remove_parts );
        }
    }

#if 1
    void Refiner::set_active_part()
    {
      set_active_part(m_eMesh);
    }

    // static
    void Refiner::set_active_part(PerceptMesh& eMesh)
    {
      // deal with parts
      stk::mesh::EntityRank part_ranks[] = {eMesh.element_rank(), eMesh.side_rank()};
      for (unsigned irank=0; irank < 2; irank++)
        {
          std::string active_part_name = "refine_active_elements_part_"+toString(part_ranks[irank]);
          std::string inactive_part_name = "refine_inactive_elements_part_"+toString(part_ranks[irank]);
          stk::mesh::Part* child_elements_part = eMesh.get_non_const_part(active_part_name);
          stk::mesh::Part* parent_elements_part = eMesh.get_non_const_part(inactive_part_name);
          stk::mesh::Selector in_child_part(*child_elements_part);
          stk::mesh::Selector in_parent_part(*parent_elements_part);

          if (child_elements_part && parent_elements_part)
            {
              std::vector<stk::mesh::Part*> child_parts(1, child_elements_part);
              std::vector<stk::mesh::Part*> parent_parts(1, parent_elements_part);
              stk::mesh::Selector on_locally_owned_part =  ( eMesh.get_fem_meta_data()->locally_owned_part() );

              std::vector<stk::mesh::Entity> child_entities;
              std::vector<stk::mesh::Entity> parent_entities;

              const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( part_ranks[irank] );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;

                  if (on_locally_owned_part(bucket))
                    {
                      bool in_child_part_element = in_child_part(bucket);
                      bool in_parent_part_element = in_parent_part(bucket);
                      const unsigned num_entity_in_bucket = bucket.size();
                      for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                        {
                          stk::mesh::Entity element = bucket[ientity];
                          if (eMesh.hasFamilyTree(element) && eMesh.isParentElement(element, true))
                            {
                              if (in_child_part_element || !in_parent_part_element)
                              {
                                parent_entities.push_back(element);
                              }
                            }
                          else
                            {
                              if (!in_child_part_element || in_parent_part_element)
                              {
                                child_entities.push_back(element);
                              }
                            }
                        }
                    }
                }

              //std::cout << "tmp Refiner::set_active_part: child_entities= " << child_entities.size() << " parent_entities= " << parent_entities.size() << std::endl;
              for (unsigned iv=0; iv < child_entities.size(); iv++)
                {
                  eMesh.get_bulk_data()->change_entity_parts( child_entities[iv],   child_parts, parent_parts );
                }
              for (unsigned iv=0; iv < parent_entities.size(); iv++)
                {
                  eMesh.get_bulk_data()->change_entity_parts( parent_entities[iv],  parent_parts, child_parts );
                }
            }
        }

    }

#else
    void Refiner::set_active_part()
    {
      // deal with parts
      stk::mesh::EntityRank part_ranks[] = {m_eMesh.element_rank(), m_eMesh.side_rank()};
      for (unsigned irank=0; irank < 2; irank++)
        {
          std::string active_part_name = "refine_active_elements_part_"+toString(part_ranks[irank]);
          std::string inactive_part_name = "refine_inactive_elements_part_"+toString(part_ranks[irank]);
          stk::mesh::Part* child_elements_part = m_eMesh.get_non_const_part(active_part_name);
          stk::mesh::Part* parent_elements_part = m_eMesh.get_non_const_part(inactive_part_name);
          stk::mesh::Selector in_child_part(*child_elements_part);
          stk::mesh::Selector in_parent_part(*parent_elements_part);

          if (child_elements_part && parent_elements_part)
            {
              std::vector<stk::mesh::Part*> empty;
              std::vector<stk::mesh::Part*> child_parts(1, child_elements_part);
              std::vector<stk::mesh::Part*> parent_parts(1, parent_elements_part);
              stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );

              std::vector<stk::mesh::Entity> child_entities_add_to_child;
              std::vector<stk::mesh::Entity> child_entities_remove_from_parent;
              std::vector<stk::mesh::Entity> parent_entities_add_to_parent;
              std::vector<stk::mesh::Entity> parent_entities_remove_from_child;

              const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( part_ranks[irank] );

              for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
                {
                  stk::mesh::Bucket & bucket = **k ;

                  if (on_locally_owned_part(bucket))
                    {
                      bool in_child_part_element = in_child_part(bucket);
                      bool in_parent_part_element = in_parent_part(bucket);
                      const unsigned num_entity_in_bucket = bucket.size();
                      for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                        {
                          stk::mesh::Entity element = bucket[ientity];
                          if (m_eMesh.hasFamilyTree(element) && m_eMesh.isParentElement(element, true))
                            {
                              if (in_child_part_element)
                                {
                                  parent_entities_remove_from_child.push_back(element);
                                }
                              if (!in_parent_part_element)
                                {
                                  parent_entities_add_to_parent.push_back(element);
                                }
                            }
                          else
                            {
                              if (!in_child_part_element)
                              {
                                child_entities_add_to_child.push_back(element);
                              }
                              if (in_parent_part_element)
                              {
                                child_entities_remove_from_parent.push_back(element);
                              }
                            }
                        }
                    }
                }

              //std::cout << "tmp Refiner::set_active_part: child_entities= " << child_entities.size() << " parent_entities= " << parent_entities.size() << std::endl;
              for (unsigned iv=0; iv < child_entities_add_to_child.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( child_entities_add_to_child[iv],   child_parts, empty);
                }
              for (unsigned iv=0; iv < parent_entities_add_to_parent.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( parent_entities_add_to_parent[iv],  parent_parts, empty);
                }
              for (unsigned iv=0; iv < child_entities_remove_from_parent.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( child_entities_remove_from_parent[iv],   empty, parent_parts );
                }
              for (unsigned iv=0; iv < parent_entities_remove_from_child.size(); iv++)
                {
                  m_eMesh.get_bulk_data()->change_entity_parts( parent_entities_remove_from_child[iv], empty, child_parts );
                }
            }
        }

    }
#endif
    //    ========================================================================================================================
    //    ========================================================================================================================
    //    ========================================================================================================================


    void Refiner::check_db(std::string msg)
    {
      std::cout << "P[" << m_eMesh.get_rank() << "] tmp check_db msg= " << msg << std::endl;
      check_db_entities_exist(msg);
      check_db_ownership_consistency(msg);
      std::cout << "P[" << m_eMesh.get_rank() << "] tmp check_db done msg= " << msg << std::endl;
    }


    void Refiner::check_db_entities_exist(std::string msg)
    {
      std::vector<stk::mesh::EntityRank> ranks_to_check;
      ranks_to_check.push_back(m_eMesh.node_rank());
      ranks_to_check.push_back(stk::topology::ELEMENT_RANK);
      for (unsigned irank=0; irank < ranks_to_check.size(); irank++)
        {

          const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( ranks_to_check[irank] );

          for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
            {
              stk::mesh::Bucket & bucket = **k ;

              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity element = bucket[ientity];
                  VERIFY_OP_ON(m_eMesh.is_valid(element), ==, true, "check_db_entities_exist bad element");
                  stk::mesh::Entity element_1 = m_eMesh.get_bulk_data()->get_entity(ranks_to_check[irank], m_eMesh.identifier(element));
                  if (element != element_1 || m_eMesh.identifier(element) != m_eMesh.identifier(element_1))
                    {
                      std::cout << "msg= " << msg << " error element, element_1, ids= "
                                << &element << " " << element_1 << " " << m_eMesh.identifier(element) << " " << m_eMesh.identifier(element_1) << std::endl;
                      throw std::logic_error("check_db_entities_exist:: error #1");
                    }

                }
            }
        }
    }

    void Refiner::check_db_ownership_consistency(std::string msg)
    {
      SubDimCellToDataMap& cell_2_data_map = m_nodeRegistry->getMap();
      stk::mesh::BulkData &bulk_data = *m_eMesh.get_bulk_data();

      for (SubDimCellToDataMap::iterator cell_iter = cell_2_data_map.begin(); cell_iter != cell_2_data_map.end(); ++cell_iter)
        {
          SubDimCellData& nodeId_elementOwnderId = (*cell_iter).second;
          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          stk::mesh::EntityRank owning_elementRank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

          if (nodeIds_onSE.size())
            {

              if (!owning_elementId)
                throw std::logic_error("check_db_ownership_consistency:: error #1 msg= "+msg);

              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);

              if (!m_eMesh.is_valid(owning_element))
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] "
                            << " error check_db_ownership_consistency: msg= " << msg << " owning_elementId= " << owning_elementId << std::endl;
                  throw std::logic_error("check_db_ownership_consistency:: error #2, msg= "+msg);
                }

              if (stk::mesh::Deleted == bulk_data.state(owning_element) )
                {
                  std::cout << "P[" << m_eMesh.get_rank() << "] "
                            << " error check_db_ownership_consistency: Deleted, msg= " << msg << " owning_elementId= " << owning_elementId << std::endl;
                  throw std::logic_error("check_db_ownership_consistency:: error #2.0, msg= "+msg);
                }


              if (m_eMesh.identifier(owning_element) != owning_elementId)
                {
                  std::cout << "msg= " << msg << " check_db_ownership_consistency error element, element_1, ids= "
                            << owning_elementId << " " << m_eMesh.identifier(owning_element) << std::endl;
                  throw std::logic_error("check_db_ownership_consistency:: error #1.1, msg= "+msg);
                }

              if (!m_eMesh.isGhostElement(owning_element))
                {

                  for (unsigned inode = 0; inode < nodeIds_onSE.size(); inode++)
                    {
                      stk::mesh::Entity node = nodeIds_onSE[inode];
                      if (!m_eMesh.is_valid(node))
                        throw std::logic_error("check_db_ownership_consistency:: error #3, msg= "+msg);

                      stk::mesh::Entity node1 = m_eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, nodeIds_onSE.m_entity_id_vector[inode]);
                      if (!m_eMesh.is_valid(node1))
                        throw std::logic_error("check_db_ownership_consistency:: error #3a, msg= "+msg);

                      stk::mesh::Entity node2 = m_eMesh.get_bulk_data()->get_entity(stk::topology::NODE_RANK, m_eMesh.identifier(node) );
                      if (!m_eMesh.is_valid(node2))
                        throw std::logic_error("check_db_ownership_consistency:: error #3b, msg= "+msg);
                      if (node != node2)
                        throw std::logic_error("check_db_ownership_consistency:: error #3c, msg= "+msg);

                    }
                }
            }
        }
    }

    void Refiner::initializeDB(bool use_rebuild_node_registry)
    {
      if (m_eMesh.getProperty("use_rebuild_node_registry") == "true")
        use_rebuild_node_registry = true;
      if (m_eMesh.getProperty("use_rebuild_node_registry") == "false")
        use_rebuild_node_registry = false;
      if (use_rebuild_node_registry)
        {
          RefinerUtil::rebuild_node_registry(m_eMesh, *m_nodeRegistry, true);
          return;
        }
      this->special_processing("refiner_pre_initializeDB");
      bool debug=false;
      m_nodeRegistry->initialize();
      //m_nodeRegistry->getMap().clear();
      m_nodeRegistry->init_comm_all();

      //std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk INFO: ranks.size()= " << m_ranks.size() << " m_breakPattern.size()= " << m_breakPattern.size() << std::endl;
      for (unsigned irank = 0; irank < m_ranks.size(); irank++)
        {
          unsigned elementType = m_breakPattern[irank]->getFromTypeKey();
          {
            EXCEPTWATCH;

            ElementRankTypeInfo& e_info = m_elementRankTypeInfo[irank];
            VERIFY_OP_ON(m_ranks[irank], ==, e_info.first,"er1");
            VERIFY_OP_ON(elementType, ==, e_info.second,"er2");

            vector<NeededEntityType> needed_entity_ranks;
            m_breakPattern[irank]->fillNeededEntities(needed_entity_ranks);

            bool count_only = false;
            bool doAllElements = true;

            unsigned num_elem_not_ghost_0_incr = doForAllElements(irank, "initializeEmpty",
                                                                  m_ranks[irank], &NodeRegistry::initializeEmpty,
                                                                  elementType, needed_entity_ranks,
                                                                  count_only, doAllElements);

            if (debug) std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk1 irank= " << irank << " ranks[irank]= " << m_ranks[irank]
                             << " nodeRegistry size= " << m_nodeRegistry->getMap().size()
                             << " num_elem_not_ghost_0_incr= " << num_elem_not_ghost_0_incr
                             << std::endl;
          }
        }

      if (use_rebuild_node_registry)
        {
          RefinerUtil::rebuild_node_registry(m_eMesh, *m_nodeRegistry, false);
        }
      else {

      SubDimCellToDataMap& cell_2_data_map = m_nodeRegistry->getMap();
      if (debug) std::cout << "tmp srk1 cell_2_data_map size= " << cell_2_data_map.size() << std::endl;

      for (SubDimCellToDataMap::iterator cell_iter = cell_2_data_map.begin(); cell_iter != cell_2_data_map.end(); ++cell_iter)
        {
          const SubDimCell_SDCEntityType& subDimEntity = (*cell_iter).first;
          SubDimCellData& nodeId_elementOwnderId = (*cell_iter).second;
          stk::mesh::EntityId owning_elementId = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().id();
          NodeIdsOnSubDimEntityType& nodeIds_onSE = nodeId_elementOwnderId.get<SDC_DATA_GLOBAL_NODE_IDS>();
          stk::mesh::EntityRank owning_elementRank = nodeId_elementOwnderId.get<SDC_DATA_OWNING_ELEMENT_KEY>().rank();

          if (0 && debug) std::cout << "P[" << m_eMesh.get_rank() << "] tmp srk1 subDimEntity.size= " << subDimEntity.size() << " owning_elementId= " << owning_elementId << " nodeIds_onSE.size= " << nodeIds_onSE.size() << std::endl;
          // if (subDimEntity.size() == 1)
          //   continue;

          if (1)
            {
              stk::mesh::Entity owning_element = m_eMesh.get_bulk_data()->get_entity(owning_elementRank, owning_elementId);
              VERIFY_OP_ON(m_eMesh.is_valid(owning_element), ==, true, "hmmm");
              //double ela = m_eMesh.edge_length_ave(owning_element);
              double centroid[3] = {0,0,0};
              int nDim = m_eMesh.get_spatial_dim();
              if (subDimEntity.size() == 1)
                {
                  computeCentroid(subDimEntity[0], centroid, *(m_eMesh.get_coordinates_field()));
                }
              else
                {
                  for (unsigned ii=0; ii < subDimEntity.size(); ++ii)
                    {
                      stk::mesh::Entity node = subDimEntity[ii];
                      double *coords = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), node);
                      for (int jj=0; jj < nDim; ++jj)
                        {
                          centroid[jj] += coords[jj]/double(subDimEntity.size());
                        }
                    }
                }

              std::vector<stk::mesh::Entity> children;
              double min_dist = std::numeric_limits<double>::max();
              stk::mesh::Entity min_node = stk::mesh::Entity();
              const percept::MyPairIterRelation parent_nodes(m_eMesh, owning_element, m_eMesh.node_rank());

              m_eMesh.getChildren(owning_element, children, false);
              if (children.size() == 0)
                continue;

              double emin = std::numeric_limits<double>::max();
              for (unsigned ich = 0; ich < children.size(); ich++)
                {
                  double min_edge_len = 0, max_edge_len = 0;
                  m_eMesh.edge_length_ave(owning_element, m_eMesh.get_coordinates_field(), &min_edge_len, &max_edge_len);
                  if (min_edge_len < emin)
                    emin = min_edge_len;

                  const percept::MyPairIterRelation child_nodes(m_eMesh, children[ich], m_eMesh.node_rank());
                  for (unsigned jj=0; jj < child_nodes.size(); ++jj)
                    {
                      stk::mesh::Entity cnode = child_nodes[jj].entity();
                      bool fnd=false;
                      for (unsigned kk=0; kk < parent_nodes.size(); ++kk)
                        {
                          if (cnode == parent_nodes[kk].entity())
                            {
                              fnd=true;
                              break;
                            }
                        }
                      if (fnd) continue;
                      double *c_coords = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), cnode);
                      double dist=0;
                      for (int k=0; k < nDim; ++k)
                        {
                          dist += (centroid[k] - c_coords[k])*(centroid[k] - c_coords[k]);
                        }
                      dist = std::sqrt(dist);
                      if (dist < min_dist)
                        {
                          min_dist = dist;
                          min_node = cnode;
                        }
                    }
                }
              bool found = min_dist < 1.e-3*emin;

              // FIXME - what about quad faces?
              if (0 && !found)
                {
                  double ela_local = 0.0;
                  double nela = 0.0;
                  for (unsigned ii=0; ii < subDimEntity.size()-1; ++ii)
                    {
                      double *coords_ii = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), subDimEntity[ii]);
                      for (unsigned jj=ii+1; jj < subDimEntity.size(); ++jj)
                        {
                          double *coords_jj = m_eMesh.field_data_inlined(m_eMesh.get_coordinates_field(), subDimEntity[jj]);
                          nela += 1.0;
                          double dist=0.0;
                          for (int k=0; k < nDim; ++k)
                            {
                              dist += (coords_ii[k] - coords_jj[k])*(coords_ii[k] - coords_jj[k]);
                            }
                          dist = std::sqrt(dist);
                          ela_local += dist;
                        }
                    }
                  ela_local /= nela;
                  found = min_dist < 1.e-4*ela_local;
                }
              if (found)
                {
                  if (0)
                    {
                      std::cout << "P[" << m_eMesh.get_rank() << "] found node= " << m_eMesh.identifier(min_node)
                                << " " << m_eMesh.print_entity_compact(min_node)
                                << std::endl;
                    }

                  nodeIds_onSE.resize(1);
                  nodeIds_onSE.m_entity_id_vector.resize(1);
                  nodeIds_onSE[0] = min_node;
                  nodeIds_onSE.m_entity_id_vector[0] = m_eMesh.identifier(min_node);
                  nodeIds_onSE.m_mark = NodeRegistry::NR_MARK;
                }
              else
                {
                  nodeIds_onSE.resize(0);
                  nodeIds_onSE.m_entity_id_vector.resize(0);
                }
            }
        }
      }

      if (1)
        {
          mod_begin();
          this->removeDanglingNodes();
          mod_end(0,"RefinerInitDB");
        }

      this->special_processing("refiner_post_initializeDB");

      //m_nodeRegistry->checkDB("after initializeDB");
    }

    void Refiner::removeFromNewNodesPart()
    {
      stk::mesh::Part* new_nodes_part = m_eMesh.get_non_const_part("refine_new_nodes_part");
      VERIFY_OP_ON(new_nodes_part, !=, 0, "new_nodes_part");

      std::vector<stk::mesh::Part*> remove_parts(1, new_nodes_part);
      std::vector<stk::mesh::Part*> add_parts;
      std::vector<stk::mesh::Entity> node_vec;

      stk::mesh::Selector removePartSelector(*new_nodes_part & m_eMesh.get_fem_meta_data()->locally_owned_part() );
      const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          stk::mesh::Bucket & bucket = **k ;
          if (removePartSelector(bucket))
            {
              const unsigned num_entity_in_bucket = bucket.size();
              for (unsigned ientity = 0; ientity < num_entity_in_bucket; ientity++)
                {
                  stk::mesh::Entity node = bucket[ientity];
                  node_vec.push_back(node);
                  if (m_eMesh.m_new_nodes_field)
                    {
                      NewNodesType_type *ndata = stk::mesh::field_data(*m_eMesh.m_new_nodes_field, node);
                      if (ndata)
                        {
                          ndata[0] = static_cast<NewNodesType_type>(0);
                        }
                    }
                }
            }
        }
      for (unsigned ii=0; ii < node_vec.size(); ii++)
        {
          m_eMesh.get_bulk_data()->change_entity_parts( node_vec[ii], add_parts, remove_parts );
        }
    }

  void Refiner::
  fix_side_sets_2(bool allow_not_found, SetOfEntities *avoid_elems, SetOfEntities *avoid_sides, RefinerSelector *sel, const std::string& msg)
  {
    if (m_avoidFixSideSets)
      return;

    FixSideSets fss(this, m_eMesh, m_excludeParts, m_side_part_map, m_geomFile, m_avoidFixSideSetChecks, sel, m_doProgress);
    fss.fix_side_sets_2(allow_not_found, avoid_elems, avoid_sides, msg);
  }

  void Refiner::
  build_side_set(SetOfEntities& side_set, bool only_roots)
  {
    FixSideSets fss(this, m_eMesh, m_excludeParts, m_side_part_map, m_geomFile, m_avoidFixSideSetChecks, 0, m_doProgress); // FIXME - RefinerSelector?
    fss.build_side_set(side_set, only_roots);
  }

  bool Refiner::bucket_acceptable(stk::mesh::Bucket& bucket, stk::mesh::EntityRank rank)
  {
    stk::mesh::PartVector const& side_parts = bucket.supersets();
    for (unsigned isp=0; isp < side_parts.size(); ++isp)
      {
        stk::mesh::Part& part = *side_parts[isp];
        bool is_auto = stk::mesh::is_auto_declared_part(part);
        const AutoPart *side_auto_part = part.attribute<AutoPart>();
        bool is_percept_auto_part = side_auto_part != 0;
        if (!is_percept_auto_part && !is_auto && part.primary_entity_rank() == rank)
          {
            return true;
          }
      }
    return false;
  }

  void Refiner::mod_begin(stk::diag::Timer *timer)
  {
    stk::diag::Timer& timerRoot_ = (getModBegEndRootTimer() ? *getModBegEndRootTimer() : (timer ? *timer :  rootTimer()));
    stk::diag::Timer timerModBeg_("percept::ModBeg", timerRoot_);
    stk::diag::TimeBlock timerModBegBlock_(timerModBeg_);

    m_eMesh.get_bulk_data()->modification_begin();

  }

  void Refiner::mod_end(stk::diag::Timer *timer, const std::string& msg)
  {
    stk::diag::Timer& timerRoot_ = (getModBegEndRootTimer() ? *getModBegEndRootTimer() : (timer ? *timer :  rootTimer()));
    stk::diag::Timer timerModEndAll_("ModEndAll", timerRoot_);
    stk::diag::Timer timerModEnd_("ModEnd"+msg, timerModEndAll_);
    stk::diag::TimeBlock timerModEndBlock_(timerModEnd_);
    stk::diag::TimeBlock timerModEndBlockAll_(timerModEndAll_);

    //stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
    m_eMesh.get_bulk_data()->modification_end();
  }

  void add_ft_nodes(PerceptMesh& eMesh, stk::mesh::Entity parent_elem, stk::mesh::Entity child_elem)
  {
    const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(eMesh.element_rank() + 1u);
    SetOfEntities to_add(*eMesh.get_bulk_data());

    MyPairIterRelation parent_ft(eMesh, parent_elem, FAMILY_TREE_RANK);
    unsigned parent_elem_ft_level_1 = 0;
    if (parent_ft.size() == 2)
      parent_elem_ft_level_1 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_1, parent_elem);

    stk::mesh::Entity family_tree = parent_ft[parent_elem_ft_level_1].entity();

    bool checkInShared = true;

    percept::MyPairIterRelation parent_elem_nodes (eMesh, parent_elem,  stk::topology::NODE_RANK );
    for (unsigned i = 0; i < parent_elem_nodes.size(); i++)
      {
        if (checkInShared && !eMesh.get_bulk_data()->in_shared(eMesh.key(parent_elem_nodes[i].entity()))) continue;

        bool found = false;
        percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
        for (unsigned j = 0; j < ft_nodes.size(); j++)
          {
            if (ft_nodes[j].entity() == parent_elem_nodes[i].entity())
              {
                found = true;
                break;
              }
          }
        if (!found)
          {
            to_add.insert(parent_elem_nodes[i].entity());
            VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(parent_elem_nodes[i].entity()), ==, true, "parent_elem_nodes bad");
            //eMesh.get_bulk_data()->declare_relation(family_tree, parent_elem_nodes[i].entity(), ft_nodes.size());
          }
      }

    percept::MyPairIterRelation child_elem_nodes (eMesh, child_elem,  stk::topology::NODE_RANK );
    if (child_elem_nodes.size() == 0)
      {
        throw std::runtime_error("child_elem has no nodes");
      }
    for (unsigned i = 0; i < child_elem_nodes.size(); i++)
      {
        if (checkInShared && !eMesh.get_bulk_data()->in_shared(eMesh.key(child_elem_nodes[i].entity()))) continue;

        bool found = false;
        percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
        for (unsigned j = 0; j < ft_nodes.size(); j++)
          {
            if (ft_nodes[j].entity() == child_elem_nodes[i].entity())
              {
                found = true;
                break;
              }
          }
        if (!found)
          {
            to_add.insert(child_elem_nodes[i].entity());
            VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(child_elem_nodes[i].entity()), ==, true, "child_elem_nodes bad");
            //eMesh.get_bulk_data()->declare_relation(family_tree, child_elem_nodes[i].entity(), ft_nodes.size());

          }
      }

    // check for second level and subsequent refinement
    if (parent_ft.size() == 2)
      {
        unsigned parent_elem_ft_level_0 = eMesh.getFamilyTreeRelationIndex(FAMILY_TREE_LEVEL_0, parent_elem);
        stk::mesh::Entity family_tree_level_0 = parent_ft[parent_elem_ft_level_0].entity();

        percept::MyPairIterRelation ft_level_0_nodes (eMesh, family_tree_level_0,  stk::topology::NODE_RANK );
        for (unsigned i = 0; i < ft_level_0_nodes.size(); i++)
          {
            if (checkInShared && !eMesh.get_bulk_data()->in_shared(eMesh.key(ft_level_0_nodes[i].entity()))) continue;

            bool found = false;
            percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
            for (unsigned j = 0; j < ft_nodes.size(); j++)
              {
                if (ft_nodes[j].entity() == ft_level_0_nodes[i].entity())
                  {
                    found = true;
                    break;
                  }
              }
            if (!found)
              {
                VERIFY_OP_ON(eMesh.get_bulk_data()->in_index_range(ft_level_0_nodes[i].entity()), ==, true, "ft_level_0_nodes bad 0");
                //eMesh.get_bulk_data()->declare_relation(family_tree, ft_level_0_nodes[i].entity(), ft_nodes.size());
                to_add.insert(ft_level_0_nodes[i].entity());
              }
          }
      }

    // add nodes to family_tree
    {
      percept::MyPairIterRelation ft_nodes (eMesh, family_tree,  stk::topology::NODE_RANK );
      unsigned ftns=ft_nodes.size();

      std::vector<stk::mesh::Entity> to_add_vec(to_add.begin(), to_add.end());

      for (unsigned ita=0; ita < to_add_vec.size(); ita++)
        {
          eMesh.get_bulk_data()->declare_relation(family_tree, to_add_vec[ita], ftns+ita);
        }
    }

  }

  void add_ft_nodes(PerceptMesh& eMesh, stk::mesh::Entity elem)
  {
    std::vector<stk::mesh::Entity> children;
    if (eMesh.hasFamilyTree(elem))
      {
        eMesh.getChildren(elem, children, true, false);
        if (children.size() == 0)
          {
            return;
          }
      }

    for (auto child : children)
      {
        add_ft_nodes(eMesh, elem, child);
      }
  }

  void delete_ft_nodes(PerceptMesh& eMesh, stk::mesh::Entity ft)
  {
    std::vector<stk::mesh::Entity> nodes(eMesh.get_bulk_data()->begin_nodes(ft), eMesh.get_bulk_data()->end_nodes(ft));
    std::vector<stk::mesh::ConnectivityOrdinal> nords( eMesh.get_bulk_data()->begin_node_ordinals(ft),  eMesh.get_bulk_data()->end_node_ordinals(ft));

    for (unsigned jj=0; jj < nodes.size(); ++jj)
      {
        bool del = eMesh.get_bulk_data()->destroy_relation( ft, nodes[jj], nords[jj]);

        VERIFY_OP_ON(del, ==, true, "bad del");
      }
  }

  void Refiner::reset_family_tree_to_node_relations()
  {

    const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(m_eMesh.element_rank() + 1u);

    std::vector<stk::mesh::Entity> ftvec;
    stk::mesh::get_selected_entities(m_eMesh.get_fem_meta_data()->locally_owned_part(), m_eMesh.get_bulk_data()->buckets(FAMILY_TREE_RANK), ftvec);
    for (auto ft : ftvec)
      {
        delete_ft_nodes(m_eMesh, ft);
      }

    for (stk::mesh::EntityRank rank_iter = m_eMesh.side_rank(); rank_iter <= m_eMesh.element_rank(); ++rank_iter)
      {
        SetOfEntities eset(*m_eMesh.get_bulk_data());
        std::vector<stk::mesh::Entity> evec;
        stk::mesh::get_selected_entities(m_eMesh.get_fem_meta_data()->locally_owned_part() , m_eMesh.get_bulk_data()->buckets(rank_iter), evec);
        for (auto elem : evec)
          {
            //eset.insert(elem);

            if (m_eMesh.numChildren(elem) == 0)
              {
                eset.insert(elem);

                stk::mesh::Entity parent = m_eMesh.getParent(elem, false);
                if (m_eMesh.is_valid(parent))
                  {
                    eset.insert(parent);
                    if (0)
                      {
                        stk::mesh::Entity grand_parent = m_eMesh.getParent(parent, false);
                        if (m_eMesh.is_valid(grand_parent))
                          {
                            eset.insert(grand_parent);
                          }
                      }
                  }
              }
          }
        for (auto elem : eset)
          {
            add_ft_nodes(m_eMesh, elem);
          }
      }
  }

} // namespace percept
