// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



#define IN_MESH_ADAPT 1
#define MESH_ADAPT_CPP11_MVI 0

#include <sys/unistd.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <stdint.h>
#include <map>
#include <list>
#include <string>

#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif
#if defined(WITH_KOKKOS)
#include <Kokkos_Core.hpp>
#endif

#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>
#include <percept/RunEnvironment.hpp>
#include <percept/ProgressMeter.hpp>
#include <percept/Histograms.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/memory_util.hpp>

#include <stk_mesh/base/MemoryUsage.hpp>

#include <adapt/RefinerUtil.hpp>

#include <adapt/main/RunAdaptRun.hpp>
#include <adapt/TransitionElementAdapter.hpp>
#include <adapt/ElementRefinePredicate.hpp>

#include <adapt/UniformRefiner.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <Ioss_Utils.h>
#include <Ioss_SerializeIO.h>

#include <adapt/SerializeNodeRegistry.hpp>
#include <adapt/AdaptedMeshVerifier.hpp>

#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <percept/mesh/geometry/stk_geom/3D/FitGregoryPatches.hpp>

#if defined( STK_PERCEPT_HAS_GEOMETRY )
#include <percept/DihedralAngleCheck.hpp>
#include <percept/mesh/geometry/kernel/GeometryKernelOpenNURBS.hpp>
#include <percept/mesh/geometry/stk_geom/LocalCubicSplineFit.hpp>
#endif
#include <percept/mesh/geometry/recovery/GeometryRecoverySplineFit.hpp>

#define ALLOW_MEM_TEST 1
#define DEBUG_ADAPT_MAIN 0

#include "AdaptMain.hpp"
#include "RunAdaptRun.hpp"
#include <stk_mesh/base/MeshUtils.hpp>
#include <adapt/main/MemoryMultipliers.hpp>

class PGeom;
#ifdef HAVE_ACIS
class PGeomACIS;
#endif

namespace percept {

  class MeshAdapt {

  public:

    MeshAdapt();

    double MegaByte(MemorySizeType x) { return  ((double)x/1024.0/1024.0); }

    enum GeometryType {
      GEOM_NONE,
      OPENNURBS,
      PGEOM_OPENNURBS,
      MESH_BASED_GREGORY_PATCH,
      PGEOM_ACIS,
      N_GeometryType
    };

    bool has_suffix(const std::string &str, const std::string &suffix) {
      return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    MemorySizeType memory_dump(int dump_level, const stk::ParallelMachine& comm, stk::mesh::BulkData& bulkData, NodeRegistry* node_reg, std::string msg);
    void test_memory(MemorySizeType n_elements, MemorySizeType n_nodes);
    // Determine the basename and refinement level from a given an ExodusII file name.
    // This is the last integer before the final '.e' or '.exo' in the file name.
    // If no level was found, returns 0 as the level.
    //
    // Return:
    //    string - the basename which is everything before the number
    //    int    - the refinement level
    std::pair<std::string,int> get_basename_and_level(const std::string& input_mesh);
    void checkInput(std::string option, std::string value, std::string allowed_values, Teuchos::CommandLineProcessor& clp);
    void print_simple_usage(int argc, char **argv);
    int check_for_simple_options(int argc, char **argv);

    int adapt_main_simple_options(int argc_in, char **argv_in);
    void dump_args(int argc, char **argv);
    void check_args(int argc, char **argv);

    struct EntitySelectorUCF {
      PerceptMesh& m_eMesh;
      int m_dim;
      EntitySelectorUCF(PerceptMesh& eMesh) : m_eMesh(eMesh), m_dim(m_eMesh.get_spatial_dim()) {}

      bool operator()(stk::mesh::Entity e_from, stk::mesh::Entity e_to)
      {
        VERIFY_OP_ON(m_eMesh.m_unprojected_coordinates, !=, 0, "bad m_unprojected_coordinates");
        double *c_from = stk::mesh::field_data( *m_eMesh.get_coordinates_field() , e_from );
        double *uc_to = stk::mesh::field_data( *m_eMesh.m_unprojected_coordinates , e_to );
        if (std::fabs(uc_to[3]) < 1.e-6)
          {
            for (int j=0; j < m_dim; j++)
              {
                uc_to[j] = c_from[j];
              }
            if (m_dim != 3)
              uc_to[2] = 0.0;
            uc_to[3] = 1.0;
          }
        // don't do further processing
        return false;
      }

    };

    // helpers (could be made virtual and specialized by derived,
    //  for now, we use if-stmts for special cases)
    void create_refine_pattern();
    void create_refiner();
    void pre_refiner_tasks(int iBreak);
    int  bucket_in_block_names(stk::mesh::Bucket& bucket, std::string& block);

    // utilities
    /// removes beams and shells in favor of node sets
    void do_convert_geometry_parts_OpenNURBS(std::string geometry_file, std::string newExodusFile);
    void setup_m2g_parts(std::string input_geometry);
    void initialize_m2g_geometry(std::string input_geometry);

    // sub-algorithms of do_run algorithms
    void do_precheck_memory_usage();
    int setup_options(Teuchos::CommandLineProcessor& clp, const stk::ParallelMachine& comm, int argc, char **argv);
    BlockNamesType process_block_names();
    void mesh_based_geometry_setup();
    void mesh_based_geometry_fitting();

    void verify_mesh_util(bool isInit);
    void verify_mesh_init();
    void verify_mesh_final();
    void do_histograms_init();
    void do_histograms_final();
    void compute_hmesh_sizes_init();
    void compute_hmesh_sizes_final();
    void do_test_memory();
    void run_adapt_run();
    void do_dihedral_angle_check();

    // do_run algorithms
    void pre_open();
    int do_run_pre_commit();
    int do_run_post_commit();
    int do_run_post_refine();
    int do_post_proc();

    // version - return true if built in Sierra
    bool get_version(std::string* v=0);

    void log_usage( bool status = true );

    // main routines
    int adapt_main(int argc, char **argv);
    int adapt_main_full_options(int argc, char **argv);
    int adapt_main_full_options_normal(int argc, char **argv);
    int adapt_main_full_options_streaming(int argc, char **argv);
    int main(int argc, char **argv);

    // member variable initialization (c++11) start
    // ====================================================================================================================================
#if MESH_ADAPT_CPP11_MVI
#define CPP11MVITYPE(x) x
#define CPP11MVI(x) {x}
#define CPP11MVI_3(x,y,z) = {x,y,z}
#define CPP11MVIEQ(x) x
#define CPP11MV1_ARRAY_INIT(type, var, sz, tmpvar, arrInit) type var[sz] arrInit;
#else
#define CPP11MVITYPE(x) x
#define CPP11MVI(x) = x
#define CPP11MVI_3(x,y,z)
#define CPP11MVIEQ(x) x
#define CPP11MV1_ARRAY_INIT(type, var, sz, tmpvar, arrInit) type var[sz];
#endif

#include "MeshAdaptMemberVarInit.hpp"

#undef CPP11MVITYPE
#undef CPP11MVI
#undef CPP11MVI_3
#undef CPP11MVIEQ
#undef CPP11MV1_ARRAY_INIT

    // ====================================================================================================================================
    // member variable initialization (c++11) end
    std::shared_ptr<FitGregoryPatches> fitter; // 3D
    std::shared_ptr<GeometryRecoverySplineFit> grsf; // 2D

    std::shared_ptr<UniformRefinerPatternBase> pattern;
    std::shared_ptr<Refiner> refiner;

    stk::mesh::FieldBase* proc_rank_field_ptr = 0;

    std::shared_ptr<PerceptMesh> eMeshP;
    std::shared_ptr<AdaptedMeshVerifier> adaptedMeshVerifier;
    PGeom * m_PGeomPntr = NULL;
#ifdef HAVE_ACIS
    PGeomACIS * m_PGA = NULL;
#endif

    BlockNamesType m_block_names;
    stk::mesh::Selector block_selector;

    std::shared_ptr<ElementRefinePredicate> element_refine_predicate;
    std::shared_ptr<stk::mesh::Selector> univ_selector;

    typedef std::map<std::string, int> StringIntMap;
    StringIntMap block_names_x_map;
    std::shared_ptr<DihedralAngleCheck> m_dihedral_angle_check;
  };


}

