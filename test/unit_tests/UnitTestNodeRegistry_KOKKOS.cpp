// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

#include <Teuchos_ScalarTraits.hpp>

#include <stk_io/IossBridge.hpp>

#include <percept/Util.hpp>
#include <percept/ExceptionWatch.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/function/ConstantFunction.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/fixtures/Fixture.hpp>

#include <adapt/NodeRegistry.hpp>
#include <adapt/NodeRegistry_KOKKOS.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/SerializeNodeRegistry.hpp>
#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <stk_util/parallel/Parallel.hpp>
#include <math.h>
#include <stk_mesh/base/MeshUtils.hpp>

namespace percept {
namespace unit_tests {

#define EXTRA_PRINT 0
//=============================================================================
//=============================================================================
//=============================================================================
TEST(nodeRegistry_KOKKOS, oldNR_v_newNR_create_new_nodes)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;

    const unsigned p_size = stk::parallel_machine_size( pm );
    if (p_size > 1) return;

/*  encore_rtest/postprocessors/integral_least_squares_fit_function/square_eight_tri3.e
    encore_rtest/postprocessors/recover_gradient_small_meshes/square_hole_tri3.e
    percept_rtest/mesh_adapt/heavy_tests/input_files._.square_tri3_uns.e
    percept_rtest/mesh_adapt/regression_tests/inputs/input_files._.square_tri3_0.e

    encore_rtest/postprocessors/recover_gradient_small_meshes/three_tet4.e
    encore_rtest/postprocessors/norms/norm_tet4/cubit_umr/cube_h0_tet4.e
    percept_rtest/mesh_adapt/regression_tests/inputs/input_files._.cube_tet4_uns.e*/


    percept::PerceptMesh eMesh_BOOST_NR(3u);
    percept::PerceptMesh eMesh_KOKKOS_NR(3u);

    eMesh_BOOST_NR.open("square_eight_tri3.e");
    eMesh_KOKKOS_NR.open("square_eight_tri3.e");

    {
      stk::mesh::Part& part = eMesh_KOKKOS_NR.get_fem_meta_data()->declare_part_with_topology("refine_new_nodes_part", stk::topology::NODE);
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      meta.declare_attribute_no_delete(part, &percept::auto_part);
    }

    {
      stk::mesh::Part& part = eMesh_BOOST_NR.get_fem_meta_data()->declare_part_with_topology("refine_new_nodes_part", stk::topology::NODE);
      stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
      meta.declare_attribute_no_delete(part, &percept::auto_part);
    }

    eMesh_KOKKOS_NR.commit();
    eMesh_BOOST_NR.commit();

    stk::mesh::EntityRank stk_mesh_Edge = stk::topology::EDGE_RANK;
    NeededEntityType needed_entity_rank( stk_mesh_Edge, 1u);
    std::vector<NeededEntityType> needed_entity_ranks(1, needed_entity_rank);

    //DO REGISTRATION FOR NEW NR
    std::cout << "About to do registration for kokkosNR\n";
    unsigned sizeHint = 3*8; //no edges per triangle X number of triangles = the max number of new nodes we would require
    NodeRegistry_KOKKOS kokkosNR(eMesh_KOKKOS_NR,0,true,sizeHint);
    kokkosNR.initialize();
    kokkosNR.beginRegistration();

    stk::mesh::MetaData * md_KOKKOS = eMesh_KOKKOS_NR.get_fem_meta_data();
    stk::mesh::Selector selector_KOKKOS (md_KOKKOS->universal_part());
    const stk::mesh::BucketVector & buckets_KOKKOS = eMesh_KOKKOS_NR.get_bulk_data()->get_buckets(
                            stk::topology::ELEMENT_RANK, selector_KOKKOS);

    stk::mesh::BucketVector::const_iterator iBucket_KOKKOS = buckets_KOKKOS.begin();
    for(;iBucket_KOKKOS!=buckets_KOKKOS.end();iBucket_KOKKOS++)
    {
        stk::mesh::Bucket & bucket_KOKKOS = **iBucket_KOKKOS;
        for(unsigned iEntity = 0; iEntity<bucket_KOKKOS.size();iEntity++){
            stk::mesh::Entity ent_KOKKOS = bucket_KOKKOS[iEntity];
            kokkosNR.doForAllSubEntities(&NodeRegistry_KOKKOS::registerNeedNewNode, ent_KOKKOS, needed_entity_ranks,0);
        }
    }

    kokkosNR.endRegistration();
    //DO REGISTRATION FOR NEW NR

    NodeRegistry boostNR(eMesh_BOOST_NR);
    boostNR.initialize();
    //DO REGISTRATION FOR OLD NR
    /*
     * 1st of three steps to create and associate new nodes - register need for new nodes, then check if node is remote, then get
     *   from remote proc if necessary; finally, the local node database is ready to be queried
     *
     * The pattern is to begin the step, loop over all elements (including ghosts) and invoke the local operation
     * The method doForAllSubEntities is a utility for performing the operation on all the sub entities.
     * If more granularity is desired, the member functions can be invoked directly for a particular sub-entity.
     */
    std::cout << "About to do registration for boostNR\n";
    boostNR.beginRegistration();

    stk::mesh::MetaData * md_BOOST = eMesh_BOOST_NR.get_fem_meta_data();
    stk::mesh::Selector selector_BOOST (md_BOOST->universal_part());
    const stk::mesh::BucketVector & buckets_BOOST = eMesh_BOOST_NR.get_bulk_data()->get_buckets(
                            stk::topology::ELEMENT_RANK, selector_BOOST);

    stk::mesh::BucketVector::const_iterator iBucket_BOOST = buckets_BOOST.begin();
    for(;iBucket_BOOST!=buckets_BOOST.end();iBucket_BOOST++)
    {
        stk::mesh::Bucket & bucket_BOOST = **iBucket_BOOST;
        for(unsigned iEntity = 0; iEntity<bucket_BOOST.size();iEntity++){
            stk::mesh::Entity ent_BOOST = bucket_BOOST[iEntity];
            boostNR.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, ent_BOOST, needed_entity_ranks,0);
        }
    }

    boostNR.endRegistration();
    //DO REGISTRATION FOR OLD NR

    unsigned no_boost_NOT_in_kokkos=0;
    SetOfEntities boost_U_kokkos;
    bool allInKokkos = boostNR.verifyAllKeysInKokkosNR(&kokkosNR,boost_U_kokkos,no_boost_NOT_in_kokkos);

    EXPECT_EQ(no_boost_NOT_in_kokkos,0u);
    EXPECT_EQ(allInKokkos,true);


    std::cout << std::endl << std::endl << std::endl;


    unsigned no_kokkos_NOT_in_boost=0;
    SetOfEntities kokkos_U_boost;
    bool allInBoost = kokkosNR.verifyAllKeysInBoostNR(&boostNR,kokkos_U_boost,no_kokkos_NOT_in_boost);

    EXPECT_EQ(no_kokkos_NOT_in_boost,0u);
    EXPECT_EQ(allInBoost,true);

    bool same_sets = kokkos_U_boost==boost_U_kokkos;

    EXPECT_EQ(same_sets,true);
}//TEST

TEST(nodeRegistry_KOKKOS, oldNR_v_newNR_create_prolongate_new_nodes)
{

}
} // namespace unit_tests
} // namespace percept

