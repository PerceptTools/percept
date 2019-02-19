// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <regression_tests/RegressionTestExampleLocalRefine.hpp>

#include <gtest/gtest.h>

namespace percept
{
namespace regression_tests
{

TEST(example, two_tri_local_refiner)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD ;

    const unsigned p_size = stk::parallel_machine_size( pm );
    const unsigned p_rank = stk::parallel_machine_rank( pm );
    if (p_size>2) return;

    PerceptMesh eMesh;
    eMesh.open("two_tri3.g");

    Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);

    eMesh.register_and_set_refine_fields();

    Teuchos::RCP<UniformRefinerPatternBase> localBreakPattern = make_local_break_pattern(eMesh);

    eMesh.commit();

    eMesh.save_as("two_tri3_post_commit.g");

    // set up for TEA refine/unrefine
    stk::mesh::Selector univ_selector(eMesh.get_fem_meta_data()->universal_part());
    ElementRefinePredicate erp(eMesh, &univ_selector, eMesh.m_refine_field, 0.0);

    TransitionElementAdapter<ElementRefinePredicate> breaker(erp, eMesh, *localBreakPattern, 0, false);

    breaker.setRemoveOldElements(false);
    breaker.setAlwaysInitializeNodeRegistry(false);
    breaker.initializeRefine();

    // first refine element 1 (lower left)
    if (p_size==1 || p_rank==1)
        set_refine_field_to_value(eMesh, {1}, 1);

    eMesh.save_as("two_tri3_pre_refine1.e");

    breaker.refine();

    eMesh.save_as("two_tri3_post_refine1.e");

    // next refine element 6 (center element) and 7 (top - transition element)

    {
        std::ofstream mesh0("two_tri_post_refine1." + std::to_string(p_size) + "." + std::to_string(p_rank) + ".txt");
        stk::mesh::BulkData &bulk_data = *eMesh.get_bulk_data();
        bulk_data.dump_all_mesh_info(mesh0);
    }

    if (p_size==1)
        set_refine_field_to_value(eMesh, {6,7}, 1);
    else if (p_rank==0)
        set_refine_field_to_value(eMesh, {3}, 1);
    else
        set_refine_field_to_value(eMesh, {10}, 1);

    eMesh.save_as("two_tri3_pre_refine2.e");

    breaker.refine();

    eMesh.save_as("two_tri3_post_refine2.e");

    // set unrefine field on fully refined elements (2 levels)
    if (p_size==1)
        set_refine_field_to_value(eMesh, {11,12,13,14,21,22,23,24}, -1);
    else if (p_rank==0)
        set_refine_field_to_value(eMesh, {15,17,19,21}, -1);
    else
        set_refine_field_to_value(eMesh, {35,37,39,41}, -1);

    eMesh.save_as("two_tri3_pre_unrefine1.e");

    // unrefine
    breaker.unrefine();

    eMesh.save_as("two_tri3_post_unrefine1.e");

    // set unrefine field on fully refined elements
    if (p_size==1)
        set_refine_field_to_value(eMesh, {3,4,5,6}, -1);
    else if (p_rank==1)
        set_refine_field_to_value(eMesh, {4,6,8,10}, -1);

    // this seems required to get the full unrefinement
    stk::mesh::communicate_field_data(eMesh.get_bulk_data()->aura_ghosting(), {eMesh.m_refine_field});

    eMesh.save_as("two_tri3_pre_unrefine2.e");

    // unrefine
    breaker.unrefine();

    eMesh.save_as("two_tri3_post_unrefine2.e");
}

TEST(example, two_tri_local_manual)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD ;

    const unsigned p_size = stk::parallel_machine_size( pm );
    if (p_size!=1) return;

    PerceptMesh eMesh;
    eMesh.open("two_tri3.g");

    eMesh.commit();

    //eMesh.print_info("dummy", 5);

    stk::mesh::BulkData &bulk_data = *eMesh.get_bulk_data();
    stk::mesh::MetaData &meta_data = *eMesh.get_fem_meta_data();

    //std::ofstream mesh0("two_tri.txt");
    //bulk_data.dump_all_mesh_info(mesh0);

    // if the next line is removed the test fails
    bulk_data.delete_face_adjacent_element_graph();

    CoordinatesFieldType* coords_field = eMesh.get_coordinates_field();

    stk::mesh::PartVector elem_part_vec(1, meta_data.get_part("block_1"));

    stk::mesh::PartVector side_part1_vec(1, meta_data.get_part("surface_tri3_edge2_1"));
    stk::mesh::PartVector side_part2_vec(1, meta_data.get_part("surface_tri3_edge2_2"));
    stk::mesh::PartVector side_part3_vec(1, meta_data.get_part("surface_tri3_edge2_3"));
    stk::mesh::PartVector side_part4_vec(1, meta_data.get_part("surface_tri3_edge2_4"));

    stk::mesh::PartVector side_part_vec = {side_part1_vec[0],side_part2_vec[0],
            side_part3_vec[0],side_part4_vec[0]};

    // first refinement

    bulk_data.modification_begin();

    create_new_nodes(bulk_data, coords_field, {5,6,7},
            {{-1,0},{0,0},{0,-1}});

    stk::mesh::EntityVector all_nodes;
    gather_entities(bulk_data, stk::topology::NODE_RANK, all_nodes, 7);

    create_new_elems(bulk_data, all_nodes,
            {3,4,5,6,7,8},
            {{4, 7, 5},
                    {7, 1, 6},
                    {5, 6, 3},
                    {6, 5, 7},
                    {3, 6, 2},
                    {6, 1, 2}});

    stk::mesh::EntityVector all_elems;
    gather_entities(bulk_data, stk::topology::ELEMENT_RANK, all_elems, 8);

    //bulk_data.modification_end();
    //bulk_data.modification_begin();

    create_new_edges(bulk_data, side_part_vec, all_elems,
            {3,5,3,4},
            {2,2,0,0},
            {3,3,4,4});

    reconnect_old_edges(bulk_data,
            {23,21},
            {8,7},
            {1,2});

    bulk_data.modification_end();

    {
        std::ofstream mesh0("two_tri_refine1.txt");
        bulk_data.dump_all_mesh_info(mesh0);
    }

    eMesh.save_as("two_tri3_manual_post_refine1.e");

    // second refinement

    bulk_data.modification_begin();

    create_new_nodes(bulk_data, coords_field, {8,9,10,11,12},
            {{1,0},{-0.5,0},{0,-0.5},{-0.5,-0.5},{0,1}});

    gather_entities(bulk_data, stk::topology::NODE_RANK, all_nodes, 12);

    unrefine_transition_elements(bulk_data, {7,8});

    create_new_elems(bulk_data, all_nodes,
            {7,8,9,10,11,12,13,14,15,16,17,18,19,20},
            {{2, 12, 8},
                    {12, 3, 6},
                    {8, 6, 1},
                    {6, 8, 12},
                    {7, 11, 4},
                    {11, 5, 4},
                    {6, 10, 1},
                    {10, 7, 1},
                    {5, 9, 3},
                    {9, 6, 3},
                    {6, 9, 10},
                    {9, 5, 11},
                    {10, 11, 7},
                    {11, 10, 9}});

    gather_entities(bulk_data, stk::topology::ELEMENT_RANK, all_elems, 20);

    create_new_edges(bulk_data, side_part_vec, all_elems,
            {7,9,7,8},
            {2,2,0,0},
            {1,1,2,2});

    reconnect_old_edges(bulk_data,
            {31,41,33,53},
            {11,14,12,15},
            {2, 1, 1, 2});

    bulk_data.modification_end();

    {
        std::ofstream mesh0("two_tri_refine2.txt");
        bulk_data.dump_all_mesh_info(mesh0);
    }

    eMesh.save_as("two_tri3_manual_post_refine2.e");
}

TEST(DISABLED_example, two_tri_uniform_node_registry)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD ;

    const unsigned p_size = stk::parallel_machine_size( pm );
    //const unsigned p_rank = stk::parallel_machine_rank( pm );
    if (p_size>2) return;

    PerceptMesh eMesh;
    eMesh.open("two_tri3.g");

    // Q: which to use?
    //Local_Tri3_Tri3_N break_tri_to_tri_N(eMesh);
    //UniformRefinerPattern<shards::Triangle<3>, shards::Triangle<3>, 4, SierraPort > refiner_pattern(eMesh);
    RefinerPattern<shards::Triangle<3>, shards::Triangle<3>, -1 > refiner_pattern(eMesh);

    eMesh.commit();

    stk::mesh::BulkData &bulk_data = *eMesh.get_bulk_data();
    stk::mesh::MetaData &meta_data = *eMesh.get_fem_meta_data();

    //CoordinatesFieldType* coords_field = eMesh.get_coordinates_field();

    stk::mesh::PartVector elem_part_vec(1, meta_data.get_part("block_1"));

    stk::mesh::PartVector side_part1_vec(1, meta_data.get_part("surface_tri3_edge2_1"));
    stk::mesh::PartVector side_part2_vec(1, meta_data.get_part("surface_tri3_edge2_2"));
    stk::mesh::PartVector side_part3_vec(1, meta_data.get_part("surface_tri3_edge2_3"));
    stk::mesh::PartVector side_part4_vec(1, meta_data.get_part("surface_tri3_edge2_4"));

    stk::mesh::PartVector side_part_vec = {side_part1_vec[0],side_part2_vec[0],
            side_part3_vec[0],side_part4_vec[0]};

    stk::mesh::Selector locally_owned_selector(eMesh.get_fem_meta_data()->locally_owned_part());

    NodeRegistry nodeRegistry(eMesh);
    nodeRegistry.initialize();
    nodeRegistry.init_comm_all();

    std::vector<NeededEntityType> needed_entities;
    refiner_pattern.fillNeededEntities(needed_entities);
    EXPECT_EQ(needed_entities.size(), 1u);

    // first refinement
    {
        //bulk_data.modification_begin();

        stk::mesh::EntityVector elems;
        bulk_data.get_entities(stk::topology::ELEMENT_RANK, locally_owned_selector, elems);

        nodeRegistry.beginRegistration();

        //const bool needNodes = true;
        for (unsigned e=0; e<elems.size(); e++) {
            const stk::mesh::Entity element = elems[e];
            const CellTopologyData * const bucket_topo_data = eMesh.get_cell_topology(element);
            nodeRegistry.doForAllSubEntities(&NodeRegistry::registerNeedNewNode, element, needed_entities, bucket_topo_data);
        }

        // this call does a mod-begin
        nodeRegistry.endRegistration();

        // or just call this which is part of endRegistration
        //nodeRegistry.createNewNodesInParallel();

        // this does not work
        //nodeRegistry.dumpDB("NR.two_tri3_uniform_1.txt");

        nodeRegistry.beginCheckForRemote();

        for (unsigned e=0; e<elems.size(); e++) {
            const stk::mesh::Entity element = elems[e];
            const CellTopologyData * const bucket_topo_data = eMesh.get_cell_topology(element);

            nodeRegistry.doForAllSubEntities(&NodeRegistry::checkForRemote, element, needed_entities, bucket_topo_data);
        }

        nodeRegistry.endCheckForRemote();

        nodeRegistry.beginGetFromRemote();

        for (unsigned e=0; e<elems.size(); e++) {
            const stk::mesh::Entity element = elems[e];
            const CellTopologyData * const bucket_topo_data = eMesh.get_cell_topology(element);

            nodeRegistry.doForAllSubEntities(&NodeRegistry::getFromRemote, element, needed_entities, bucket_topo_data);
        }

        // this does a mod-end
        nodeRegistry.endGetFromRemote();

        //bulk_data.modification_end();
    }
    eMesh.save_as("two_tri3_uniform_NR_post_refine1.e");
}

} // namespace regression_tests
} // namespace percept

