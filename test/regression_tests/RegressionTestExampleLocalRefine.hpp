// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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

#include <percept/Percept.hpp>
#include <percept/Util.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/PerceptUtils.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>
#include <adapt/AdaptHelperFunctions.hpp>

#include <adapt/ElementRefinePredicate.hpp>
#include <adapt/TransitionElementAdapter.hpp>

namespace percept
{
namespace regression_tests
{

void set_refine_field_to_value(PerceptMesh &eMesh, std::vector<stk::mesh::EntityId> entity_list, const int refine_value)
{
  RefineFieldType *refine_field = eMesh.m_refine_field;

  for (unsigned i=0; i<entity_list.size(); i++) {
    stk::mesh::Entity entity = eMesh.get_bulk_data()->get_entity(stk::topology::ELEMENT_RANK,
                                                                 stk::mesh::EntityId(entity_list[i]));

    if (!eMesh.is_valid(entity)) {
        stk::ParallelMachine pm = MPI_COMM_WORLD ;

        const unsigned p_rank = stk::parallel_machine_rank( pm );
        std::cout << "P" << p_rank << ": cannot set value on entity ID " << entity_list[i] << std::endl;
        continue;
    }

    int *val= stk::mesh::field_data(*refine_field, entity);

    *val = refine_value;
  }
}

void create_new_nodes(stk::mesh::BulkData &bulk_data,
        CoordinatesFieldType *coords_field,
        const std::vector<stk::mesh::EntityId>  node_ids,
        const std::vector<std::vector<double>> node_coords)
{
    const unsigned num_new_nodes=node_ids.size();
    stk::mesh::PartVector empty;

    for (unsigned i=0; i<num_new_nodes; i++) {
        stk::mesh::Entity new_node = bulk_data.declare_node(node_ids[i], empty);
        double * coord = stk::mesh::field_data(*coords_field, new_node);
        coord[0] = node_coords[i][0];
        coord[1] = node_coords[i][1];
        // TODO call add_node_sharing on both procs
    }
}

void create_new_elems(stk::mesh::BulkData &bulk_data,
        const std::vector<stk::mesh::Entity>  all_nodes,
        const std::vector<stk::mesh::EntityId>  elem_ids,
        const std::vector<std::vector<int>> elem_connect)
{
    const unsigned num_new_elems=elem_ids.size();
    stk::mesh::PartVector elem_part_vec(1, bulk_data.mesh_meta_data().get_part("block_1"));

    {
        for (unsigned e=0; e<num_new_elems; e++) {
            stk::mesh::Entity new_elem = bulk_data.declare_element(elem_ids[e], elem_part_vec);
            for (int n=0; n<3; n++) {
                bulk_data.declare_relation(new_elem, all_nodes[elem_connect[e][n]-1], n);
            }
        }
    }
}

void create_new_edges(stk::mesh::BulkData &bulk_data,
        const stk::mesh::PartVector &side_part_vec,
        const std::vector<stk::mesh::Entity>  &all_elems,
        const std::vector<stk::mesh::EntityId> elem_ids,
        const std::vector<unsigned> elem_ordinals,
        const std::vector<unsigned> sideset_ids)
{
    const unsigned num_new_edges = elem_ids.size();

    for (unsigned e=0; e<num_new_edges; e++)
    {
        stk::mesh::Entity elem = bulk_data.get_entity(stk::topology::ELEMENT_RANK,
                stk::mesh::EntityId(elem_ids[e]));
        const unsigned ordinal = elem_ordinals[e];
        stk::mesh::PartVector partvec(1,side_part_vec[sideset_ids[e]-1]);

        bulk_data.declare_element_side(elem, ordinal, partvec);
    }
}

void reconnect_old_edges(stk::mesh::BulkData &bulk_data,
        std::vector<stk::mesh::EntityId> edge_ids,
        std::vector<stk::mesh::EntityId> new_elem_ids,
        std::vector<unsigned> new_ordinals)
{
    const unsigned num_old_edges = edge_ids.size();

    for (unsigned e=0; e<num_old_edges; e++) {

        const unsigned old_ordinal = (edge_ids[e] % 10) - 1; // 0-based
        const stk::mesh::EntityId old_elem_id = (edge_ids[e] - (edge_ids[e] % 10)) / 10; // 1-based

        stk::mesh::Entity edge     = bulk_data.get_entity(stk::topology::EDGE_RANK,
                stk::mesh::EntityId(    edge_ids[e]));
        stk::mesh::Entity old_elem = bulk_data.get_entity(stk::topology::ELEMENT_RANK,
                stk::mesh::EntityId(old_elem_id));

        bulk_data.destroy_relation(old_elem, edge, old_ordinal);

        stk::mesh::Entity new_elem = bulk_data.get_entity(stk::topology::ELEMENT_RANK,
                stk::mesh::EntityId(new_elem_ids[e]));

        bulk_data.declare_relation(new_elem, edge, new_ordinals[e]);
    }
}

void unrefine_transition_elements(stk::mesh::BulkData &bulk_data,
        std::vector<stk::mesh::EntityId> elems_to_delete)
{
    for (unsigned e=0; e<elems_to_delete.size(); e++) {
        stk::mesh::Entity elem = bulk_data.get_entity(stk::topology::ELEMENT_RANK,
                elems_to_delete[e]);
        bulk_data.destroy_entity(elem);
    }
}

void gather_entities(stk::mesh::BulkData &bulk_data, const stk::mesh::EntityRank rank, stk::mesh::EntityVector &entities, const unsigned max_id)
{
    entities.clear();

    for (unsigned i=0; i<max_id; i++) {
        entities.push_back(bulk_data.get_entity(rank, stk::mesh::EntityId(i+1)));
    }
}

}
}
