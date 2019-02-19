// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CreateAdjacentEntities.hpp>

#include <percept/stk_rebalance/Rebalance.hpp>
#include <percept/stk_rebalance/Partition.hpp>
#include <percept/stk_rebalance/ZoltanPartition.hpp>

#include <percept/stk_rebalance_utils/RebalanceUtils.hpp>
#include <gtest/gtest.h>
#include <stk_mesh/base/MeshUtils.hpp>

static const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

typedef stk::mesh::Field<double> ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorField ;

enum { nx = 2, ny = 2 };

static const bool constraint_is_nodes = false;

TEST(UnitTestZoltanSimple, testUnit)
{
#ifdef STK_HAS_MPI
  stk::ParallelMachine comm(MPI_COMM_WORLD);
#else
  stk::ParallelMachine comm(0);
#endif

  unsigned spatial_dimension = 2;
  std::vector<std::string> rank_names = stk::mesh::entity_rank_names();
  const stk::mesh::EntityRank constraint_rank = static_cast<stk::mesh::EntityRank>(rank_names.size());
  rank_names.push_back("Constraint");

  stk::mesh::MetaData fem_meta;
  fem_meta.initialize( spatial_dimension, rank_names );
  stk::mesh::MetaData & meta_data = fem_meta;
  stk::mesh::BulkData bulk_data( meta_data , comm , stk::mesh::BulkData::AUTO_AURA );
  const stk::mesh::EntityRank element_rank    = stk::topology::ELEMENT_RANK;

  shards::CellTopology quad_top(shards::getCellTopologyData<shards::Quadrilateral<4> >());
  stk::mesh::Part & quad_part( fem_meta.declare_part("quad", stk::topology::ELEMENT_RANK ) );
  stk::mesh::set_cell_topology(quad_part, quad_top);

  VectorField & coord_field = fem_meta.declare_field< VectorField >( stk::topology::NODE_RANK,  "coordinates" ) ;
  ScalarField & weight_field = fem_meta.declare_field< ScalarField >( stk::topology::ELEMENT_RANK,  "element_weights" ) ;

  stk::mesh::put_field_on_mesh( coord_field ,  fem_meta.universal_part() , nullptr);
  stk::mesh::put_field_on_mesh(weight_field ,  fem_meta.universal_part() , nullptr);

  fem_meta.commit();

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  bulk_data.modification_begin();

  if ( p_rank == 0 ) {

    std::vector<std::vector<stk::mesh::Entity> > quads(nx);
    for ( unsigned ix = 0 ; ix < nx ; ++ix ) quads[ix].resize(ny);

    const unsigned nnx = nx + 1 ;
    const unsigned nny = ny + 1 ;
    stk::mesh::EntityIdVector nodes(4) ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk::mesh::Entity q = stk::mesh::declare_element( bulk_data , quad_part , elem , nodes );
        quads[ix][iy] = q;
      }
    }

    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::Entity  e = bulk_data.get_entity( element_rank, elem );
        double * const e_weight = stk::mesh::field_data( weight_field , e );
        *e_weight = 1.0;
      }
    }

    for ( unsigned iy = 0 ; iy <= ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix <= nx ; ++ix ) {
        stk::mesh::EntityId nid = 1 + ix + iy * nnx ;
        stk::mesh::Entity  n = bulk_data.get_entity( NODE_RANK, nid );
        double * const coord = stk::mesh::field_data( coord_field , n );
        coord[0] = .1*ix;
        coord[1] = .1*iy;
        coord[2] = 0;
      }
    }

    if (constraint_is_nodes)
      {
        const unsigned iy_left  =  0;
        const unsigned iy_right = ny;
        stk::mesh::PartVector add(1, &fem_meta.locally_owned_part());
        for ( unsigned ix = 0 ; ix <= nx ; ++ix ) {
          stk::mesh::EntityId nid_left  = 1 + ix + iy_left  * nnx ;
          stk::mesh::EntityId nid_right = 1 + ix + iy_right * nnx ;
          stk::mesh::Entity  n_left  = bulk_data.get_entity( NODE_RANK, nid_left  );
          stk::mesh::Entity  n_right = bulk_data.get_entity( NODE_RANK, nid_right );
          const stk::mesh::EntityId constraint_entity_id =  1 + ix + nny * nnx;
          stk::mesh::Entity  c = bulk_data.declare_constraint(constraint_entity_id, add );
          bulk_data.declare_relation( c , n_left  , 0 );
          bulk_data.declare_relation( c , n_right , 1 );
        }
      }
    else
      {
        stk::mesh::PartVector add(1, &fem_meta.locally_owned_part());
        for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
          for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
            stk::mesh::EntityId elem_id = 1 + ix + iy * nx ;
            stk::mesh::Entity  elem  = bulk_data.get_entity( element_rank, elem_id);
            const stk::mesh::EntityId constraint_entity_id =  elem_id + (ny+1) * (nx+1);
            stk::mesh::Entity  c = bulk_data.declare_constraint(constraint_entity_id, add );
            bulk_data.declare_relation( c , elem  , 0 );
          }
        }
      }

  }

  // Only P0 has any nodes or elements
  if ( p_rank == 0 ) {
    ASSERT_TRUE( ! bulk_data.buckets( NODE_RANK ).empty() );
    ASSERT_TRUE( ! bulk_data.buckets( element_rank ).empty() );
  }
  else {
    ASSERT_TRUE( bulk_data.buckets( NODE_RANK ).empty() );
    ASSERT_TRUE( bulk_data.buckets( element_rank ).empty() );
  }


  stk::mesh::fixup_ghosted_to_shared_nodes(bulk_data);
  bulk_data.modification_end();

  // create some sides and faces to rebalance.
  stk::mesh::PartVector add_parts;
  stk::mesh::create_adjacent_entities(bulk_data, add_parts);

  // Zoltan partition is specialized fomm a virtual base class, stk::rebalance::Partition.
  // Other specializations are possible.
  Teuchos::ParameterList emptyList;
  stk::rebalance::Zoltan zoltan_partition(bulk_data, comm, spatial_dimension, emptyList);

  //stk::mesh::EntityRank rank_to_rebalance = constraint_rank;
  stk::mesh::EntityRank rank_to_rebalance = element_rank;
  const double imbalance_threshold_before = stk::rebalance::check_balance(bulk_data, &weight_field, rank_to_rebalance);
  {
    stk::mesh::Selector selector(fem_meta.universal_part());

    stk::rebalance::Rebalance rb;
    rb.rebalance(bulk_data, selector, &coord_field, &weight_field, zoltan_partition, rank_to_rebalance);
  }

  const double imbalance_threshold = stk::rebalance::check_balance(bulk_data, &weight_field, rank_to_rebalance);
  const bool do_rebal = 1.5 < imbalance_threshold;

  // Check that we satisfy our threshhold
  std::cout << "imbalance_threshold_before= " << imbalance_threshold_before
            << " imbalance_threshold after= " << imbalance_threshold
            << std::endl;
  if (0)
    {
  ASSERT_TRUE( !do_rebal );
  if( (2 == p_size) || (4 == p_size) )
  {
    ASSERT_NEAR(imbalance_threshold, 1.0, 1.e-8);
  }
  else
  {
    ASSERT_LE(imbalance_threshold, 1.5);
  }
    }

  // And verify that all dependent entities are on the same proc as their parent element
  {
    stk::mesh::EntityVector entities;
    stk::mesh::Selector selector = fem_meta.locally_owned_part();

    get_selected_entities(selector, bulk_data.buckets(NODE_RANK), entities);
    if (rank_to_rebalance == element_rank)
      {
        bool result = stk::rebalance::verify_dependent_ownership(bulk_data, element_rank, entities);
        ASSERT_TRUE( result );
      }
    if (rank_to_rebalance == constraint_rank)
      {
        if (!constraint_is_nodes)
          get_selected_entities(selector, bulk_data.buckets(element_rank), entities);
        //std::cout << "num constraints= " << entities.size() << std::endl;
        bool result = stk::rebalance::verify_dependent_ownership(bulk_data, constraint_rank, entities);
        ASSERT_TRUE( result );
      }
  }
}

/// \page percept/stk_rebalance_unit_test_zoltan
///  \ingroup percept/stk_rebalance_unit_test_module
/// \section percept/stk_rebalance_unit_test_zoltan_description Simple Zoltan Unit Test
///
/// This unit test creates a 2D quad mesh on proc 0 with coordinates and
/// Parts associated with edges, nodes, and constraints and then moves these entities
/// as determined by calling Zoltan's load balancing capability using
/// default settings.
///
/// Using the following mesh of 4 quads,
///
///  Global node and element numbering
/// <pre>
///
///   7       8       9
///   +-------+-------+
///   |       |       |
///   |  e3   |  e4   |
///   |       |5      |
///  4+-------+-------+6
///   |       |       |     Y
///   |  e1   |  e2   |     |
///   |       |       |     |
///   +-------+-------+     *--> X
///   1       2      3
/// </pre>
///
///  Local node numbering
///
/// <pre>
///
///   3       4
///   +-------+
///   |       |
///   |  e1   |
///   |       |
///   +-------+
///   1       2
/// </pre>
///
/// Use of Zoltan with default settings is achieved by instantiating the
/// appropriate Partition class with an empty parameter list as follows,
/// \dontinclude UnitTestZoltanSimple.cpp
/// \skip Zoltan partition is
/// \until zoltan_partition(
///
/// An initial assessment of imbalance is made using an element weight field
/// followed by a call to actually do the rebalance as follows,
/// \skip Force a rebalance
/// \until rebalance::rebalance(
///
/// Perfect balancing should result using 2 or 4 procs, and
/// on 3 procs, the imbalance threshold should be below 1.5.
/// The test passes if these criteria are satisfied.
///
/// See \ref UnitTestZoltanSimple.cpp for the complete source listing.
///
