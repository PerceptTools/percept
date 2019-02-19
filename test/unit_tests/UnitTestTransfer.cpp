// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_io/StkMeshIoBroker.hpp>

#include <gtest/gtest.h>
#include <stk_util/diag/StringUtil.hpp>

#include <percept/xfer/STKMeshTransferSetup.hpp>
#include <percept/function/StringFunction.hpp>
#include <percept/function/FieldFunction.hpp>

#include <iostream>

namespace percept {
namespace unit_tests {

struct TransferFixture
{
  stk::io::StkMeshIoBroker mesh_data;
  stk::mesh::BulkData * bulkData;
  stk::mesh::Field<double, stk::mesh::Cartesian> * coordinates;
  stk::mesh::Field<double> * scalar_field;

  TransferFixture(stk::ParallelMachine comm,
		  size_t num_xyz = 4,
		  size_t num_y_arg=0,
		  size_t num_z_arg=0,
		  bool sidesets=false)
    : mesh_data(comm),
      bulkData(NULL),
      coordinates(NULL),
      scalar_field(NULL)
  {
    const size_t num_x = num_xyz;
    const size_t num_y = num_y_arg? num_y_arg : num_xyz;
    const size_t num_z = num_z_arg? num_z_arg : num_xyz;
    std::string config_mesh =
      sierra::to_string(num_x) + "x" +
      sierra::to_string(num_y) + "x" +
      sierra::to_string(num_z) + "|bbox:-0.5,-0.5,-0.5,0.5,0.5,0.5";
    if (sidesets) config_mesh += "|sideset:xXyYzZ";

    mesh_data.add_mesh_database(config_mesh, "generated", stk::io::READ_MESH);

    // Open, read, filter meta data from the input mesh file:
    // The coordinates field will be set to the correct dimension.
    mesh_data.create_input_mesh();

    // This defines all fields found on the input mesh as stk fields
    mesh_data.add_all_mesh_fields_as_input_fields();

    stk::mesh::MetaData * metaData = & mesh_data.meta_data();
    scalar_field = & (metaData->declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "scalar"));
    stk::mesh::put_field_on_mesh( *scalar_field, metaData->universal_part(), nullptr);

    mesh_data.populate_bulk_data();

    bulkData = & mesh_data.bulk_data();

    {
      StringFunction sf_scalar( "1+x+y+z+x*y+y*z+z*x" , Name("scalar"), 3, 1);
      FieldFunction ff_scalar("scalar", scalar_field, bulkData,
			      Dimensions(3), Dimensions(1));
      ff_scalar.interpolateFrom(sf_scalar);
    }

    coordinates = metaData->get_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");
  }

};

TEST(transfer, one_hex8_self)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if (1 < stk::parallel_machine_size(comm)) return;

  TransferFixture coarse_mesh (comm, 1);
  TransferFixture fine_mesh (comm, 1);

  {
    boost::shared_ptr<STKMeshTransfer> transfer =
      buildSTKMeshTransfer<STKMeshTransfer>(*coarse_mesh.bulkData,
			coarse_mesh.coordinates,
			coarse_mesh.scalar_field,
			*fine_mesh.bulkData,
			fine_mesh.coordinates,
			fine_mesh.scalar_field,
			"transfer_coords");

    initializeSTKMeshTransfer(&*transfer);

    transfer->apply();

    // TODO verify result
  }

  // TODO add vector/tensor fields for testing

}

}
}
