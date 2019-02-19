// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//new
#include <adapt/UniformRefinerPattern_ShellQuad4_ShellQuad4_4_sierra.hpp>

#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <percept/mesh/geometry/kernel/GeometryKernelPGEOM.hpp>
#include <percept/PerceptMesh.hpp>

#include "CubitVector.hpp"
#include "PGeomAssoc.hpp"
#include "PGeom.hpp"
#ifdef HAVE_ACIS
#include "PGeomACIS.hpp"
#endif

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/CreateEdges.hpp>

// this is for testing the local-refine refactoring
#define UNIFORM_REFINER UniformRefiner

namespace percept {
namespace regression_tests {

#include "RegressionTestFileLoc.hpp"

std::string to_string(int x) {
	std::string translated = "";
	std::string nums = "0123456789";

	int y = x;
	int divisor = 10;
	int index = 0;

	while (y != 0) {
		index = (y % divisor);
		y = ((y - index) / 10);
		std::string num = nums.substr(index, 1);
		translated = num + translated;
	}
	return translated;
}

bool validate_node_owner_ship (const stk::mesh::BulkData *mesh, const std::vector<int> &nodeIDsToCheck)
{
	stk::mesh::Selector localOrShared(mesh->mesh_meta_data().locally_owned_part());
	localOrShared |=mesh-> mesh_meta_data().globally_shared_part();

	for(unsigned i=0;i<nodeIDsToCheck.size(); i++)
	{
		stk::mesh::Entity nodeToCheck = mesh->get_entity(stk::topology::NODE_RANK, (stk::mesh::EntityId) nodeIDsToCheck[i]);
		if (nodeToCheck==stk::mesh::Entity::InvalidEntity){
			return false;
		}

		stk::mesh::Bucket &node_bucket = mesh->bucket(nodeToCheck);
		if (!localOrShared(node_bucket))
			return false;

	}

	return true;
}

//==================================
bool get_node_from_id(const stk::mesh::BulkData *mesh, int node_id,
		stk::mesh::Entity &node) {
	node = mesh->get_entity(stk::topology::NODE_RANK,
			(stk::mesh::EntityId) node_id);
	return (node == stk::mesh::Entity::InvalidEntity) ? false : true;


}
//===================================
bool get_beam_from_ids(const stk::mesh::BulkData *mesh,
		std::vector<int> edge_node_ids, stk::mesh::Entity &beam) {

	std::vector<stk::mesh::Entity> edgeNodes(2);
	for (unsigned i = 0; i < edgeNodes.size(); i++) {
		get_node_from_id(mesh, edge_node_ids[i], edgeNodes[i]);
	}

	stk::mesh::Selector localAndShared(
			mesh->mesh_meta_data().locally_owned_part());
	localAndShared |= mesh->mesh_meta_data().globally_shared_part();
	const stk::mesh::BucketVector & buckets = mesh->get_buckets(
			stk::topology::ELEMENT_RANK, localAndShared);


	stk::mesh::BucketVector::const_iterator k = buckets.begin();


	for (; k != buckets.end(); k++) {

		stk::mesh::Bucket & bucket = **k;

		if(bucket.topology()!=stk::topology::BEAM_2)
			continue;

		const unsigned num_elements_in_bucket = bucket.size();
		unsigned iEntity = (num_elements_in_bucket - 1);

		for (iEntity = 0; iEntity < num_elements_in_bucket; iEntity++) {

			stk::mesh::Entity beamToInspect = bucket[iEntity];

			const MyPairIterRelation beam_nodes(*mesh, beamToInspect,
					stk::topology::NODE_RANK);

			bool found = true;
			for (unsigned jj = 0; jj < beam_nodes.size(); jj++) {
				if (beam_nodes[jj].entity() != edgeNodes[0]
						&& beam_nodes[jj].entity() != edgeNodes[1]) {
					found = false;
					break;
				}
			}

			if (found) {

				beam = beamToInspect;
				return true;
			}
		}
	}
	beam = stk::mesh::Entity::InvalidEntity;
	return false;
}

//===================================
bool get_shell_from_ids(const stk::mesh::BulkData *mesh,
		std::vector<int> face_node_ids, stk::mesh::Entity &shell) {
	std::vector<stk::mesh::Entity> faceNodes(4);

	for (unsigned i = 0; i < faceNodes.size(); i++) {
		get_node_from_id(mesh, face_node_ids[i], faceNodes[i]);
	}

	//does this reverse the filtering of ghosted nodes? As in, will we
	unsigned numElements1 = mesh->num_elements(faceNodes[0]); //i assume: number of elements that contain a specific node
	const stk::mesh::Entity * elements_1 = mesh->begin_elements(faceNodes[0]); //i assume: pointer to list of elements that contain all specific node


	//QUESITON: IS THERE AN EFFICIENT MEANS TO ONLY GET SELECTIONS THAT CONTAIN ALL FOUR NODES?

	//THIS IS FASTER!
	for(unsigned e = 0; e < numElements1; ++e) { //loop over elements
	   stk::mesh::Entity element = elements_1[e];
	   // test if this matches the nodes ...


	   //filter out topology WITH BUCKET?, nah just check the size of face_nodes. If it's too big, you got the wrong kind of element.
	   //However, my pair relation is expensive, so we might want to try and find a better means of figuring out element topology
	   //not sure how expensive bucketing is, through.

	   stk::mesh::Bucket & node_bucket = mesh->bucket(element);
	   if (node_bucket.topology()!=stk::topology::SHELL_QUAD_4 && node_bucket.topology()!=stk::topology::SHELL_TRI_3)
		   continue;

	   const MyPairIterRelation face_nodes(*mesh, element, stk::topology::NODE_RANK); //want to compute this as little as possible. I think it's cheaper than bucketing
//	   if (face_nodes.size()!=4)//4 would be the number of faces in a quad_sell
//		   continue;

		bool found = true;
		for (unsigned jj = 0; jj < face_nodes.size(); jj++) {

//			if(face_nodes.size()==3)
//				if

			if (	face_nodes[jj].entity() != faceNodes[0]
					&& face_nodes[jj].entity() != faceNodes[1]
					&& face_nodes[jj].entity() != faceNodes[2]
					&& face_nodes[jj].entity() != faceNodes[3]) {
				found = false; //if you can't find it, look at the next element
				break;
			}
		}

		if (found) {

			shell = element;
			return true;
		}


	}

	shell = stk::mesh::Entity::InvalidEntity;
	return false;
} //end getshell

bool get_EDGE_from_ids(const stk::mesh::BulkData *mesh,
		std::vector<int> edge_node_ids, stk::mesh::Entity &face) {
	std::vector<stk::mesh::Entity> edgeNodes(2);
	for (unsigned i = 0; i < edgeNodes.size(); i++) {
		get_node_from_id(mesh, edge_node_ids[i], edgeNodes[i]);
	}
	//iterate through nodes, see if they'ye in a beam

	stk::mesh::Selector localAndShared(
			mesh->mesh_meta_data().locally_owned_part());
	localAndShared |= mesh->mesh_meta_data().globally_shared_part();
	const stk::mesh::BucketVector & buckets = mesh->get_buckets(
			stk::topology::EDGE_RANK, localAndShared);


	stk::mesh::BucketVector::const_iterator k = buckets.begin();


	for (; k != buckets.end(); k++) {

		stk::mesh::Bucket & bucket = **k;

		const unsigned num_elements_in_bucket = bucket.size();
		unsigned iEntity = (num_elements_in_bucket - 1);

		for (iEntity = 0; iEntity < num_elements_in_bucket; iEntity++) {

			stk::mesh::Entity faceToInspect = bucket[iEntity];

			const MyPairIterRelation beam_nodes(*mesh, faceToInspect,
					stk::topology::NODE_RANK);

			bool found = true;
			for (unsigned jj = 0; jj < beam_nodes.size(); jj++) {
				if (	beam_nodes[jj].entity() != edgeNodes[0]
						&& beam_nodes[jj].entity() != edgeNodes[1]) {
					found = false;
					break;
				}
			}

			if (found) {

				face = faceToInspect;
				return true;
			}
		}
	}

	face = stk::mesh::Entity::InvalidEntity;
	return false;
} //end getEDGE

bool get_FACE_from_ids(const stk::mesh::BulkData *mesh,
		std::vector<int> face_node_ids, stk::mesh::Entity &shell) {
	std::vector<stk::mesh::Entity> faceNodes(4);
	for (unsigned i = 0; i < faceNodes.size(); i++) {
		get_node_from_id(mesh, face_node_ids[i], faceNodes[i]);
	} //is issue is that it's giving me IDs from the .m2g file that don't correspond to anything on this particular mesh.


	stk::mesh::Selector localAndShared(
			mesh->mesh_meta_data().locally_owned_part());
	localAndShared |= mesh->mesh_meta_data().globally_shared_part(); //if I take this off change_enity parts errors out because I have face I'm not allowed to access
	const stk::mesh::BucketVector & bucketsOfFaces = mesh->get_buckets(
		stk::topology::FACE_RANK, localAndShared);

	stk::mesh::BucketVector::const_iterator faceBucketpntr =
			bucketsOfFaces.begin();
	for (; faceBucketpntr != bucketsOfFaces.end(); faceBucketpntr++) { //maybe optimize by creating a diminshing part

		stk::mesh::Bucket & bucket = **faceBucketpntr;

		const unsigned num_elements_in_bucket = bucket.size();

		for (unsigned iEntity = 0; iEntity < num_elements_in_bucket;
				iEntity++) {

			stk::mesh::Entity faceToInspect = bucket[iEntity];

			const MyPairIterRelation face_nodes(*mesh, faceToInspect,
					stk::topology::NODE_RANK);

			bool found = true;
			for (unsigned jj = 0; jj < face_nodes.size(); jj++) {
				if (face_nodes[jj].entity() != faceNodes[0]
						&& face_nodes[jj].entity() != faceNodes[1]
						&& face_nodes[jj].entity() != faceNodes[2]
						&& face_nodes[jj].entity() != faceNodes[3]) {
					found = false;
					break;
				}
			}

			if (found) {

				shell = faceToInspect;
				return true;
			}
		}
	}



	shell = stk::mesh::Entity::InvalidEntity;
	return false;



}//end getface


//===================================
bool get_hex_from_id(const stk::mesh::BulkData *mesh, int hex_id,
		stk::mesh::Entity &hex) {
	hex = mesh->get_entity(stk::topology::ELEM_RANK,
			(stk::mesh::EntityId) hex_id);
	return (hex == stk::mesh::Entity::InvalidEntity) ? false : true;
} //end gethex

TEST(DISABLED_DISABLED_ParallelTest, CCS){
	std::string meshName = "coarseCylinderSplitCurved.g";


	percept::PerceptMesh eMesh(3u);
	eMesh.open(meshName);
	stk::mesh::MetaData * md = eMesh.get_fem_meta_data();


	stk::mesh::Part & skin_part = md->declare_part("SkinPart", stk::topology::FACE_RANK); //doing this before the break pattern should make this work

	int scalarDimension = 0;
	stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank",
			stk::topology::ELEMENT_RANK, scalarDimension);
	URP_Heterogeneous_3D break_pattern(eMesh);



	eMesh.commit();
//	eMesh.get_bulk_data()->modification_begin();doesn't help
	std::vector< stk::mesh::Part * > partToPutSidesInto(1,&skin_part);
	stk::mesh::Selector blocksToSkinOrConsider(md->locally_owned_part());
	blocksToSkinOrConsider |= (md->globally_shared_part());
	stk::mesh::create_exposed_block_boundary_sides(*eMesh.get_bulk_data(), blocksToSkinOrConsider, partToPutSidesInto);
	stk::mesh::create_interior_block_boundary_sides(*eMesh.get_bulk_data(),blocksToSkinOrConsider, partToPutSidesInto);
//	eMesh.get_bulk_data()->modification_end();doesn't help


	std::cout<<"Refining..."<<std::endl;
	UNIFORM_REFINER breaker(eMesh,break_pattern,proc_rank_field);
	breaker.setRemoveOldElements(true);
	breaker.setAvoidFixSideSets(true);
	breaker.doBreak();



	eMesh.save_as("output_ParallelTest_CCS.g");



}

TEST(DISABLED_numFaces, TEST){ //Seems to run fine in serial. STILL NEED TO FIGURE OUT HOW TO SET DEFAULT VALIDATE CALLBACK!
	percept::PerceptMesh eMesh(3u);
	std::string meshToOpen;
	std::cin >> meshToOpen;
	eMesh.open(meshToOpen);


	stk::mesh::MetaData * md = eMesh.get_fem_meta_data();

	stk::mesh::Part & skin_part = md->declare_part("SkinPart", stk::topology::FACE_RANK); //doing this before the break pattern should make this work
	eMesh.commit();

	std::vector< stk::mesh::Part * > partToPutSidesInto(1,&skin_part);
	stk::mesh::Selector blocksToSkinOrConsider(md->locally_owned_part());
	blocksToSkinOrConsider |= (md->globally_shared_part());
	stk::mesh::create_exposed_block_boundary_sides(*eMesh.get_bulk_data(), blocksToSkinOrConsider, partToPutSidesInto);
	stk::mesh::create_interior_block_boundary_sides(*eMesh.get_bulk_data(),blocksToSkinOrConsider, partToPutSidesInto);
	create_edges(*eMesh.get_bulk_data());
	stk::mesh::Selector selector (md->universal_part());
	const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->get_buckets(
						stk::topology::FACE_RANK, selector);

	stk::mesh::BucketVector::const_iterator k2 = buckets.begin();

	size_t numFace =0;
	for(;k2 != buckets.end();k2++){
		stk::mesh::Bucket & bucket = **k2;
		numFace = bucket.size() + numFace;
	}


	std::cout << numFace <<std::endl;
}


#ifdef HAVE_ACIS
TEST(DISABLED_PGeomAssoc, NEWIMP){ //Seems to run fine in serial. STILL NEED TO FIGURE OUT HOW TO SET DEFAULT VALIDATE CALLBACK!
	percept::PerceptMesh eMesh(3u);
	eMesh.open("multipleGeomNoBlocksNoSS_2.g");
	std::string geomfile = "multipleGeomNoBlocksNoSS_2.sat";
	std::string assocfile = "multipleGeomNoBlocksNoSS_2.m2g";
	PGeomACIS pg;
	pg.initialize(ACIS_GEOMETRY_ENGINE);
	bool acis_read = pg.import_acis_file(geomfile.c_str());
	EXPECT_EQ(true, acis_read);
	PGeom* pgeom = &pg;
	stk::mesh::MetaData * md = eMesh.get_fem_meta_data();

	stk::mesh::Part & skin_part = md->declare_part("SkinPart", stk::topology::FACE_RANK); //doing this before the break pattern should make this work
	eMesh.commit();

	std::vector< stk::mesh::Part * > partToPutSidesInto(1,&skin_part);
	stk::mesh::Selector blocksToSkinOrConsider(md->locally_owned_part());
	blocksToSkinOrConsider |= (md->globally_shared_part());
	stk::mesh::create_exposed_block_boundary_sides(*eMesh.get_bulk_data(), blocksToSkinOrConsider, partToPutSidesInto);
	stk::mesh::create_interior_block_boundary_sides(*eMesh.get_bulk_data(),blocksToSkinOrConsider, partToPutSidesInto);
	create_edges(*eMesh.get_bulk_data());


	PGeomAssoc<stk::mesh::BulkData, stk::mesh::Entity, stk::mesh::Entity,
			stk::mesh::Entity, stk::mesh::Entity> geom_assoc(pgeom); //
	geom_assoc.set_node_callback(get_node_from_id);
	geom_assoc.set_edge_callback(get_EDGE_from_ids);
	geom_assoc.set_face_callback(get_FACE_from_ids);
	geom_assoc.set_elem_callback(get_hex_from_id);

	geom_assoc.set_mesh(eMesh.get_bulk_data());

	std::cout<<"Importing..." <<std::endl;
        const bool geometry_exists = true;
	bool read_assoc = geom_assoc.import_m2g_file(assocfile.c_str(), geometry_exists );
	std::cout<<"     Done importing" << std::endl;
	EXPECT_EQ(true, read_assoc);




	std::vector<int> surf_ids;
	pgeom->get_surfaces(surf_ids);
	std::vector<int> curve_ids;
	pgeom->get_curves(curve_ids);
	std::vector<int> vol_ids;
	pgeom->get_volumes(vol_ids);
	std::vector<int> vert_ids;
	pgeom->get_vertices(vert_ids);


	std::cout<<"Here are your vol/elems: " << std::endl;
	for (unsigned i=0;i<vol_ids.size();i++){
		std::vector<stk::mesh::Entity> elems= geom_assoc.get_vol_elems(vol_ids[i]);
		std::cout << "    " << elems.size() << std::endl;
	}
	std::cout<<"Here are your surfaces/faces: " << std::endl;
	for (unsigned i=0;i<surf_ids.size();i++){
		std::vector<stk::mesh::Entity> faces= geom_assoc.get_surf_faces(surf_ids[i]);
		std::cout << "    " << faces.size() << std::endl;
	}
	std::cout<<"Here are your curves/edges: " << std::endl;
	for (unsigned i=0;i<curve_ids.size();i++){
		std::vector<stk::mesh::Entity> edges= geom_assoc.get_curv_edges(curve_ids[i]);
		std::cout << "    " << edges.size() << std::endl;
	}
	std::cout<<"Here are your vol/nodes: " << std::endl;
	for (unsigned i=0;i<vol_ids.size();i++){
		std::vector<stk::mesh::Entity> nodes= geom_assoc.get_vol_nodes(vol_ids[i]);
		std::cout << "    " << nodes.size() << std::endl;
	}
	std::cout<<"Here are your vert/nodes: " << std::endl;
	for (unsigned i=0;i<vert_ids.size();i++){
		std::vector<stk::mesh::Entity> verts= geom_assoc.get_vert_nodes(vol_ids[i]);
		std::cout << "    " << verts.size() << std::endl;
	}
}




TEST(DISABLED_ManualBeamsShells, BMS) {

//	std::cout<< "Enter a file name with NO extension: " << std::endl;
//	std::string fileName = "cylinderSplitHollow"; //ERROR: mpirun noticed that process rank 0 with PID 23020 on node ceerws2901c exited on signal 9 (Killed). USING TOO MANY RESOURCES
//	std::string fileName = "coarseCylinderSplitCurved_8_ELEMS";
	std::string fileName = "coarseBigMesh";
//	std::string fileName = "TETMESHcylinderSplit";
//	std::string fileName = "TETMESHcube";
//	std::string fileName = "coarseCylinderSplit";
//	std::string fileName = "coarseCylinderSplitCurved";
//	std::string fileName = "multipleGeomNoBlocksNoSS_2";
//	std::string fileName = "coarseCylinderSplitHollow"; //this works though
//	cin >> fileName;


//	int THIS_PROC_ID = stk::parallel_machine_rank( MPI_COMM_WORLD);
//	std::ofstream output_debug("ManualBeamShells_BMS"+std::to_string(THIS_PROC_ID));
//	output_debug << "PROCESSOR " << THIS_PROC_ID<< std::endl;
//	std::ofstream output_debug_2("call_back_debug"+std::to_string(THIS_PROC_ID), std::ofstream::out | std::ofstream::app);
	std::string meshName = fileName + ".g";
	std::string assocfile = fileName + ".m2g";
	std::string geomfile = fileName + ".sat";

	percept::PerceptMesh eMesh(3u);
	eMesh.open(meshName);
	stk::mesh::MetaData * md = eMesh.get_fem_meta_data();


	PGeomACIS pg;
	pg.initialize(ACIS_GEOMETRY_ENGINE);
	bool acis_read = pg.import_acis_file(geomfile.c_str());
	EXPECT_EQ(true, acis_read);
	PGeom* pgeom = &pg;


	std::vector<int> surfIDs;
	pgeom->get_surfaces(surfIDs);
	std::vector<std::string> quadNames(surfIDs.size());
	std::vector<std::string> triNames(surfIDs.size());

	for(unsigned i = 0; i < surfIDs.size(); i++) { //make parts that we'll store new mesh on
		std::string name = "geom_surface_QUAD_";
		name = name + std::to_string(surfIDs[i]);
		stk::mesh::Part& part = md->declare_part_with_topology(name,
				stk::topology::SHELL_QUAD_4);
		stk::io::put_io_part_attribute(part); //set them up so they show up on output
		quadNames[i] = name;
	} //endfor


	for(unsigned i = 0; i<surfIDs.size();i++){
		std::string name = "geom_surface_TRI_";
		name = name + std::to_string(surfIDs[i]);
		stk::mesh::Part& part = md->declare_part_with_topology(name,
				stk::topology::SHELL_TRI_3);
		stk::io::put_io_part_attribute(part); //set them up so they show up on output
		triNames[i] = name;

	}



	std::vector<int> curveIDs;
	pgeom->get_curves(curveIDs);
	std::vector<std::string> curveNames(curveIDs.size());
	for (unsigned i = 0; i < curveIDs.size(); i++) { //make parts that we'll store new mesh on
		std::string name = "geom_curve_";
		name = name + std::to_string(curveIDs[i]);
		stk::mesh::Part& part = md->declare_part_with_topology(name,
				stk::topology::BEAM_2);
		stk::io::put_io_part_attribute(part); //set them up so they show up on output
		curveNames[i] = name;
	} //endfor

	//setup refinement
	int scalarDimension = 0;
	stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank",
			stk::topology::ELEMENT_RANK, scalarDimension);
	URP_Heterogeneous_3D break_pattern(eMesh);


	eMesh.commit();


	stk::mesh::BulkData * bd = eMesh.get_bulk_data();
	eMesh.initializeIdServer(); //initialize the ID server for creation of unique IDs for beam/shell parts

	stk::mesh::Selector localAndSharedPart = (md->locally_owned_part()); //these still contain ghosted_shared peices of mesh along boundary
			localAndSharedPart |= md->globally_shared_part();


	//BEGIN ASSOCATION
	PGeomAssoc<stk::mesh::BulkData, stk::mesh::Entity, stk::mesh::Entity,
			stk::mesh::Entity, stk::mesh::Entity> geom_assoc(pgeom); //
	geom_assoc.set_node_callback(get_node_from_id);
	geom_assoc.set_edge_callback(get_beam_from_ids);
	geom_assoc.set_face_callback(get_shell_from_ids);
	geom_assoc.set_elem_callback(get_hex_from_id);
	geom_assoc.set_validate_nodes_callback(validate_node_owner_ship);
	geom_assoc.set_mesh(eMesh.get_bulk_data());
	geom_assoc.set_fill_curve_and_surface_maps_during_import(false);
	bool read_assoc = geom_assoc.import_m2g_file(assocfile.c_str(), true);
	EXPECT_EQ(true, read_assoc);
	std::cout<<"FINISHED FIRST PART OF ASSOCIATION"<<std::endl;



	eMesh.get_bulk_data()->modification_begin();
	/////////////////////////////////////////////
	int num_beams=0;
//	output_debug <<"Num curves is: "<< curveIDs.size() <<std::endl;
	for(unsigned i = 0; i < curveIDs.size(); i++) { //create beams and put them into corresponding curve parts

		std::vector<stk::mesh::Part *> add_parts_beams(1,
				static_cast<stk::mesh::Part*>(0));

		add_parts_beams[0] = md->get_part(curveNames[i]);

		std::vector<std::vector<int>> edge_node_ids;
		geom_assoc.get_curve_edge_nodes(curveIDs[i], edge_node_ids);

		for (unsigned ii = 0; ii < edge_node_ids.size(); ii++) {

			bool toDeclare = true;

			std::vector<stk::mesh::EntityId> beamNodeIDs;
			for (unsigned j = 0; j < edge_node_ids[ii].size(); j++)
				beamNodeIDs.push_back((stk::mesh::EntityId) edge_node_ids[ii][j]);

			for (unsigned j = 0; j < edge_node_ids[ii].size(); j++) {



						stk::mesh::Entity cur_node = bd->get_entity(
								stk::topology::NODE_RANK, beamNodeIDs[j]); //CHECK IF NODE IS INVALID
						if ( !bd->is_valid(cur_node) ){
							toDeclare=false;
							break;
						}


						stk::mesh::Bucket &node_bucket = eMesh.bucket(cur_node);
						if (!localAndSharedPart(node_bucket)) //errors (sefaults) out here on proc 0 on first iteration
							toDeclare = false;

					}

					if (toDeclare) {

						stk::mesh::EntityId id2 = eMesh.getNextId(
								stk::topology::ELEMENT_RANK);
						stk::mesh::declare_element(*eMesh.get_bulk_data(),
								add_parts_beams, id2, beamNodeIDs);
						num_beams++;



					}

				} //endfor

			} //endfor
	std::cout<<"Created beams"<< std::endl;




	stk::mesh::Selector sharedOnly(md->globally_shared_part());
	sharedOnly -= md->locally_owned_part();


	for (unsigned i = 0; i < surfIDs.size(); i++) {
		std::vector<stk::mesh::Part *> add_parts_shells(1,
				static_cast<stk::mesh::Part*>(0));



		std::vector<std::vector<int>> face_node_ids;
		geom_assoc.get_surface_face_nodes(surfIDs[i], face_node_ids);

		for (unsigned ii = 0; ii < face_node_ids.size(); ii++) {




			std::vector<stk::mesh::EntityId> shellNodeIDs;
			for (unsigned j = 0; j < face_node_ids[ii].size(); j++)
				shellNodeIDs.push_back(
						(stk::mesh::EntityId) face_node_ids[ii][j]);

			bool isTri = true;
			if(shellNodeIDs.size()==3)
				add_parts_shells[0] = md->get_part(triNames[i]); //store these parts in a partvector for FASTER access, later
			else{
				add_parts_shells[0] = md->get_part(quadNames[i]);
				isTri=false;
			}
			bool toDeclare = true;
			//call function to see if you own/share the nodes only then create the beams

			int num_shared_nodes = 0;
			for (unsigned j = 0; j < face_node_ids[ii].size(); j++) {

				stk::mesh::Entity cur_node = bd->get_entity( //CHECK IF NODE IS INVALID
						stk::topology::NODE_RANK, shellNodeIDs[j]);
				if ( !bd->is_valid(cur_node) ){
					toDeclare=false; //if the node resides on a different mesh, don't try to create a shell out of it
					break;
				}
				stk::mesh::Bucket &node_bucket = eMesh.bucket(cur_node);
				if(sharedOnly(node_bucket)){
					num_shared_nodes++;
					if (( num_shared_nodes == 4 && !isTri ) || (num_shared_nodes == 3 && isTri) )
						toDeclare = false; //if all of the nodes are shared, don't create the shell
				}
				if (!localAndSharedPart(node_bucket))
					toDeclare = false;
			}



			if (toDeclare) {

				stk::mesh::EntityId id2 = eMesh.getNextId(
						stk::topology::ELEMENT_RANK);

				stk::mesh::declare_element(*eMesh.get_bulk_data(),
						add_parts_shells, id2, shellNodeIDs);



			}

		} //endfor

	}
	std::cout<<"Created faces" << std::endl;

	eMesh.get_bulk_data()->modification_end();

	geom_assoc.fill_curve_and_surface_maps();
	std::cout<<"filled maps" << std::endl;

	eMesh.save_as("output_ManualBeamsShells_BMS_INTERMEDIATEMESH.g");//this is a valid mesh that refines when saved out



	UNIFORM_REFINER breaker(eMesh,break_pattern,proc_rank_field);
	breaker.setRemoveOldElements(true);
	breaker.doBreak();


	// FOR SOME REASON THIS SNAPPING LOGIC SEVERELY DEFORMS SINGLE HEX CASE!!!!


	eMesh.save_as("output_ManualBeamsShells_BMSR1_NOTSNAPPED.g");

		const stk::mesh::BucketVector & buckets4 =
				eMesh.get_bulk_data()->get_buckets(stk::topology::NODE_RANK,
						localAndSharedPart);

		stk::mesh::BucketVector::const_iterator k4 = buckets4.begin();

	for (; k4 != buckets4.end(); k4++) {
		std::vector<int> isInCurves;
		std::vector<int> isInSurfs;
		stk::mesh::Bucket & bucket2 = **k4;

		for (unsigned i = 0; i < curveIDs.size(); i++) {

			stk::mesh::Selector boundarySelectorInner(
					*md->get_part(curveNames[i]));

			if (boundarySelectorInner(bucket2)) //if the bucket is in our part
			{
				isInCurves.push_back(curveIDs[i]); //add it to our evalutors
			}

		}

		for (unsigned i = 0; i < surfIDs.size(); i++) {
			stk::mesh::Selector boundarySelectorInner(
					*md->get_part(quadNames[i]));

			if (boundarySelectorInner(bucket2))
				isInSurfs.push_back(surfIDs[i]);
		}

		for (unsigned i = 0; i < surfIDs.size(); i++) {
			stk::mesh::Selector boundarySelectorInner(
					*md->get_part(triNames[i]));

			if (boundarySelectorInner(bucket2))
				isInSurfs.push_back(surfIDs[i]);
		}

		if (isInCurves.size() == 0 && isInSurfs.size() == 0)
			continue;

		const unsigned num_elements_in_bucket2 = bucket2.size();
		for (unsigned iEntity = 0; iEntity < num_elements_in_bucket2;
				iEntity++) {

			stk::mesh::Entity node = bucket2[iEntity];

			double * data = eMesh.field_data(*eMesh.get_coordinates_field(),
					node);

			CubitVector old_locations(data[0], data[1], data[2]);

			if (isInCurves.size() > 0) {
				double smallest_dist_sqr = std::numeric_limits<double>::max();
				for (unsigned ii = 0; ii < isInCurves.size(); ii++) {
					CubitVector new_locations = pgeom->get_curv_closest_point(
							isInCurves[ii], old_locations);

					double dist_sqr = pow(
							(old_locations.x() - new_locations.x()), 2)
							+ pow((old_locations.y() - new_locations.y()), 2)
							+ pow((old_locations.z() - new_locations.z()), 2);

					if (dist_sqr < smallest_dist_sqr) {
						smallest_dist_sqr = dist_sqr;

						data[0] = new_locations.x();
						data[1] = new_locations.y();
						data[2] = new_locations.z();
					}

				}

			}

			else {
				double smallest_dist_sqr = std::numeric_limits<double>::max();
				for (unsigned ii = 0; ii < isInSurfs.size(); ii++) {
					CubitVector new_locations = pgeom->get_surf_closest_point(
							isInSurfs[ii], old_locations);

					double dist_sqr = pow(
							(old_locations.x() - new_locations.x()), 2)
							+ pow((old_locations.y() - new_locations.y()), 2)
							+ pow((old_locations.z() - new_locations.z()), 2);

					if (dist_sqr < smallest_dist_sqr) {
						smallest_dist_sqr = dist_sqr;
						data[0] = new_locations.x();
						data[1] = new_locations.y();
						data[2] = new_locations.z();
					}

				}

			}

		}

	}
	eMesh.save_as("output_ManualBeamsShells_BMSR1.g");

}//manualBeamShells



TEST(DISABLED_ManualBeamsFaces, BMS) {

//	std::cout<< "Enter a file name with NO extension: " << std::endl;
	std::string fileName = "multipleGeomNoBlocksNoSS_2";
//	cin >> fileName;


	int THIS_PROC_ID = stk::parallel_machine_rank( MPI_COMM_WORLD);
	std::ofstream output_debug("ManualBeamFaces_BMS"+std::to_string(THIS_PROC_ID));
	output_debug << "PROCESSOR " << THIS_PROC_ID<< std::endl;
	std::string meshName = fileName + ".g";
	std::string assocfile = fileName + ".m2g";
	std::string geomfile = fileName + ".sat";

	percept::PerceptMesh eMesh(3u);
	eMesh.open(meshName);
	stk::mesh::MetaData * md = eMesh.get_fem_meta_data();


	PGeomACIS pg;
	pg.initialize(ACIS_GEOMETRY_ENGINE);
	bool acis_read = pg.import_acis_file(geomfile.c_str());
	EXPECT_EQ(true, acis_read);
	PGeom* pgeom = &pg;


	std::vector<int> surfIDs;
	pgeom->get_surfaces(surfIDs);
	std::vector<std::string> partNames(surfIDs.size());
	//POTENTIAL ISSUE WITH HYBRID MESHES:
	//These parts must be created beforehand. Thus, the knowledge of whether or not this particular part needs to be a SHELL_TRI_3 or SHELL_QUAD_4 needs to occur here but the
	//knowledge of whether or not this surface contains quads or tris doesn't come until AFTER the meta data's commit. If all the surfaces will be quads or tris, a boolean branching to
	//alternative (quad or tri) creation loops could be implemented
	for (unsigned i = 0; i < surfIDs.size(); i++) { //make parts that we'll store new mesh on
		std::string name = "geom_surface_";
		name = name + std::to_string(surfIDs[i]);
		stk::mesh::Part& part = md->declare_part(name,
				stk::topology::FACE_RANK);
		stk::io::put_io_part_attribute(part); //set them up so they show up on output
		partNames[i] = name;
	} //endfor




	std::vector<int> curveIDs;
	pgeom->get_curves(curveIDs);
	std::vector<std::string> curveNames(curveIDs.size());
	for (unsigned i = 0; i < curveIDs.size(); i++) { //make parts that we'll store new mesh on
		std::string name = "geom_curve_";
		name = name + std::to_string(curveIDs[i]);
		stk::mesh::Part& part = md->declare_part_with_topology(name,
				stk::topology::BEAM_2);
		stk::io::put_io_part_attribute(part); //set them up so they show up on output
		curveNames[i] = name;
	} //endfor

	//setup refinement

	stk::mesh::Part & skin_part = md->declare_part("SkinPart", stk::topology::FACE_RANK); //doing this before the break pattern should make this work

	int scalarDimension = 0;
	stk::mesh::FieldBase* proc_rank_field = eMesh.add_field("proc_rank",
			stk::topology::ELEMENT_RANK, scalarDimension);
	URP_Heterogeneous_3D break_pattern(eMesh);


	eMesh.commit();
	stk::mesh::BulkData * bd = eMesh.get_bulk_data();
	eMesh.initializeIdServer(); //initialize the ID server for creation of unique IDs for beam/shell parts

	stk::mesh::Selector localAndSharedPart = (md->locally_owned_part()); //these still contain ghosted_shared peices of mesh along boundary
			localAndSharedPart |= md->globally_shared_part();



	std::vector< stk::mesh::Part * > partToPutSidesInto(1,&skin_part);
	stk::mesh::Selector blocksToSkinOrConsider(md->locally_owned_part());
	blocksToSkinOrConsider |= (md->globally_shared_part());
	stk::mesh::create_exposed_block_boundary_sides(*eMesh.get_bulk_data(), blocksToSkinOrConsider, partToPutSidesInto);
//	stk::mesh::create_interior_block_boundary_sides(*eMesh.get_bulk_data(),blocksToSkinOrConsider, partToPutSidesInto); // comment this out to get to work




	//BEGIN ASSOCATION
	PGeomAssoc<stk::mesh::BulkData, stk::mesh::Entity, stk::mesh::Entity,
			stk::mesh::Entity, stk::mesh::Entity> geom_assoc(pgeom); //
	geom_assoc.set_node_callback(get_node_from_id);
	geom_assoc.set_edge_callback(get_beam_from_ids);
	geom_assoc.set_face_callback(get_FACE_from_ids);
	geom_assoc.set_elem_callback(get_hex_from_id);
	geom_assoc.set_validate_nodes_callback(validate_node_owner_ship);
	geom_assoc.set_mesh(eMesh.get_bulk_data());
	geom_assoc.set_fill_curve_and_surface_maps_during_import(false);
	bool read_assoc = geom_assoc.import_m2g_file(assocfile.c_str(), true);//searching for faces messes up because it grabs shared stuff
	EXPECT_EQ(true, read_assoc);
	std::cout<<"FINISHED FIRST PART OF ASSOCIATION"<<std::endl;
	output_debug<<"FINISHED FIRST PART OF ASSOCIATON"<<std::endl;


	eMesh.get_bulk_data()->modification_begin();
	/////////////////////////////////////////////
	int num_beams=0;
	output_debug <<"Num curves is: "<< curveIDs.size() <<std::endl;
	for(unsigned i = 0; i < curveIDs.size(); i++) { //create beams and put them into corresponding curve parts

		std::vector<stk::mesh::Part *> add_parts_beams(1,
				static_cast<stk::mesh::Part*>(0));

		add_parts_beams[0] = md->get_part(curveNames[i]);

		std::vector<std::vector<int>> edge_node_ids;
		geom_assoc.get_curve_edge_nodes(curveIDs[i], edge_node_ids);
//		output_debug << "got nodes in curve" << std::endl;
		for (unsigned ii = 0; ii < edge_node_ids.size(); ii++) {

			bool toDeclare = true;

			std::vector<stk::mesh::EntityId> beamNodeIDs;
			for (unsigned j = 0; j < edge_node_ids[ii].size(); j++)
				beamNodeIDs.push_back((stk::mesh::EntityId) edge_node_ids[ii][j]);
			//call function to see if you own/share the nodes only then create the beams

//			output_debug << "got all nodes in this curve" << std::endl;


			for (unsigned j = 0; j < edge_node_ids[ii].size(); j++) {




						stk::mesh::Entity cur_node = bd->get_entity(
								stk::topology::NODE_RANK, beamNodeIDs[j]); //CHECK IF NODE IS INVALID
						if ( !bd->is_valid(cur_node) ){
							toDeclare=false;
							break;
						}

//						output_debug << "got a node" << std::endl;
						stk::mesh::Bucket &node_bucket = eMesh.bucket(cur_node);
						if (!localAndSharedPart(node_bucket)) //errors (sefaults) out here on proc 0 on first iteration
							toDeclare = false;
//						output_debug <<"checked my node" <<std::endl;
					}

					if (toDeclare) {
						output_debug << "creating a beam..." << std::endl;
						stk::mesh::EntityId id2 = eMesh.getNextId(
								stk::topology::ELEMENT_RANK);
						stk::mesh::declare_element(*eMesh.get_bulk_data(),
								add_parts_beams, id2, beamNodeIDs);
						num_beams++;

						output_debug << "beam id is: "<< id2 << std::endl;

					}

				} //endfor

			} //endfor


	output_debug << "The Number of beams on this processor is: " << num_beams << "     "<< eMesh.get_number_edges() << std::endl;
	output_debug << "end beam creation" << std::endl;
	eMesh.get_bulk_data()->modification_end();
	geom_assoc.fill_curve_and_surface_maps();

	int num_faces = 0;
	eMesh.get_bulk_data()->modification_begin();


	for (unsigned i = 0; i < surfIDs.size(); i++) {

		std::vector<stk::mesh::Part*> add_parts(1,
				static_cast<stk::mesh::Part*>(0));
		std::vector<stk::mesh::Part*> remove_parts;

		add_parts[0] = md->get_part(partNames[i]);
		remove_parts.push_back(md->get_part("SkinPart"));

		std::vector<stk::mesh::Entity> faces = geom_assoc.get_surf_faces(
				surfIDs[i]);


		for (unsigned ii = 0; ii < faces.size(); ++ii) {
			stk::mesh::Entity element = faces[ii];
			if(eMesh.parallel_owner_rank(element)!=THIS_PROC_ID)
				continue;
			eMesh.get_bulk_data()->change_entity_parts(element, add_parts,
					remove_parts);
			num_faces++;

		}
		stk::mesh::fixup_ghosted_to_shared_nodes(*eMesh.get_bulk_data());

	}

	output_debug << "Number faces on proc: " << num_faces << std::endl;
	eMesh.get_bulk_data()->modification_end();



	output_debug << "geom and mesh associated" << std::endl;
//	eMesh.get_bulk_data()->modification_end();
	output_debug<<"done with modifications" << std::endl;

	eMesh.save_as("output_ManualBeamsShells_BMS_INTERMEDIATEMESH.g");//this is a valid mesh that refines when saved out


	UNIFORM_REFINER breaker(eMesh,break_pattern,proc_rank_field);
	breaker.setRemoveOldElements(true);
	std::cout << "THIS PROC NUM IS: " << THIS_PROC_ID <<std::endl;
	breaker.doBreak();
	std::cout << "done with breaker" << std::endl;

	// FOR SOME REASON THIS SNAPPING LOGIC SEVERELY DEFORMS SINGLE HEX CASE!!!!


	eMesh.save_as("output_ManualBeamsFaces_BMSR1_NOTSNAPPED.g");

		const stk::mesh::BucketVector & buckets4 =
				eMesh.get_bulk_data()->get_buckets(stk::topology::NODE_RANK,
						localAndSharedPart);

		stk::mesh::BucketVector::const_iterator k4 = buckets4.begin();

	for (; k4 != buckets4.end(); k4++) {
		std::vector<int> isInCurves;
		std::vector<int> isInSurfs;
		stk::mesh::Bucket & bucket2 = **k4;

		for (unsigned i = 0; i < curveIDs.size(); i++) {

			stk::mesh::Selector boundarySelectorInner(
					*md->get_part(curveNames[i]));

			if (boundarySelectorInner(bucket2)) //if the bucket is in our part
			{
				isInCurves.push_back(curveIDs[i]); //add it to our evalutors
			}

		}

		for (unsigned i = 0; i < surfIDs.size(); i++) {
			stk::mesh::Selector boundarySelectorInner(
					*md->get_part(partNames[i]));

			if (boundarySelectorInner(bucket2))
				isInSurfs.push_back(surfIDs[i]);
		}

		if (isInCurves.size() == 0 && isInSurfs.size() == 0)
			continue;

		const unsigned num_elements_in_bucket2 = bucket2.size();
		for (unsigned iEntity = 0; iEntity < num_elements_in_bucket2;
				iEntity++) {

			stk::mesh::Entity node = bucket2[iEntity];

			double * data = eMesh.field_data(*eMesh.get_coordinates_field(),
					node);
			// std::cout << data[0] << std::endl;
			CubitVector old_locations(data[0], data[1], data[2]);

			if (isInCurves.size() > 0) {
				double smallest_dist_sqr = std::numeric_limits<double>::max();
				for (unsigned ii = 0; ii < isInCurves.size(); ii++) {
					CubitVector new_locations = pgeom->get_curv_closest_point(
							isInCurves[ii], old_locations);

					double dist_sqr = pow(
							(old_locations.x() - new_locations.x()), 2)
							+ pow((old_locations.y() - new_locations.y()), 2)
							+ pow((old_locations.z() - new_locations.z()), 2);

					if (dist_sqr < smallest_dist_sqr) { //BEAMS DO NOT SNAP CORRECTLY
						smallest_dist_sqr = dist_sqr;

						data[0] = new_locations.x();
						data[1] = new_locations.y();
						data[2] = new_locations.z();

						output_debug << "curve Entity: " << isInCurves[ii] << " and node id is: " << eMesh.identifier(node) << "it moved :" << dist_sqr << std::endl;
						output_debug <<"     " << "new_locations: "<< "x :" << new_locations.x() << "y: " << new_locations.y() << "z: " <<new_locations.z() <<std::endl;
						output_debug <<"     " << "new_locations: "<< "x :" << new_locations.x() << "y: " << new_locations.y() << "z: " <<new_locations.z() <<std::endl;
					}

				}

			}

			else {
				double smallest_dist_sqr = std::numeric_limits<double>::max();
				for (unsigned ii = 0; ii < isInSurfs.size(); ii++) {
					CubitVector new_locations = pgeom->get_surf_closest_point(
							isInSurfs[ii], old_locations);

					double dist_sqr = pow(
							(old_locations.x() - new_locations.x()), 2)
							+ pow((old_locations.y() - new_locations.y()), 2)
							+ pow((old_locations.z() - new_locations.z()), 2);

					if (dist_sqr < smallest_dist_sqr) { //SURFACES SNAP CORRECTLY
						smallest_dist_sqr = dist_sqr;
						data[0] = new_locations.x();
						data[1] = new_locations.y();
						data[2] = new_locations.z();
					}

				}

			}

		}

	}




	std::cout<<"Refined and Snapped" << std::endl;



	eMesh.save_as("output_ManualBeamsFaces_BMSR1.g");

}//manualBeamFaces
#endif //endacis

    }//    namespace regression_tests
  }//  namespace percept


