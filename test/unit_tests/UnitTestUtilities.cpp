// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <stk_util/parallel/Parallel.hpp>

#if !STK_PERCEPT_LITE
#  if defined(STK_BUILT_IN_SIERRA)

#    include <init/Ionit_Initializer.h>

#  else

#    include <Ionit_Initializer.h>

#  endif
#endif


#include <Ioss_SubSystem.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>

#include <Shards_BasicTopologies.hpp>

#include <stk_io/IossBridge.hpp>

#include <percept/RunEnvironment.hpp>
#include <percept/PerceptMesh.hpp>

#include <iostream>
#include <assert.h>
#include <vector>

using namespace stk ;

/** \addtogroup stk_io_module
 * \{
 */

/**
 * Example code showing a basic, but complete, mesh to results output
 * coding including subsetting and periodic field input and output.
 * Includes handling of nodeblocks, element blocks, nodesets,
 * and sidesets.  Attribute fields and distribution factor
 * fields are also supported.
 *
 * This example can serve as the basis for adding binary IO support to
 * an application.  The code here uses the Ioss to/from stk::mesh
 * bridge functions in the stk::io namespace defined in IossBridge.hpp
 * include file.
 */
namespace stk_example_io {

/// Declare "coordinates" field and put it on the universal part. This
/// example also defines all Ioss::Field::TRANSIENT fields that exist on the
/// Ioss::Nodeblock as fields on the universal part.
void process_nodeblocks    (Ioss::Region &region, stk::mesh::MetaData &meta);

/// Declare a part for each element block on the Ioss::Region
/// 'region' unless the element block has the "omitted" property set
/// to the value 1. The example then iterates each element block and
/// defines any Ioss::Field::ATTRIBUTE and Ioss::Field::TRANSIENT fields that exist on the
/// Ioss::ElementBlock as fields on the corresponding part.
void process_elementblocks (Ioss::Region &region, stk::mesh::MetaData &meta);

/// Declare a part for each Ioss::NodeSet on the Ioss::Region
/// 'region' unless the nodeset has the "omitted" property set
/// to the value 1. The example then iterates each nodeset and
/// defines any "distribution factor" and Ioss::Field::TRANSIENT fields that
/// exist on the Ioss::NodeSet as fields on the corresponding
/// part.
void process_nodesets      (Ioss::Region &region, stk::mesh::MetaData &meta);

/// Declare a part for each Ioss::SideSet on the Ioss::Region
/// 'region' unless the sideset has the "omitted" property set
/// to the value 1. The example then iterates each sideset and
/// defines any "distribution factor" and Ioss::Field::TRANSIENT fields that
/// exist on the Ioss::SideSet as fields on the corresponding
/// part.
///
/// Each sideblock in the active sidesets is then processed by
/// defining a part for each Ioss::SideBlock on the Ioss::SideSet
/// unless the sideblock has the "omitted" property set to the value
/// 1. The example then iterates each sideblock and defines any
/// "distribution factor" and Ioss::Field::TRANSIENT fields that exist on the
/// Ioss::SideBlock as fields on the corresponding part.
void process_sidesets      (Ioss::Region &region, stk::mesh::MetaData &meta);

/// NOTE: This must be called after the process_elementblocks() call
/// since there may be nodes that exist in the database that are
/// not part of the analysis mesh due to subsetting of the element
/// blocks.
///
/// Populates  the "coordinates" field for all active nodes in the model.
void process_nodeblocks    (Ioss::Region &region, stk::mesh::BulkData &bulk);

/// NOTE: This should be the first function called of any of the
/// "process_X" type functions that take an stk::mesh::BulkData
/// argument, especially if the input Ioss::Region mesh is going to
/// be subsetted (have element blocks omitted).
///
/// This function iterates all non-omitted element blocks and
/// declares each element (and the corresponding nodes) in the
/// element block. If there are any Ioss::Field::ATTRIBUTE fields on the element
/// block (for example, shell thickness or particle radius), then
/// that field data is alse read and the corresponding
/// stk::mesh::Field populated.
void process_elementblocks (Ioss::Region &region, stk::mesh::BulkData &bulk);

/// Iterates each non-omitted Ioss::NodeSet and then iterates each
/// node in the Ioss::NodeSet.  If the node exists (that is, it is
/// connected to a non-omitted Ioss::ElementBlock), then that node
/// is associated with the part corresponding to this
/// Ioss::NodeSet. If the "distribution_factor" field exists, then
/// that data is also associated with the field.
void process_nodesets      (Ioss::Region &region, stk::mesh::BulkData &bulk);

/// Process each non-omitted Ioss::SideSet and the contained
/// non-omitted Ioss::SideBlock and associate each element-side pair with
/// the corresponding part if the underlying element is active.  If
/// the "distribution_factor" field exists, then that data is also
/// associated with the corresponding field.
void process_sidesets      (Ioss::Region &region, stk::mesh::BulkData &bulk);

/// A minimal example function showing how field data on the
/// Ioss::Region entities can be periodically transferred to the
/// corresponding field(s) on the stk::mesh entities. This would be
/// used to bring in initial condition data or interpolation data or
/// any other scenario in which data on the mesh file needs to be
/// transferred to the stk::mesh fields.
void process_input_request (Ioss::Region &region, stk::mesh::BulkData &bulk, int step);

/// A minimal example function showing how stk::mesh field data can
/// periodically be output to a results, history, heartbeat, or
/// restart database.  The scheduling would be done either in this
/// function or at a higher level and is not shown here. The
/// function iterates all parts and if there is a corresponding Ioss
/// part on the Ioss::Region, all fields defined to be output are
/// iterated and their data output to the corresponding
/// Ioss::Field. The function calls the
/// stk::io::is_valid_part_field() function to determine whether the
/// field should be output and then calls the
/// stk::io::field_data_to_ioss() function to do the actual output
/// of the field.
void process_output_request(Ioss::Region &region, stk::mesh::BulkData &bulk, int step);

/// This function shows the basic calls needed to perform definition
/// and input of the mesh model and definition and periodic output
/// of a results database. The function is given the mesh filename
/// and the output filename and goes through all steps of
/// associating the filename with an Ioss::DatabaseIO object of the
/// correct type ("exodusII" in this example); creating an
/// Ioss::Region and then defining an stk::mesh corresponding to
/// this mesh.  The function also provides an example of how
/// specific element blocks existing in the mesh database could be
/// omitted from the analysis model.
///
/// The example then shows how to define a results database
/// corresponding to the analysis model and periodically output the
/// results in an execute loop.
///
/// A true application would have to provide additional
/// functionality and robustness, but the example shows how the
/// basic functionality can be provided by an application.
///
/// Note that the paradigm illustrated here is different than the
/// mesh input and output paradigm provided in the current
/// framework.  In this case, the application is responsible for the
/// majority of the IO behavior and the toolkit only provides some
/// helper functions to bridge between the Ioss and the stk::mesh.
/// It is hoped that this paradigm will result in more functionality
/// for the application with less complication and overhead.

typedef mesh::Field<double>                    ScalarFieldType ;
typedef mesh::Field<double,mesh::Cartesian>    VectorFieldType ;

// Specification for the aggressive gather pointer-field for elements.

void my_test(
  mesh::BulkData & M ,
  const stk::mesh::EntityRank elem_type ,
  const VectorFieldType & coord_field ,
  const VectorFieldType & elem_centroid_field );

void io_example( stk::ParallelMachine comm,
                 const std::string& in_filename,
                 const std::string& out_filename)
{
  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  std::cout << "========================================================================\n"
            << " Use Case: Subsetting with df and attribute field input/output          \n"
            << "========================================================================\n";

  std::string dbtype("exodusII");
  Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(dbtype, in_filename, Ioss::READ_MODEL,
                                                  comm);
  if (dbi == NULL || !dbi->ok()) {
    std::cerr  << "ERROR: Could not open database '" << in_filename
               << "' of type '" << dbtype << "'\n";
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'in_region' owns 'dbi' pointer at this time...
  Ioss::Region in_region(dbi, "input_model");

  // SUBSETTING PARSING/PREPROCESSING...
  // Just an example of how application could control whether an
  // entity is subsetted or not...


  // Example command line in current code corresponding to behavior below:
  std::cout << "\nWhen processing file multi-block.g for use case 2, the blocks below will be omitted:\n";
  std::cout << "\tOMIT BLOCK Cblock Eblock I1 I2\n\n";
  Ioss::ElementBlock *eb = in_region.get_element_block("cblock");
  if (eb != NULL)
    eb->property_add(Ioss::Property(std::string("omitted"), 1));

  eb = in_region.get_element_block("eblock");
  if (eb != NULL)
    eb->property_add(Ioss::Property(std::string("omitted"), 1));

  eb = in_region.get_element_block("i1");
  if (eb != NULL)
    eb->property_add(Ioss::Property(std::string("omitted"), 1));

  eb = in_region.get_element_block("i2");
  if (eb != NULL)
    eb->property_add(Ioss::Property(std::string("omitted"), 1));

#if 0
  // Example for subsetting -- omit "odd" blocks
  if (entity->type() == Ioss::ELEMENTBLOCK) {
    int id = entity->get_property("id").get_int();
    if (id % 2) {
      entity->property_add(Ioss::Property(std::string("omitted"), 1));
      std::cout << "Skipping " << entity->type_string() << ": "  << entity->name() << "\n";
    }
  }
#endif

#define IOTEST 1

  //----------------------------------
  // Process Entity Types. Subsetting is possible.
  stk::mesh::MetaData meta_data(3);
  stk::mesh::Part & universal = meta_data.universal_part();
#if IOTEST
  process_elementblocks(in_region, meta_data);
  process_nodeblocks(in_region,    meta_data);
  process_sidesets(in_region,      meta_data);
  process_nodesets(in_region,      meta_data);
#endif

  //--------------------------------
  // Declare coordinates field on all nodes with 3D:

  VectorFieldType & coordinates_field =
    meta_data.declare_field< VectorFieldType >( stk::topology::NODE_RANK, "coordinates" );

  //int SpatialDim = 3;
  enum { SpatialDim = 3 };

  stk::mesh::put_field_on_mesh(
    coordinates_field ,  universal , SpatialDim , nullptr);

  //--------------------------------

  VectorFieldType & elem_centroid_field =
    *meta_data.get_field< VectorFieldType >( stk::topology::ELEMENT_RANK, "ind_coarse" );

  VectorFieldType & elem_centroid_field2 =
    meta_data.declare_field< VectorFieldType >( stk::topology::ELEMENT_RANK, "ind_coarse2" );

  mesh::Part * block_1 = meta_data.get_part("block_1");
  stk::mesh::put_field_on_mesh(
                       elem_centroid_field2 , *block_1 , SpatialDim , nullptr);



  //----------------------------------
  // Done populating meta data, commit and create bulk data
  meta_data.commit();

  //----------------------------------
  // Process Bulkdata for all Entity Types. Subsetting is possible.
  stk::mesh::BulkData bulk_data(meta_data, comm);
#if IOTEST
  process_elementblocks(in_region, bulk_data);
  process_nodeblocks(in_region,    bulk_data);
  process_sidesets(in_region,      bulk_data);
  process_nodesets(in_region,      bulk_data);
#endif


  //----------------------------------
  // OUTPUT...Create the output "mesh" portion

  Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, out_filename,
                                                  Ioss::WRITE_RESULTS,
                                                  comm);
  if (dbo == NULL || !dbo->ok()) {
    std::cerr << "ERROR: Could not open results database '" << out_filename
              << "' of type '" << dbtype << "'\n";
    std::exit(EXIT_FAILURE);
  }

  // NOTE: 'out_region' owns 'dbo' pointer at this time...
  Ioss::Region out_region(dbo, "results_output");
  stk::io::OutputParams params(out_region, bulk_data);
  stk::io::define_output_db(params, {}, &in_region);
  stk::io::write_output_db(params);

  // ------------------------------------------------------------------------
  /** \todo REFACTOR A real app would register a subset of the
   * fields on the mesh database as fields that the app would want
   * read at one or all or specified steps.  In this example, all
   * fields existing on the input mesh database are defined on the
   * parts in the stk::mesh.
   *
   * The real app would also only register a subset of the stk::mesh
   * fields as output fields and would probably have a mapping from
   * the internally used name to some name picked by the user. In
   * this example, all Ioss::Field::TRANSIENT fields defined on the stk::mesh are
   * output to the results database and the internal stk::mesh field
   * name is used as the name on the database....
   */

  out_region.begin_mode(Ioss::STATE_DEFINE_TRANSIENT);

  // Special processing for nodeblock (all nodes in model)...
  stk::io::ioss_add_fields(meta_data.universal_part(), stk::topology::NODE_RANK,
                           out_region.get_node_blocks()[0],
                           Ioss::Field::TRANSIENT);

  const stk::mesh::PartVector & all_parts = meta_data.get_parts();
  for ( stk::mesh::PartVector::const_iterator
          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;

    // Check whether this part should be output to results database.
    if (stk::io::is_part_io_part(*part)) {
      // Get Ioss::GroupingEntity corresponding to this part...
      Ioss::GroupingEntity *entity = out_region.get_entity(part->name());
      if (entity != NULL) {
        if (entity->type() == Ioss::SIDESET) {
	  Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
          int block_count = sset->block_count();
          for (int i=0; i < block_count; i++) {
            Ioss::EntityBlock *fb = sset->get_block(i);
            stk::io::ioss_add_fields(*part,
                                     percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                                     fb, Ioss::Field::TRANSIENT);
          }
        } else {
          stk::io::ioss_add_fields(*part,
                                   percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                                   entity, Ioss::Field::TRANSIENT);
        }
      } else {
        /// \todo IMPLEMENT handle error... Possibly an assert since
        /// I think the corresponding entity should always exist...
      }
    }
  }
  out_region.end_mode(Ioss::STATE_DEFINE_TRANSIENT);
  // ------------------------------------------------------------------------

  // Read and Write transient fields...
  out_region.begin_mode(Ioss::STATE_TRANSIENT);
  int timestep_count = in_region.get_property("state_count").get_int();
  for (int step = 1; step <= timestep_count; step++) {
    double time = in_region.get_state_time(step);

    // Read data from the io input mesh database into stk::mesh fields...
    process_input_request(in_region, bulk_data, step);

    // execute()
    my_test (bulk_data, stk::topology::ELEMENT_RANK , coordinates_field, elem_centroid_field);

    // Write data from the stk::mesh fields out to the output database.a
    int out_step = out_region.add_state(time);
    process_output_request(out_region, bulk_data, out_step);
  }
  out_region.end_mode(Ioss::STATE_TRANSIENT);
}

// ========================================================================
void process_nodeblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  assert(nb->field_exists("mesh_model_coordinates"));
  Ioss::Field coordinates = nb->get_field("mesh_model_coordinates");
  int spatial_dim = coordinates.transformed_storage()->component_count();

  stk::mesh::Field<double,stk::mesh::Cartesian> & coord_field =
    meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");

  stk::mesh::put_field_on_mesh( coord_field, meta.universal_part(),
                        spatial_dim, nullptr);

  /** \todo IMPLEMENT truly handle fields... For this case we are
   * just defining a field for each transient field that is present
   * in the mesh...
   */
  stk::io::define_io_fields(nb, Ioss::Field::TRANSIENT, meta.universal_part(),stk::topology::NODE_RANK);
}

// ========================================================================
void process_elementblocks(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();
  mesh::EntityRank stk_mesh_Element = stk::topology::ELEMENT_RANK;
  stk::io::default_part_processing(elem_blocks, meta, stk_mesh_Element);

  // Parts were created above, now handle element block specific
  // information (topology, attributes, ...);
  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());
      assert(part != NULL);

      // Element Block attributes (if any)...
      /** \todo IMPLEMENT truly handle attribute fields... For this
       * case we are just defining a field for each attribute field
       * that is present in the mesh...
       */
      stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE,
                                *part,
                                percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ) );

      /** \todo IMPLEMENT truly handle fields... For this case we
       * are just defining a field for each transient field that is
       * present in the mesh...
       */
      stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                *part,
                                percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ) );

      const CellTopologyData* cell_topo = meta.get_cell_topology(*part).getCellTopologyData();
      std::string cell_topo_name = "UNKNOWN";
      if (cell_topo != NULL)
        cell_topo_name = cell_topo->name;

      std::cout << entity->type_string() << ": " << entity->name()
                << " , celltop = " << cell_topo_name
                << std::endl ;
    }
  }
}

// ========================================================================
void process_nodesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();
  stk::io::default_part_processing(node_sets, meta, stk::topology::NODE_RANK);

  /** \todo REFACTOR should "distribution_factor" be a default field
   * that is automatically declared on all objects that it exists
   * on as is done in current framework?
   */
  stk::mesh::Field<double> & distribution_factors_field =
    meta.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "distribution_factors");

  /** \todo REFACTOR How to associate distribution_factors field
   * with the nodeset part if a node is a member of multiple
   * nodesets
   */

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      stk::mesh::Part* const part = meta.get_part(entity->name());
      assert(part != NULL);
      assert(entity->field_exists("distribution_factors"));

      stk::mesh::put_field_on_mesh(distribution_factors_field, *part, nullptr);

      /** \todo IMPLEMENT truly handle fields... For this case we
       * are just defining a field for each transient field that is
       * present in the mesh...
       */
      stk::io::define_io_fields(entity, Ioss::Field::TRANSIENT,
                                *part,
                                percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ) );
    }
  }
}

// ========================================================================
void process_surface_entity(Ioss::SideSet *entity, stk::mesh::MetaData &meta,
                            stk::mesh::EntityRank entity_rank)
{
  assert(entity->type() == Ioss::SIDESET);
  const Ioss::SideBlockContainer& blocks = entity->get_side_blocks();
  stk::io::default_part_processing(blocks, meta, entity_rank);

  stk::mesh::Part* const fs_part = meta.get_part(entity->name());
  assert(fs_part != NULL);

  stk::mesh::Field<double, stk::mesh::ElementNode> *distribution_factors_field = NULL;
  bool surface_df_defined = false; // Has the surface df field been defined yet?

  int block_count = entity->block_count();
  for (int i=0; i < block_count; i++) {
    Ioss::SideBlock *fb = entity->get_block(i);
    if (stk::io::include_entity(fb)) {
      std::cout << fb->type_string() << " " << fb->name() << "\n";
      stk::mesh::Part * const fb_part = meta.get_part(fb->name());
      assert(fb_part != NULL);
      meta.declare_part_subset(*fs_part, *fb_part);

      stk::topology::rank_t rank = static_cast<stk::topology::rank_t>(percept::PerceptMesh::fem_entity_rank( fb_part->primary_entity_rank() ));

      if (fb->field_exists("distribution_factors")) {
        if (!surface_df_defined) {
          std::string field_name = entity->name() + "_df";
          distribution_factors_field =
          &meta.declare_field<stk::mesh::Field<double, stk::mesh::ElementNode> >(rank, field_name);
          stk::io::set_distribution_factor_field(*fs_part, *distribution_factors_field);
          surface_df_defined = true;
        }
        stk::io::set_distribution_factor_field(*fb_part, *distribution_factors_field);
        int face_node_count = fb->topology()->number_nodes();
        stk::mesh::put_field_on_mesh(*distribution_factors_field,
                             *fb_part, face_node_count, nullptr);
      }

      /** \todo IMPLEMENT truly handle fields... For this case we
       * are just defining a field for each transient field that is
       * present in the mesh...
       */
      stk::io::define_io_fields(fb, Ioss::Field::TRANSIENT,
                                *fb_part,
                                percept::PerceptMesh::fem_entity_rank( fb_part->primary_entity_rank() ) );
    }
  }
}

// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::MetaData &meta)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();
  //!<
  mesh::EntityRank stk_mesh_Side = stk::topology::FACE_RANK;
  stk::io::default_part_processing(side_sets, meta, stk_mesh_Side);

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, meta, stk_mesh_Side);  // FIXME
    }
  }
}

// ========================================================================
// Bulk Data
// ========================================================================
void process_nodeblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // This must be called after the "process_element_blocks" call
  // since there may be nodes that exist in the database that are
  // not part of the analysis mesh due to subsetting of the element
  // blocks.

  const Ioss::NodeBlockContainer& node_blocks = region.get_node_blocks();
  assert(node_blocks.size() == 1);

  Ioss::NodeBlock *nb = node_blocks[0];

  std::vector<stk::mesh::Entity> nodes;
  stk::io::get_input_entity_list(nb, stk::topology::NODE_RANK, bulk, nodes);

  /** \todo REFACTOR Application would probably store this field
   * (and others) somewhere after the declaration instead of
   * looking it up each time it is needed.
   */
  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
  stk::mesh::Field<double,stk::mesh::Cartesian> *coord_field =
    meta.get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");

  stk::io::field_data_from_ioss(bulk, coord_field, nodes, nb, "mesh_model_coordinates");
}

// ========================================================================
void process_elementblocks(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const Ioss::ElementBlockContainer& elem_blocks = region.get_element_blocks();

  for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
      it != elem_blocks.end(); ++it) {
    Ioss::ElementBlock *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string &name = entity->name();
      const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);

      const CellTopologyData* cell_topo = stk::mesh::MetaData::get(bulk).get_cell_topology(*part).getCellTopologyData();
      if (cell_topo == NULL) {
        std::ostringstream msg ;
        msg << " UnitTestUtilities::process_elementblocks::INTERNAL_ERROR: Part " << part->name() << " returned NULL from m_eMesh.get_cell_topology()";
        throw std::runtime_error( msg.str() );
      }

      std::vector<int> elem_ids ;
      std::vector<int> connectivity ;

      entity->get_field_data("ids", elem_ids);
      entity->get_field_data("connectivity", connectivity);

      size_t element_count = elem_ids.size();
      int nodes_per_elem = cell_topo->node_count ;
      stk::mesh::EntityIdVector e_connectivity(nodes_per_elem) ;

      std::vector<stk::mesh::Entity> elements(element_count);
      for(size_t i=0; i<element_count; ++i) {
        int *conn = &(connectivity[i*nodes_per_elem]);
        for (int j=0; j < nodes_per_elem; j++)
          {
            e_connectivity[j] = (mesh::EntityId)conn[j];
          }
        elements[i] = stk::mesh::declare_element(bulk, *part, elem_ids[i], e_connectivity);
      }

      // For this example, we are just taking all attribute fields
      // found on the io database and populating fields on the
      // corresponding mesh part.  In practice, would probably be
      // selective about which attributes to use...
      Ioss::NameList names;
      entity->field_describe(Ioss::Field::ATTRIBUTE, &names);
      for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
        if (*I == "attribute" && names.size() > 1)
          continue;
        stk::mesh::FieldBase *field = stk::mesh::get_field_by_name(*I,meta);
        stk::io::field_data_from_ioss(bulk, field, elements, entity, *I);

      }
    }
  }
}

// ========================================================================
void process_nodesets(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  // Should only process nodes that have already been defined via the element
  // blocks connectivity lists.
  const Ioss::NodeSetContainer& node_sets = region.get_nodesets();

  for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
      it != node_sets.end(); ++it) {
    Ioss::NodeSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      const std::string & name = entity->name();
      const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);
      stk::mesh::Part* const part = meta.get_part(name);
      assert(part != NULL);
      stk::mesh::PartVector add_parts( 1 , part );

      std::vector<int> node_ids ;
      int node_count = entity->get_field_data("ids", node_ids);

      std::vector<stk::mesh::Entity> nodes(node_count);
      for(int i=0; i<node_count; ++i) {
        nodes[i] = bulk.get_entity( stk::topology::NODE_RANK, node_ids[i] );
        if (bulk.is_valid(nodes[i]))
          bulk.declare_node(node_ids[i], add_parts );
      }


      /** \todo REFACTOR Application would probably store this field
       * (and others) somewhere after the declaration instead of
       * looking it up each time it is needed.
       */
      stk::mesh::Field<double> *df_field =
        meta.get_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, "distribution_factors");

      if (df_field != NULL) {
        stk::io::field_data_from_ioss(bulk, df_field, nodes, entity, "distribution_factors");
      }
    }
  }
}

// ========================================================================
void process_surface_entity(const Ioss::SideSet* io ,
                            stk::mesh::BulkData & bulk)
{
  assert(io->type() == Ioss::SIDESET);
  const stk::mesh::MetaData& meta = stk::mesh::MetaData::get(bulk);

  int block_count = io->block_count();
  for (int i=0; i < block_count; i++) {
    Ioss::SideBlock *block = io->get_block(i);
    if (stk::io::include_entity(block)) {
      std::vector<int> elem_side ;

      stk::mesh::Part * const fb_part = meta.get_part(block->name());

      block->get_field_data("element_side", elem_side);

      stk::mesh::PartVector add_parts( 1 , fb_part );

      size_t side_count = elem_side.size() / 2;
      std::vector<stk::mesh::Entity> sides(side_count);
      for(size_t is=0; is<side_count; ++is) {

        mesh::EntityRank stk_mesh_Element = stk::topology::ELEMENT_RANK;

        stk::mesh::Entity const elem = bulk.get_entity(stk_mesh_Element, elem_side[is*2]);
        // If NULL, then the element was probably assigned to an
        // Ioss uses 1-based side ordinal, stk::mesh uses 0-based.
        // Hence the '-1' in the following line.
        int side_ordinal = elem_side[is*2+1] - 1 ;

        // element block that appears in the database, but was
        // subsetted out of the analysis mesh. Only process if
        // non-null.
        if (bulk.is_valid(elem)) {
          stk::mesh::Entity side = bulk.declare_element_side(elem, side_ordinal, stk::mesh::ConstPartVector{});
          bulk.change_entity_parts( side, add_parts );
          sides[is] = side;
        } else {
          sides[is] = stk::mesh::Entity();
        }
      }

      const stk::mesh::FieldBase *df_field = stk::io::get_distribution_factor_field(*fb_part);
      if (df_field != NULL) {
        stk::io::field_data_from_ioss(bulk, df_field, sides, block, "distribution_factors");
      }
    }
  }
}

// ========================================================================
void process_sidesets(Ioss::Region &region, stk::mesh::BulkData &bulk)
{
  const Ioss::SideSetContainer& side_sets = region.get_sidesets();

  for(Ioss::SideSetContainer::const_iterator it = side_sets.begin();
      it != side_sets.end(); ++it) {
    Ioss::SideSet *entity = *it;

    if (stk::io::include_entity(entity)) {
      process_surface_entity(entity, bulk);
    }
  }
}

// ========================================================================
// ========================================================================
void get_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    Ioss::Field::RoleType filter_role)
{
  std::vector<stk::mesh::Entity> entities;
  stk::io::get_input_entity_list(io_entity, part_type, bulk, entities);

  stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
  const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

  std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const stk::mesh::FieldBase *f = *I; ++I;
    if (stk::io::is_valid_part_field(f, part_type, part, filter_role)) {
      stk::io::field_data_from_ioss(bulk, f, entities, io_entity, f->name());
    }
  }
}

void process_input_request(Ioss::Region &region,
                           stk::mesh::BulkData &bulk,
                           int step)
{
  region.begin_state(step);

  // Special processing for nodeblock (all nodes in model)...
  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(bulk);

  // ??? Get field data from nodeblock...
  get_field_data(bulk, meta.universal_part(), stk::topology::NODE_RANK,
                 region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

  const stk::mesh::PartVector & all_parts = meta.get_parts();
  for ( stk::mesh::PartVector::const_iterator
          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;

    // Check whether this part should be output to results database.
    if (stk::io::is_part_io_part(*part)) {
      // Get Ioss::GroupingEntity corresponding to this part...
      Ioss::GroupingEntity *entity = region.get_entity(part->name());
      if (entity != NULL) {
        if (entity->type() == Ioss::SIDESET) {
	  Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
          int block_count = sset->block_count();
          for (int i=0; i < block_count; i++) {
            Ioss::SideBlock *fb = sset->get_block(i);
            get_field_data(bulk, *part,
                           percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                           fb, Ioss::Field::TRANSIENT);
          }
        } else {
          get_field_data(bulk, *part,
                         percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                         entity, Ioss::Field::TRANSIENT);
        }
      } else {
        /// \todo IMPLEMENT handle error... Possibly an assert since
        /// I think the corresponding entity should always exist...
      }
    }
  }

  region.end_state(step);
}

void put_field_data(stk::mesh::BulkData &bulk, stk::mesh::Part &part,
                    stk::mesh::EntityRank part_type,
                    Ioss::GroupingEntity *io_entity,
                    Ioss::Field::RoleType filter_role)
{
  std::vector<stk::mesh::Entity> entities;
  stk::io::OutputParams params(bulk);
  stk::io::get_output_entity_list(io_entity, part_type, params, entities);

  stk::mesh::MetaData & meta = stk::mesh::MetaData::get(part);
  const std::vector<stk::mesh::FieldBase*> &fields = meta.get_fields();

  std::vector<stk::mesh::FieldBase *>::const_iterator I = fields.begin();
  while (I != fields.end()) {
    const stk::mesh::FieldBase *f = *I; ++I;
    if (stk::io::is_valid_part_field(f, part_type, part, filter_role)) {
      stk::io::field_data_to_ioss(bulk, f, entities, io_entity, f->name(), filter_role);
    }
  }
}

void process_output_request(Ioss::Region &region,
                            stk::mesh::BulkData &bulk,
                            int step)
{
  region.begin_state(step);
  // Special processing for nodeblock (all nodes in model)...
  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(bulk);

  put_field_data(bulk, meta.universal_part(), stk::topology::NODE_RANK,
                 region.get_node_blocks()[0], Ioss::Field::TRANSIENT);

  const stk::mesh::PartVector & all_parts = meta.get_parts();
  for ( stk::mesh::PartVector::const_iterator
          ip = all_parts.begin(); ip != all_parts.end(); ++ip ) {

    stk::mesh::Part * const part = *ip;

    // Check whether this part should be output to results database.
    if (stk::io::is_part_io_part(*part)) {

      // Get Ioss::GroupingEntity corresponding to this part...
      Ioss::GroupingEntity *entity = region.get_entity(part->name());
      if (entity != NULL) {

        if (entity->type() == Ioss::SIDESET) {
	  Ioss::SideSet *sset = dynamic_cast<Ioss::SideSet*>(entity);
          int block_count = sset->block_count();

          for (int i=0; i < block_count; i++) {
            Ioss::SideBlock *fb = sset->get_block(i);
            put_field_data(bulk, *part,
                           percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                           fb, Ioss::Field::TRANSIENT);
          }
        } else {
          put_field_data(bulk, *part,
                         percept::PerceptMesh::fem_entity_rank( part->primary_entity_rank() ),
                         entity, Ioss::Field::TRANSIENT);
        }
      } else {
        /// \todo IMPLEMENT handle error... Possibly an assert since
        /// I think the corresponding entity should always exist...
      }
    }
  }
  region.end_state(step);
}

void my_test(
  mesh::BulkData & M ,
  const stk::mesh::EntityRank elem_type ,
  const VectorFieldType & coord_field ,
  const VectorFieldType & elem_centroid_field )
{
  const mesh::MetaData & meta_data = stk::mesh::MetaData::get(M);

  // Get vector of buckets ( entities and field data)
  // for which the sides are all locally owned.

  mesh::Selector select_owned( meta_data.locally_owned_part() );

  const std::vector<mesh::Bucket*> & buckets = M.buckets( elem_type );
  //double sum[3];

  for ( std::vector<mesh::Bucket *>::const_iterator
          ik = buckets.begin() ; ik != buckets.end() ; ++ik ) if ( select_owned( **ik ) ) {

      const mesh::Bucket & bucket = **ik ;

      // Number of elems in this bucket of elems and elem field data

      const int number = bucket.size();
      ThrowAssert(number > 0);

      //double * elem_node_data = field_data( coord_field , bucket.begin() );
      double * elem_centroid_data = stk::mesh::field_data( elem_centroid_field , *bucket.begin() );

      for ( int i = 0 ; i < number ; ++i, elem_centroid_data += 3) {

        mesh::Entity elem = bucket[i] ;

        unsigned num_nodes = M.num_connectivity(elem, stk::topology::NODE_RANK);
        const stk::mesh::Entity *elem_nodes = M.begin(elem, stk::topology::NODE_RANK);

        if ( num_nodes == 8 ) {

          for (int ic=0; ic < 3; ic++) elem_centroid_data[ic]=0;

          for (int inode=0; inode < 8; inode++)
          {
            mesh::Entity node = elem_nodes[inode];

            double * const node_data = stk::mesh::field_data(coord_field , node );
            for (int ic=0; ic < 3; ic++) {
              elem_centroid_data[ic] += node_data[ic]/8.;
            }
            //double * const elem2_data = field_data( elem_centroid , elem2 );

            // Which way is the side oriented, natural for #1 or #2 ?
            //node_data[0] = 1;
            //elem_centroid_data[ic] = sum[ic];
          }
        }
      }
    }
}

}

// ========================================================================

#include <stk_util/parallel/BroadcastArg.hpp>

int myMain(int argc, char** argv)
{
  //----------------------------------
  // Broadcast argc and argv to all processors.

  stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

  stk::BroadcastArg b_arg(comm, argc, argv);

  //----------------------------------
  // Process the broadcast command line arguments

  percept::RunEnvironment run_environment(&argc, &argv);

  run_environment.clp.setDocString("options");

  std::string mesh = "";
  run_environment.clp.setOption("mesh",         &mesh, "mesh file" );
  percept::processCommandLine(run_environment.clp, argc, argv);

  //----------------------------------

  if ( mesh.length() ) {
    std::string in_filename = mesh;
    std::string out_filename = in_filename + ".out";
    std::cout << "file: " << in_filename << "\n";
    stk_example_io::io_example(comm, in_filename, out_filename );
  } else {
    std::cout << "OPTION ERROR: The '--mesh <filename>' option is required!\n";
    std::exit(EXIT_FAILURE);
  }
  stk::parallel_machine_finalize();

  return 0;
}

/**
 * \}
 */

