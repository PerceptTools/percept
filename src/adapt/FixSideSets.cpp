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

#include <percept/Percept.hpp>
#include <percept/PerceptUtils.hpp>

#include <adapt/Refiner.hpp>
#include <adapt/RefinerUtil.hpp>
#include <adapt/FixSideSets.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/MeshUtil.hpp>
#include <adapt/AdaptedMeshVerifier.hpp>

#include <stk_util/environment/memory_util.hpp>
#include <stk_util/util/human_bytes.hpp>

#define DEBUG_GSPR 0
//#define LTRACE (m_eMesh.get_rank() == 0)
#define LTRACE 0


namespace percept {


  static bool s_debugCSF = false;
#define S_ONLYLEAFTOLEAF_DEFAULT true
  bool s_onlyLeafToLeaf = S_ONLYLEAFTOLEAF_DEFAULT;
  bool s_allow_duplicates = false;


    struct LocalTimer {
      const std::string& m_name;
      Refiner *m_refiner;
      stk::diag::Timer *m_timer;
      LocalTimer(const std::string& name, Refiner *ref) : m_name(name), m_refiner(ref), m_timer(0)
      {
#if 1
        if (m_refiner)
          {
            m_timer = new stk::diag::Timer(m_name, m_refiner->rootTimer());
            m_timer->start();
          }
#endif
      }
      ~LocalTimer() {
        if (m_timer)
          {
#if 1
            m_timer->stop();
            delete m_timer;
#endif
          }
      }
    };

    void FixSideSets::fix_permutation(SetOfEntities& side_set)
    {
      LocalTimer timer_bss("FSS0_fix_perm", m_refiner);

      std::vector<stk::mesh::Entity> side_nodes;
      std::vector<stk::mesh::ConnectivityOrdinal> nords;

      for (auto side : side_set)
        {
          //stk::mesh::Entity side = *it_side;
          stk::mesh::EntityId cid = 0; //637; //15608; //3650;
          bool ldb = m_eMesh.id(side) == cid && m_eMesh.entity_rank(side) == m_eMesh.side_rank();

          if (m_eMesh.aura(side))
            continue;

          if (m_eMesh.has_default_perm(side))
            continue;

          // skip edges in 3D
          if (m_eMesh.get_spatial_dim() == 3 && m_eMesh.entity_rank(side) == m_eMesh.edge_rank())
            continue;

          percept::MyPairIterRelation side_elements (m_eMesh, side, m_eMesh.element_rank());

          stk::topology side_topo = m_eMesh.topology(side);
          bool found = false;
          stk::mesh::Entity element = stk::mesh::Entity();
          unsigned side_ord = 0;
          stk::mesh::EntityId minId = std::numeric_limits<stk::mesh::EntityId>::max();
          std::ostringstream str;
          if (!m_avoidFixSideSetChecks) VERIFY_OP_ON(side_elements.size(), >, 0, "no connected elements to this side");
          bool hasPosPerm = false;
          for (unsigned ii = 0; ii < side_elements.size(); ++ii)
            {
              stk::mesh::Entity elem = side_elements[ii].entity();
              bool sameOwner = m_eMesh.owner_rank(elem) == m_eMesh.owner_rank(side);
              stk::mesh::Permutation perm = m_eMesh.find_permutation(elem, side, side_elements[ii].relation_ordinal());
              if (perm < side_topo.num_positive_permutations())
                hasPosPerm = true;
              if (ldb)
                {
                  bool isPos = m_eMesh.is_positive_perm(elem, side, side_elements[ii].relation_ordinal());
                  str << "P[" << m_eMesh.get_rank() << "] FP TMP srk elem= " << m_eMesh.id(elem) << " side_elements.size= " << side_elements.size()
                      << " sameOwner= " << sameOwner << " isPos= " << isPos
                      << std::endl;
                  m_eMesh.print(str,side,true,true);
                  m_eMesh.print(str,elem,true,true);
                }

              if (sameOwner)
                {
                  if (m_eMesh.id(elem) < minId)
                    {
                      minId = m_eMesh.id(elem);
                      element = elem;
                      side_ord = side_elements[ii].relation_ordinal();
                    }
                  found = true;
                }
            }
          if (hasPosPerm)
            continue;

          if (!found && m_avoidFixSideSetChecks)
            continue;
          if (!found && !m_avoidFixSideSetChecks)
            {
              str << "\n" << m_eMesh.rank() << " FP not found side= " << m_eMesh.id(side)
                        << std::endl;
              m_eMesh.print(str, side, true, true);
              std::cerr << str.str()  << std::endl;
              VERIFY_MSG("couldn't find an element to reattach on same proc, side_elements.size= "+toString(side_elements.size()));
            }

          // reorient to use permIndex = 0 always
          stk::topology element_topo = m_eMesh.topology(element);
          bool isShell = element_topo.is_shell();
          const stk::mesh::Entity *element_nodes = m_eMesh.get_bulk_data()->begin_nodes(element);
          stk::mesh::Entity expected_nodes[100];
          switch (side_topo.rank())
            {
            case stk::topology::EDGE_RANK:
              element_topo.edge_nodes(element_nodes, side_ord, expected_nodes);
              break;
            case stk::topology::FACE_RANK:
              element_topo.face_nodes(element_nodes, side_ord, expected_nodes);
              break;
            default:
              VERIFY_MSG("bad side rank");
            }

          if (ldb)
            {
              str << "P[" << m_eMesh.get_rank() << "] TMP srk isShell= " << isShell << " element= " << m_eMesh.id(element)
                        << " element_topo.name= " << element_topo
                        << " erank = " << m_eMesh.entity_rank(element)
                //<< " perm= " << perm
                        << std::endl;
              m_eMesh.print(str,side,true,true);
              m_eMesh.print(str,element,true,true);
            }

          bool ldeb = ldb;
          side_nodes.assign(m_eMesh.get_bulk_data()->begin_nodes(side), m_eMesh.get_bulk_data()->end_nodes(side));
          nords.assign( m_eMesh.get_bulk_data()->begin_node_ordinals(side),  m_eMesh.get_bulk_data()->end_node_ordinals(side));
          VERIFY_OP_ON(side_nodes.size(), ==, side_topo.num_nodes(), "bad side_nodes size");

          for (unsigned jj=0; jj < side_topo.num_nodes(); ++jj)
            {
              if (ldeb) str << " node0= " << m_eMesh.get_bulk_data()->entity_key(side_nodes[jj]) << std::endl;

              bool del = m_eMesh.get_bulk_data()->destroy_relation( side, side_nodes[jj], nords[jj]);

              VERIFY_OP_ON(del, ==, true, "fix_side_sets_2:: destroy_relation failed 4");
            }

          for (unsigned ir = 0; ir < side_topo.num_nodes(); ++ir)
            {
              if (ldeb) str << " node1= " << m_eMesh.get_bulk_data()->entity_key(expected_nodes[ir]) << std::endl;
              m_eMesh.get_bulk_data()->declare_relation(side, expected_nodes[ir], ir);
            }
          if (ldb)
            {
              stk::mesh::Permutation perm = stk::mesh::INVALID_PERMUTATION;
              bool is_bad_new = m_eMesh.is_perm_bad(element, side, side_ord, perm);
              VERIFY_OP_ON(is_bad_new, ==, false, "bad node connections after renumber");
              std::cout << str.str() << std::endl;
            }
        }
    }

    bool FixSideSets::connect(stk::mesh::Entity side, bool& valid_side_part_map, SetOfEntities* avoid_elems, bool onlyPosPerm)
    {
      const stk::mesh::Entity * side_nodes = m_eMesh.get_bulk_data()->begin(side, m_eMesh.node_rank());
      const unsigned side_nodes_size = m_eMesh.get_bulk_data()->num_nodes(side);
      bool found = false;
      stk::topology side_topo = m_eMesh.topology(side);

#if 0
      SetOfEntities elem_set;
      for (unsigned isnode=0; isnode < side_nodes_size; isnode++)
        {
          //percept::MyPairIterRelation node_elements ( m_eMesh, side_nodes[isnode], m_eMesh.element_rank());
          const stk::mesh::Entity * node_elements = m_eMesh.get_bulk_data()->begin_elements( side_nodes[isnode] );
          const unsigned node_elements_size = m_eMesh.get_bulk_data()->num_elements( side_nodes[isnode] );
          for (unsigned ienode=0; ienode < node_elements_size; ienode++)
            {
              stk::mesh::Entity element = node_elements[ienode];

              if (avoid_elems && avoid_elems->find(element) != avoid_elems->end())
                continue;

              if (m_eMesh.get_bulk_data()->num_connectivity(element, m_eMesh.node_rank()) == 0)
                continue;

              if (!m_eMesh.owned(element))
                continue;

              stk::topology elem_topo = m_eMesh.topology(element);
              if (m_eMesh.get_spatial_dim() == 3 && m_eMesh.entity_rank(side) == m_eMesh.edge_rank())
                {
                  stk::topology elem_edge_topo = elem_topo.edge_topology();
                  if (side_topo != elem_edge_topo)
                    {
                      continue;
                    }
                }
              else
                {
                  stk::topology elem_side_topo = elem_topo.side_topology();
                  if (side_topo != elem_side_topo)
                    {
                      // check for wedge or pyramid
                      if (elem_topo.has_homogeneous_faces())
                        continue;
                      bool fnd = false;
                      for (unsigned ii=0; ii < elem_topo.num_sides(); ++ii)
                        {
                          if (elem_topo.side_topology(ii) == side_topo)
                            {
                              fnd = true;
                              break;
                            }
                        }
                      if (!fnd)
                        continue;
                    }
                }
              elem_set.insert(element);
            }
        }
#else
      std::vector<SetOfEntities> elem_set_nodal(side_nodes_size);
      SetOfEntities& elem_set = elem_set_nodal[side_nodes_size - 1];

      for (unsigned isnode=0; isnode < side_nodes_size; isnode++)
        {
          //percept::MyPairIterRelation node_elements ( m_eMesh, side_nodes[isnode], m_eMesh.element_rank());
          const stk::mesh::Entity * node_elements = m_eMesh.get_bulk_data()->begin_elements( side_nodes[isnode] );
          const unsigned node_elements_size = m_eMesh.get_bulk_data()->num_elements( side_nodes[isnode] );
          for (unsigned ienode=0; ienode < node_elements_size; ienode++)
            {
              stk::mesh::Entity element = node_elements[ienode];

              if (avoid_elems && avoid_elems->find(element) != avoid_elems->end())
                continue;

              if (m_eMesh.get_bulk_data()->num_connectivity(element, m_eMesh.node_rank()) == 0)
                continue;

              if (!m_eMesh.owned(element))
                continue;

              stk::topology elem_topo = m_eMesh.topology(element);
              if (m_eMesh.get_spatial_dim() == 3 && m_eMesh.entity_rank(side) == m_eMesh.edge_rank())
                {
                  stk::topology elem_edge_topo = elem_topo.edge_topology();
                  if (side_topo != elem_edge_topo)
                    {
                      continue;
                    }
                }
              else
                {
                  stk::topology elem_side_topo = elem_topo.side_topology();
                  if (side_topo != elem_side_topo)
                    {
                      // check for wedge or pyramid
                      if (elem_topo.has_homogeneous_faces())
                        continue;
                      bool fnd = false;
                      for (unsigned ii=0; ii < elem_topo.num_sides(); ++ii)
                        {
                          if (elem_topo.side_topology(ii) == side_topo)
                            {
                              fnd = true;
                              break;
                            }
                        }
                      if (!fnd)
                        continue;
                    }
                }
              // put all elements in first node's set
              if (isnode == 0)
                {
                  elem_set_nodal[0].insert(element);
                }
              else
                {
                  // only put element in current set if it is also in the previous node's set
                  if (elem_set_nodal[isnode-1].find(element) != elem_set_nodal[isnode-1].end())
                    {
                      elem_set_nodal[isnode].insert(element);
                    }
                }
            }
        }
#endif

      static std::vector<stk::mesh::Entity> elems_to_connect_to;
      static std::vector<stk::mesh::ConnectivityOrdinal> ordinals_to_connect_to;
      elems_to_connect_to.resize(0);
      ordinals_to_connect_to.resize(0);

      for (auto eit : elem_set)
        {
          stk::mesh::Entity element = eit;

          if (0)
            {
              int permIndex = 0, permPolarity = 0, k_side=0;
              bool sc = m_eMesh.should_connect(side, element, &permIndex, &permPolarity, &k_side);
              VERIFY_OP_ON(sc, ==, true, "should_connect");
            }

          if (s_debugCSF)
            {
              bool elementIsLeaf = m_eMesh.isLeafElement(element);
              bool sideIsLeaf = m_eMesh.isLeafElement(side);
              VERIFY_OP_ON(sideIsLeaf, ==, (m_eMesh.numChildren(side) == 0), "bad sideIsLeaf");
              VERIFY_OP_ON(elementIsLeaf, ==, (m_eMesh.numChildren(element) == 0), "bad elementIsLeaf");

              if (m_eMesh.should_connect(side, element))
                {

                  std::cout << "P[" << m_eMesh.get_rank() << "] side/element pair found: {" << m_eMesh.id(side) << ", " << m_eMesh.id(element) << "}"
                            << " sideIsLeaf= " << sideIsLeaf << " elementIsLeaf= " << elementIsLeaf
                            << " parent(side) = " << m_eMesh.printParent(side, true)
                            << " pv= " << m_eMesh.print_part_vector_string(m_eMesh.bucket(side).supersets(), "\n")
                            << std::endl;
                }
            }
          std::pair<bool,bool> didConnect(false,false);
          stk::mesh::ConnectivityOrdinal k_element_side;
          if (s_onlyLeafToLeaf)
            {
              //if ((sideIsLeaf && elementIsLeaf) || (!sideIsLeaf && !elementIsLeaf))
              bool elementIsLeaf = m_eMesh.isLeafElement(element);
              bool sideIsLeaf = m_eMesh.isLeafElement(side);
              if ((sideIsLeaf && elementIsLeaf) || !sideIsLeaf)
                {
                  didConnect = connectSidesForced(element, side, valid_side_part_map, &k_element_side, onlyPosPerm);
                }
            }
          else
            {
              // allow any element to connect to this side
              didConnect = connectSidesForced(element, side, valid_side_part_map, &k_element_side, onlyPosPerm);
            }

          if (didConnect.first)
            {
              elems_to_connect_to.push_back(element);
              ordinals_to_connect_to.push_back(k_element_side);
              if (!didConnect.second)
                {
                  found = true;
                }
            }
        }

      if (found)
        {
          disconnect_entity(side);
          for (unsigned irel=0; irel < elems_to_connect_to.size(); ++irel)
            {
              m_eMesh.get_bulk_data()->declare_relation(elems_to_connect_to[irel], side, ordinals_to_connect_to[irel]);
            }
        }

      return found;
    }


    // if the element (element) has a side that matches  the given side (side), connect them but first delete old connections
    std::pair<bool,bool> FixSideSets::connectSidesForced(stk::mesh::Entity element, stk::mesh::Entity side, bool& valid_side_part_map,
                                         stk::mesh::ConnectivityOrdinal *k_element_side_ret, bool onlyPosPerm )
    {
      //EXCEPTWATCH;
      bool use_coordinate_compare=false;

      valid_side_part_map = true;
      stk::mesh::EntityId cid = 0;
      s_debugCSF = m_eMesh.id(side) == cid;
      bool debug = s_debugCSF;

      // TODO
      shards::CellTopology element_topo(m_eMesh.get_cell_topology(element));
      unsigned element_nsides = (unsigned)element_topo.getSideCount();

      if (debug) {
        std::cout << "\n\ntmp srk connectSidesForced element= "; m_eMesh.print(element, true, true);
        std::cout << " side= "; m_eMesh.print(side, true, true);

        stk::mesh::PartVector const& elem_parts = m_eMesh.bucket(element).supersets();
        if (0)
          {
            for (unsigned isp = 0; isp < elem_parts.size(); isp++)
              {
                std::cout << "elem parts= " << elem_parts[isp]->name() << std::endl;
              }
          }
        std::cout << std::endl;
      }

      // check validity of connection
      // TODO
      if (m_side_part_map.size())
        {
          //std::cout << "m_side_part_map.size= " << m_side_part_map.size() << std::endl;
          bool valid = false;
          stk::mesh::PartVector const& elem_parts = m_eMesh.bucket(element).supersets();
          stk::mesh::PartVector const& side_parts = m_eMesh.bucket(side).supersets();
          for (unsigned isp = 0; isp < side_parts.size(); isp++)
            {
              if ( stk::mesh::is_auto_declared_part(*side_parts[isp]) )
                continue;

              const AutoPart *side_auto_part = side_parts[isp]->attribute<AutoPart>();
              if (side_auto_part)
                continue;

              SidePartMap::iterator found = m_side_part_map.find(side_parts[isp]);
              if (found == m_side_part_map.end())
                {
                  std::cout << "side_part = " << side_parts[isp]->name() << std::endl;
                  throw std::runtime_error("FixSideSets::connectSidesForced: couldn't find side map part");
                  //continue;
                }
              for (unsigned iep = 0; iep < elem_parts.size(); iep++)
                {
                  if ( stk::mesh::is_auto_declared_part(*elem_parts[iep]) )
                    continue;

                  stk::mesh::PartVector::iterator found_elem_part = std::find(found->second.begin(), found->second.end(), elem_parts[iep]);
                  if (found_elem_part != found->second.end())
                    {
                      if (DEBUG_GSPR & 0)
                        {
                          std::cout << "connectSidesForced: found side/elem parts = " << found->first->name() << " " << (*found_elem_part)->name() << std::endl;
                        }
                      valid = true;
                      break;
                    }
                }
            }
          valid_side_part_map = valid;
          if (!valid) return std::pair<bool,bool>(false,false);
        }

      // special case for shells
      int topoDim = UniformRefinerPatternBase::getTopoDim(element_topo);

      bool isShell = false;
      if (topoDim < (int)m_eMesh.entity_rank(element))
        {
          isShell = true;
        }
      int spatialDim = m_eMesh.get_spatial_dim();

      if (spatialDim == 3 && m_eMesh.entity_rank(side) == m_eMesh.edge_rank())
        {
          element_nsides = (unsigned) element_topo.getEdgeCount();
        }

      int permIndex = -1;
      int permPolarity = 1;

      unsigned k_element_side = 0;

      // try search
      for (unsigned j_element_side = 0; j_element_side < element_nsides; j_element_side++)
        {
          m_eMesh.element_side_permutation(element, side, j_element_side, permIndex, permPolarity, use_coordinate_compare, debug && false);
          if (permIndex >= 0)
            {
              k_element_side = j_element_side;
              break;
            }
        }

      if (DEBUG_GSPR && permPolarity < 0)
        {
          std::cout << "permPolarity < 0" << std::endl;
        }
      bool ldb = false;

      if (ldb && permIndex >= 0 && (onlyPosPerm == (permPolarity > 0)))
        {
          std::cout << "P[" << m_eMesh.get_rank() << "] TMP1 srk element= " << m_eMesh.id(element) << " G[elem]= " << m_eMesh.isGhostElement(element) << " isShell= " << isShell << " element_nsides= " << element_nsides
                    << " topoDim= " << topoDim << " element_topo.name= " << element_topo.getName()
                    << " erank = " << m_eMesh.entity_rank(element)
                    << " permPolarity = " << permPolarity << " permIndex= " << permIndex
                    << std::endl;
        }


      bool allowNegPolarity = true;
      if (!isShell && (!allowNegPolarity && permPolarity < 0))
        return std::pair<bool,bool>(false,false);

      if (permIndex >= 0)
        {

          // special case for shells
          if (isShell)
            {
              // FIXME for 2D
              if (m_eMesh.entity_rank(side) == m_eMesh.face_rank())
                {
                  percept::MyPairIterRelation elem_sides (m_eMesh, element, m_eMesh.entity_rank(side));
                  unsigned elem_sides_size= elem_sides.size();
                  if (debug) {
                    std::cout << "tmp srk found shell, elem_sides_size= " << elem_sides_size << " k_element_side= " << k_element_side << std::endl;
                  }
                  if (elem_sides_size == 1)
                    {
                      stk::mesh::RelationIdentifier rel_id = elem_sides[0].relation_ordinal();
                      if (rel_id > 1)
                        throw std::logic_error("connectSidesForced:: logic 1");
                      k_element_side = (rel_id == 0 ? 1 : 0);
                      if (debug) std::cout << "tmp srk k_element_side= " << k_element_side << " rel_id= " << rel_id << std::endl;
                    }
                }
            }

          int exists=0;
          //percept::MyPairIterRelation elem_sides (m_eMesh, element, m_eMesh.entity_rank(side));
          const stk::mesh::Entity *elem_sides = m_eMesh.get_bulk_data()->begin(element, m_eMesh.entity_rank(side));
          const unsigned elem_sides_size = m_eMesh.get_bulk_data()->num_connectivity(element, m_eMesh.entity_rank(side));
          const stk::mesh::ConnectivityOrdinal *elem_sides_ordinals = m_eMesh.get_bulk_data()->begin_ordinals(element, m_eMesh.entity_rank(side));

          for (unsigned iside=0; iside < elem_sides_size; iside++)
            {
              stk::mesh::Entity existing_side = elem_sides[iside];
              if (existing_side == side)
                {
                  ++exists;
                }

              if (elem_sides_ordinals[iside] == k_element_side ) {

                // for debugging only, will be deleted after declare_element_side changes
#if 0
                if (s_allow_duplicates)
                  std::cout << "WARNING: ";
                else
                  std::cout << "ERROR: ";
                std::cout << " Relation already exists: connectSidesForced elem_sides_size= " << elem_sides_size << " element= "; m_eMesh.print(element, true, true);
#if  defined(STK_PERCEPT_HAS_GEOMETRY)
                std::cout << " side= " << m_eMesh.identifier(side) << " in_geom= " << m_eMesh.is_in_geometry_parts(m_geomFile, m_eMesh.bucket(side)); m_eMesh.print(side, true, true);
                std::cout << " existing_side= " << m_eMesh.identifier(existing_side) << " in_geom= " << m_eMesh.is_in_geometry_parts(m_geomFile, m_eMesh.bucket(existing_side)); m_eMesh.print(existing_side, true, true);
#endif
                if (!s_allow_duplicates)
                  VERIFY_OP_ON(elem_sides[iside].relation_ordinal(), !=, k_element_side, "Relation already exists!");
#endif
              }

            }

          if (k_element_side_ret)
            *k_element_side_ret = static_cast<stk::mesh::ConnectivityOrdinal> (k_element_side);

          return std::pair<bool,bool>(true, exists != 0);
        }
      else
        {
          return std::pair<bool,bool>(false, false);
        }
    }

    void FixSideSets::disconnect_entity(stk::mesh::Entity entity)
    {
      stk::mesh::BulkData & mesh = *m_eMesh.get_bulk_data();
      stk::mesh::EntityRank entity_rank = mesh.entity_rank(entity);
      static std::vector<stk::mesh::Entity> relatives, rel;
      static std::vector<stk::mesh::ConnectivityOrdinal> relative_ordinals, relo;

      for (stk::mesh::EntityRank irank = stk::topology::ELEMENT_RANK; irank != entity_rank; --irank)
        {
          // Previously this attempted to delete forward or backward and still the list got corrupted,
          // so just copy into vector and delete from there.
          rel.assign(mesh.begin(entity, irank),mesh.end(entity, irank));
          relo.assign(mesh.begin_ordinals(entity, irank), mesh.end_ordinals(entity, irank));
          relatives.resize(0);
          relative_ordinals.resize(0);
          for (size_t irel = 0; irel < rel.size(); ++irel)
            {
              if (mesh.bucket(rel[irel]).owned())
                {
                  relatives.push_back(rel[irel]);
                  relative_ordinals.push_back(relo[irel]);
                }
            }

          for (size_t irel = 0; irel < relatives.size(); ++irel)
            {
              mesh.destroy_relation( relatives[irel], entity, relative_ordinals[irel]);
            }
        }
    }


    void FixSideSets::delete_unattached_sides(SetOfEntities& side_set, SetOfEntities *avoid_sides)
    {
      LocalTimer timer_bss("FSS0_delete_unattached_sides", m_refiner);

      SetOfEntities sides_to_delete(*m_eMesh.get_bulk_data());

      for (SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
        {
          if (m_eMesh.owned(*it_side) && m_eMesh.get_bulk_data()->num_connectivity(*it_side, m_eMesh.element_rank()) == 0)
            {
              sides_to_delete.insert(*it_side);
            }
        }

      print_set(m_eMesh, sides_to_delete, "sides_to_delete");

      // add pseudo-sides to ensure STK doesn't delete Field values on nodes
      if (avoid_sides)
        {
          std::vector<stk::mesh::Entity> new_elements;
          m_eMesh.getEntitiesUsingIdServer(m_eMesh.side_rank(), sides_to_delete.size(), new_elements);
          size_t ielem = 0;
          for (SetOfEntities::iterator siter = sides_to_delete.begin(); siter != sides_to_delete.end(); ++siter, ++ielem)
            {
              stk::mesh::Entity element = *siter;
              stk::mesh::Entity newElement = new_elements[ielem];
              const stk::mesh::PartVector& super = m_eMesh.bucket(element).supersets();
              stk::mesh::PartVector add, rem;
              for (unsigned ipart = 0; ipart < super.size(); ++ipart)
                {
                  stk::mesh::Part *  part = super[ipart];

                  if ( stk::mesh::is_auto_declared_part(*part) )
                    continue;
                  bool is_auto_part = part->attribute<AutoPart>() != 0;
                  if (is_auto_part)
                    continue;
                  add.push_back(part);
                }
              const percept::MyPairIterRelation elem_nodes (m_eMesh, element, m_eMesh.node_rank());
              for (unsigned ii=0; ii < elem_nodes.size(); ++ii)
                {
                  m_eMesh.get_bulk_data()->declare_relation(newElement, elem_nodes[ii].entity(), ii);
                }

              m_eMesh.get_bulk_data()->change_entity_parts( newElement, add, rem );

            }
          avoid_sides->insert(new_elements.begin(), new_elements.end());
        }

      const stk::mesh::EntityRank FAMILY_TREE_RANK = static_cast<stk::mesh::EntityRank>(stk::topology::ELEMENT_RANK + 1u);
      SetOfEntities family_trees(*m_eMesh.get_bulk_data());
      for (SetOfEntities::iterator siter=sides_to_delete.begin(); siter != sides_to_delete.end(); ++siter)
        {
          stk::mesh::Entity side = *siter;

          while (true)
            {
              percept::MyPairIterRelation rels (m_eMesh, side, FAMILY_TREE_RANK);
              if (!rels.size())
                break;
              stk::mesh::Entity to_rel = rels[0].entity();
              family_trees.insert(to_rel);
              stk::mesh::RelationIdentifier to_id = rels[0].relation_ordinal();

              bool del = m_eMesh.get_bulk_data()->destroy_relation( to_rel, side, to_id);
              if (!del)
                throw std::runtime_error("fix_side_sets_2:: destroy_relation failed 4");
            }

          if ( ! m_eMesh.get_bulk_data()->destroy_entity( side ) )
            {
              throw std::runtime_error("fix_side_sets_2 error 4 - couldn't delete");
            }
        }
      for (SetOfEntities::iterator fiter=family_trees.begin(); fiter != family_trees.end(); fiter++)
        {
          stk::mesh::Entity family_tree = *fiter;
          percept::MyPairIterRelation rels (m_eMesh, family_tree, m_eMesh.side_rank());
          if (rels.size() == 1)
            {
              if ( ! m_eMesh.get_bulk_data()->destroy_entity( family_tree ) )
                {
                  throw std::runtime_error("fix_side_sets_2 error 4.1 - couldn't delete family_tree");
                }
            }
        }
    }


    bool FixSideSets::bucket_acceptable(stk::mesh::Bucket& bucket, stk::mesh::EntityRank rank)
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

    void FixSideSets::build_side_set(SetOfEntities& side_set, bool only_roots)
    {
      //LocalTimer timer_bss("FSS0_build_side_set", m_refiner);

      side_set.clear();

      stk::mesh::Selector excludeSelector;
      if (m_excludeParts.size())
        {
          excludeSelector = stk::mesh::selectUnion(m_excludeParts);
        }

      stk::mesh::EntityRank side_rank_iter_begin = m_eMesh.side_rank();
      stk::mesh::EntityRank side_rank_iter_end = m_eMesh.side_rank();
      if (m_eMesh.get_spatial_dim() == 3)
        {
          side_rank_iter_begin = m_eMesh.edge_rank();
        }

      for (stk::mesh::EntityRank side_rank_iter = side_rank_iter_begin; side_rank_iter <= side_rank_iter_end; side_rank_iter++)
        {
          const stk::mesh::BucketVector & side_buckets = m_eMesh.get_bulk_data()->buckets( side_rank_iter );
          for ( stk::mesh::BucketVector::const_iterator it_side_bucket = side_buckets.begin() ; it_side_bucket != side_buckets.end() ; ++it_side_bucket )
            {
              stk::mesh::Bucket & side_bucket = **it_side_bucket ;

              if (!bucket_acceptable(side_bucket, side_rank_iter))
                {
                  continue;
                }

              if (m_excludeParts.size() && excludeSelector(side_bucket))
                {
                  continue;
                }
#if  defined(STK_PERCEPT_HAS_GEOMETRY)
              if (m_eMesh.is_in_geometry_parts(m_geomFile, side_bucket))
                {
                  continue;
                }
#endif
              const unsigned num_elements_in_side_bucket = side_bucket.size();
              for (unsigned i_side = 0; i_side < num_elements_in_side_bucket; i_side++)
                {
                  stk::mesh::Entity side = side_bucket[i_side];

                  if (only_roots)
                    {
                      if (m_eMesh.hasFamilyTree(side))
                        {
                          stk::mesh::Entity parent = m_eMesh.getParent(side, false);
                          if (m_eMesh.is_valid(parent))
                            continue;
                        }
                    }

                  bool has_connected_nodes = (0 != m_eMesh.get_bulk_data()->num_connectivity(side, m_eMesh.node_rank()));

                  if (has_connected_nodes)
                    {
                      if (m_buildSideSetSelector)
                        {
                          if (m_buildSideSetSelector->use_batch_filter() || (*m_buildSideSetSelector)(side))
                            {
                              side_set.insert(side);
                            }
                        }
                      else
                        {
                          side_set.insert(side);
                        }
                    }
                }
            }
          if (m_buildSideSetSelector && m_buildSideSetSelector->use_batch_filter())
            {
              m_buildSideSetSelector->batch_filter(side_rank_iter, side_set);
            }
        }
    }

    void FixSideSets::reconnect_sides(SetOfEntities& side_set, SetOfEntities *avoid_elems, bool onlyPosPerm)
    {
      LocalTimer timer_bss("FSS0_reconn_sides", m_refiner);

      size_t count = 0;
      for (SetOfEntities::iterator it_side=side_set.begin(); it_side != side_set.end(); ++it_side)
        {
          stk::mesh::Entity side = *it_side;

          bool valid_side_part_map = false;

          VERIFY_OP_ON(m_eMesh.is_valid(side), ==, true, "bad side");
          bool didConnect = connect(side, valid_side_part_map, avoid_elems, onlyPosPerm);
          if (didConnect)
            ++count;
        }

      if (0)
        {
          size_t sz = side_set.size();
          stk::all_reduce( m_eMesh.parallel(), stk::ReduceSum<1>( &count ) );
          stk::all_reduce( m_eMesh.parallel(), stk::ReduceSum<1>( &sz ) );

          if (!m_eMesh.get_rank())
            std::cout << m_eMesh.rank() << " count= " << count << " size= " << sz << std::endl;
        }

    }

    void FixSideSets::check_connect(SetOfEntities& side_set, SetOfEntities *avoid_elems)
    {
      LocalTimer timer_bss("FSS0_check_connect", m_refiner);

      std::vector<stk::mesh::Entity> node_vec;
      for (SetOfEntities::iterator it_side = side_set.begin(); it_side != side_set.end(); ++it_side)
        {
          stk::mesh::Entity side = *it_side;

          if (m_eMesh.owned(side) && m_eMesh.get_bulk_data()->num_connectivity(side, m_eMesh.element_rank()) == 0)
            {
              std::cout << "P[" << m_eMesh.get_rank() << "] ERROR: from: "
                        << m_eMesh.getProperty("FixSideSets::fix_side_sets_2") << " side = " << m_eMesh.id(side) << " side-is-leaf? " << m_eMesh.isLeafElement(side)
                        << std::endl;
              m_eMesh.print(side);
              s_debugCSF = true;
              bool valid_side_part_map = false;
              //bool found =
              connect(side, valid_side_part_map, avoid_elems);

              stk::mesh::PartVector const& side_parts = m_eMesh.bucket(side).supersets();
              for (unsigned isp = 0; isp < side_parts.size(); isp++)
                {
                  std::cout << "side parts= " << side_parts[isp]->name() << std::endl;
                }
#if 0
              stk::mesh::PartVector elem_parts;
              m_eMesh.bucket(elem).supersets(elem_parts);
              for (unsigned isp = 0; isp < elem_parts.size(); isp++)
                {
                  std::cout << "elem parts= " << elem_parts[isp]->name() << std::endl;
                }
#endif
              //percept::MyPairIterRelation side_nodes (m_eMesh, side, m_eMesh.node_rank());
              stk::mesh::Entity const *side_nodes = m_eMesh.get_bulk_data()->begin_nodes(side);
              unsigned side_nodes_size = m_eMesh.get_bulk_data()->num_nodes(side);
              node_vec.resize(0);
              for (unsigned isnode=0; isnode < side_nodes_size; isnode++)
                {
                  node_vec.push_back(side_nodes[isnode]);
                }
              m_eMesh.dump_vtk(node_vec, "fix_side_sets_2_error2.vtk");
              std::set<stk::mesh::Entity> side_set_0_debug;
              side_set_0_debug.insert(side);
              m_eMesh.dump_vtk("fix_side_sets_2_error2_side.vtk", true, &side_set_0_debug);
              std::set<stk::mesh::Entity> elems_debug;
              m_eMesh.get_node_neighbors(side, elems_debug);
              std::cout << "elems.size= " << elems_debug.size() << std::endl;
              m_eMesh.filter_active_only(elems_debug);
              std::cout << "elems.size= " << elems_debug.size() << std::endl;
              elems_debug.insert(side);
              m_eMesh.dump_vtk("fix_side_sets_2_error2_all.vtk", true, &elems_debug);
              std::cout << PerceptMesh::demangled_stacktrace(30, true, "Refiner") << std::endl;
              VERIFY_MSG("fix_side_sets_2 error 2");
            }
        }

    }
    void FixSideSets::end_begin(const std::string& msg)
    {
      LocalTimer timer_bss("FSS0_end_begin", m_refiner);

      if (m_refiner)
        {
          m_refiner->mod_end(0,"FSS0"+msg);
          m_refiner->mod_begin();
        }
      else
        {
          stk::mesh::fixup_ghosted_to_shared_nodes(*m_eMesh.get_bulk_data());
          m_eMesh.get_bulk_data()->modification_end();
          m_eMesh.get_bulk_data()->modification_begin();
        }
    }

    void FixSideSets::doProgressPrint(PerceptMesh& eMesh, const std::string& msg)
    {
      if (m_doProgress && eMesh.get_rank() == 0)
        {
          //std::cout << msg << " cpu: " << eMesh.cpu_time() << " [sec] mem= " << eMesh.print_memory_both() << std::endl;
          size_t now=0, hwm=0;
          stk::get_memory_usage(now, hwm);
          std::cout << msg << " cpu: " << eMesh.cpu_time() << " [sec] mem= " << stk::human_bytes(now) << " [now_proc0] " << stk::human_bytes(hwm) << " [hwm_proc0]" << std::endl;
        }
    }

    static void print_surface_blocks_map(PerceptMesh& m_eMesh)
    {
      std::vector<const stk::mesh::Part*> surfaces = m_eMesh.get_fem_meta_data()->get_surfaces_in_surface_to_block_map();
      for (auto surface : surfaces)
        {
          std::vector<const stk::mesh::Part*> blocks = m_eMesh.get_fem_meta_data()->get_blocks_touching_surface(surface);
          for (auto block : blocks)
            {
              std::cout << "print_surface_blocks_map: surface= " << surface->name() << " block= " << block->name() << std::endl;
            }
        }

    }

    /**
     * 1. loop over sides - for each volume connected to side, find parts, find leaf surfaces it is in
     * 2. if it is in one of the urpconv surfaces, remove it from there
     * 3. add to parts as below
     */

    void FixSideSets::move_sides_to_correct_surfaces()
    {
      if (m_eMesh.get_spatial_dim() == 2)
        return;

      const std::string& append_conv_string = UniformRefinerPatternBase::getAppendConvertString();
      (void)append_conv_string;

      const bool debug = false;
      if (debug)
        {
          print_surface_blocks_map(m_eMesh);
        }

      stk::mesh::PartVector pv = m_eMesh.get_fem_meta_data()->get_mesh_parts();
      for (auto partp : pv)
        {
          if (debug) std::cout << "move_sides_to_correct_surfaces: processing surface= " << partp->name()
                               << " topo= " << partp->topology() << " primary_entity_rank= " << partp->primary_entity_rank() << std::endl;

          if (partp->topology() == stk::topology::INVALID_TOPOLOGY)
            continue;
          if (partp->primary_entity_rank() != m_eMesh.side_rank())
            continue;

          stk::mesh::EntityVector sides;
          stk::mesh::Selector sel = stk::mesh::Selector(*partp) & m_eMesh.get_fem_meta_data()->locally_owned_part();
          stk::mesh::get_selected_entities(sel , m_eMesh.get_bulk_data()->buckets(m_eMesh.side_rank()), sides);

          for (auto side : sides)
            {
              MyPairIterRelation volumes(m_eMesh, side, m_eMesh.element_rank());
              for (unsigned iv = 0; iv < volumes.size(); ++iv)
                {
                  move_side_to_correct_surface(*partp, side, volumes[iv].entity());
                }
            }
        }
    }

    void FixSideSets::move_side_to_correct_surface(stk::mesh::Part& surface, stk::mesh::Entity side, stk::mesh::Entity volume)
    {
      const bool debug = false;
      const std::string& append_conv_string = UniformRefinerPatternBase::getAppendConvertString();
      (void)append_conv_string;

      stk::topology side_topo = m_eMesh.get_bulk_data()->bucket(side).topology();
      stk::topology elem_topo = m_eMesh.get_bulk_data()->bucket(volume).topology();

      if (debug) std::cout << "side,elem topo = " << side_topo << "," << elem_topo << std::endl;

      //std::string part_name = add_parts[0]->name(); // surface_hex8_quad4_1

      std::vector<stk::mesh::Part*> add_parts(1, static_cast<stk::mesh::Part *>(0));
      std::vector<stk::mesh::Part*> remove_parts(1, &surface);

      std::string new_part_name;

      std::vector<std::string> tokens;
      stk::util::tokenize(surface.name(), "_", tokens);

      // this can happen, eg. for generated meshes - with a surface name like surface_1_quad4
      if (tokens.size() <= 3)
        return;

      if (debug) {
          std::cout << "tokens = ";
          for (unsigned i=0; i<tokens.size(); i++)
              std::cout << tokens[i] << " ";
          std::cout << std::endl;
      }
      unsigned dot_pos = tokens[3].find(".");
      unsigned nl = tokens[3].length();
      if (dot_pos != std::string::npos)
        nl = dot_pos;

      // FIXME: substr usage - only works for single digit

      convert_stk_topology_to_ioss_name(elem_topo, tokens[1]);
      convert_stk_topology_to_ioss_name(side_topo, tokens[2]);

      tokens[3] = tokens[3].substr(0, nl);
      new_part_name = tokens[0];
      for (unsigned i = 1; i < 4; i++)
        new_part_name += "_" + tokens[i];

      if (debug)
          std::cout << "new_part_name=" << new_part_name << std::endl;

      add_parts[0] = m_eMesh.get_fem_meta_data()->get_part(new_part_name);
      if (!add_parts[0])
        {
          if (debug) std::cout << " new part name not found, skipping" << std::endl;
          //std::cout << "add_parts[0] = null, new_part_name= " << new_part_name << " elem_topo= " << elem_topo << " side_topo= " << side_topo << std::endl;
          //VERIFY_MSG("bad add_parts");
          return;
        }

      if (debug)
        {
          std::cout << "moving side= " << m_eMesh.id(side) << " side_topo= " << side_topo << " attached to vol= " << m_eMesh.id(volume)
                    << " of topo= " << elem_topo << " from surface= " << surface.name() << " to: " << add_parts[0]->name() << " remove_parts.size= " << remove_parts.size()
                    << std::endl;
          std::cout << "before parts= " << m_eMesh.print_entity_parts_string(side, "\n") << std::endl;
        }

      stk::mesh::PartVector supersets = m_eMesh.bucket(side).supersets();
      for (auto superset : supersets)
        {

          bool isConvertedPart = ( superset->name().find(append_conv_string) != std::string::npos);
          if (isConvertedPart && superset->primary_entity_rank() == m_eMesh.side_rank())
            {
              remove_parts[0] = superset;
            }
        }

      if (add_parts[0] == remove_parts[0])
        remove_parts.resize(0);

      if (debug && remove_parts.size()) std::cout << "found remove_parts = " << remove_parts[0]->name() << std::endl;

      m_eMesh.get_bulk_data()->change_entity_parts( side, add_parts, remove_parts );

      if (debug) std::cout << "after parts= " << m_eMesh.print_entity_parts_string(side, "\n") << std::endl;

    }


    // fast reconnector

    void FixSideSets::
    fix_side_sets_2(bool allow_not_found, SetOfEntities *avoid_elems, SetOfEntities *avoid_sides, const std::string& msg)
    {
      doProgressPrint(m_eMesh, "Stage: FixSideSets 1");
      s_allow_duplicates = false;
      if (m_eMesh.getProperty("Refiner_connect_allow_duplicate_sides") == "true")
        s_allow_duplicates = true;

      s_onlyLeafToLeaf = S_ONLYLEAFTOLEAF_DEFAULT;
      if (m_eMesh.getProperty("percept_Refiner_connect_any_to_any") == "true")
        s_onlyLeafToLeaf = false;

      bool reduced_mod_end = true;
      if (m_eMesh.getProperty("percept_reduced_mod_end") == "false")
        reduced_mod_end = false;
      (void)reduced_mod_end;

      if (LTRACE) std::cout << m_eMesh.rank() << "tmp fix_side_sets_2 start... from: " << m_eMesh.getProperty("FixSideSets::fix_side_sets_2") << std::endl;

      // loop over all sides that are leaves (not parent or have no family tree),
      //   loop over their nodes and their associated elements,
      //     connect element and side if they share a face

      SetOfEntities side_set(*m_eMesh.get_bulk_data());

      if (avoid_sides)
        VERIFY_OP_ON(avoid_sides->size(), ==, 0, "bad avoid_sides");

      doProgressPrint(m_eMesh, "Stage: FixSideSets 2");

      if (LTRACE) std::cout << m_eMesh.rank() << "tmp fix_side_sets_2 1st end_begin: " << m_eMesh.getProperty("FixSideSets::fix_side_sets_2") << std::endl;

      build_side_set(side_set);
      if (LTRACE) std::cout << m_eMesh.rank() << "reconnect_sides side_set.size= " << side_set.size() << std::endl;
      reconnect_sides(side_set, avoid_elems, false);

      end_begin(msg+"-Re-Side");

      build_side_set(side_set);
      if (LTRACE) std::cout << m_eMesh.rank() << "after reconnect_sides, allow_not_found= " << allow_not_found << " side_set.size= " << side_set.size() << std::endl;

      m_eMesh.initializeIdServer();

      if (!allow_not_found)
        check_connect(side_set, avoid_elems);
      else
        delete_unattached_sides(side_set, avoid_sides);

      if (LTRACE) std::cout << m_eMesh.rank() << "fix_permutation side_set.size= " << side_set.size() << std::endl;

      doProgressPrint(m_eMesh, "Stage: FixSideSets 3");

      if (1)
        {
          build_side_set(side_set);
          fix_permutation(side_set);
        }

      if (LTRACE) std::cout << m_eMesh.rank() << "tmp fix_side_sets_2 ...end from: " << m_eMesh.getProperty("FixSideSets::fix_side_sets_2") << std::endl;

      end_begin(msg+"moveSides");
      move_sides_to_correct_surfaces();
    }

}
