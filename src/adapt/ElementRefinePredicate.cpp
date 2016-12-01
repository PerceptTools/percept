// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <adapt/ElementRefinePredicate.hpp>

  namespace percept {

    /// Return DO_REFINE, DO_UNREFINE, DO_NOTHING
    int ElementRefinePredicate::operator()(const stk::mesh::Entity entity)
    {
      int mark = 0;

      if (getMarkNone())
        return DO_NOTHING;

      stk::mesh::EntityRank rank = m_eMesh.entity_rank(entity);
      if (rank != m_eMesh.element_rank())
        return DO_NOTHING;

      bool isParent = m_eMesh.hasFamilyTree(entity) && m_eMesh.isParentElement(entity, false);

      RefineFieldType_type *fdata = 0;
      if (m_field && m_field->entity_rank() == m_eMesh.entity_rank(entity))
        {
          //fdata = m_eMesh.field_data( *static_cast<const TransitionElementType *>(m_field) , entity );
          fdata = stk::mesh::field_data( *dynamic_cast<const RefineFieldType *>(m_field) , entity );
        }
      bool selected = (m_eb_selector==0 || (*m_eb_selector)(m_eMesh.bucket(entity)));

      bool ref_field_criterion = false;
      ref_field_criterion = (fdata  && fdata[0] == 1);

      if (getRefineStage() == -10)
        {
          if (m_eMesh.entity_rank(entity) != m_eMesh.element_rank())
            return DO_NOTHING;

          ref_field_criterion = false;
          if (isParent)
            {
              ref_field_criterion = (fdata  && fdata[0] == 2);
            }
        }
      else
        {
          if (isParent) return DO_NOTHING;
        }


      bool unref_field_criterion = (fdata && fdata[0] < 0);
      if (isParent && getRefineStage() == -10)
        {
          unref_field_criterion = false;
        }
      if (unref_field_criterion)
        {
          // can only unrefine elements with parents
          if (!m_eMesh.hasFamilyTree(entity))
            unref_field_criterion = false;
          else
            {
              stk::mesh::Entity parent = m_eMesh.getParent(entity, true);
              if (!m_eMesh.is_valid(parent))
                {
                  unref_field_criterion = false;
                }
            }
        }

      if (selected && ref_field_criterion)
        {
          mark |= DO_REFINE;
        }
      if (selected && unref_field_criterion)
        {
          mark |= DO_UNREFINE;
        }

      // testing
      if (0 && fdata)
        {
          //ScalarFieldType *refine_field = m_eMesh.get_fem_meta_data()->get_field<ScalarFieldType>(stk::topology::ELEMENT_RANK, "refine_field");
          RefineLevelType *refine_level = m_eMesh.get_refine_level_field();

          int *refine_level_elem = 0;
          if (refine_level && refine_level->entity_rank() == m_eMesh.entity_rank(entity))
            refine_level_elem = stk::mesh::field_data( *refine_level , entity );
          RefineLevelType_type *refine_field_elem = fdata;
          if (refine_field_elem[0] > 0)
            {
              if (refine_level_elem[0] > 2)
                {
                  throw std::runtime_error("bad level");
                }
            }
        }

      // FIXME - why is this necessary?
      if (0 && rank < m_eMesh.element_rank())
        {
          //if (m_will_refine[rank])
          if (m_eMesh.side_rank() == rank)
            {
              const percept::MyPairIterRelation elemNeighbors (m_eMesh, entity, m_eMesh.element_rank());
              if (elemNeighbors.size())
                {
                  for (unsigned i=0; i < elemNeighbors.size(); ++i)
                    {
                      int markElemNeighbor = this->operator()(elemNeighbors[i].entity());
                      if (markElemNeighbor & DO_REFINE)  // FIXME - was DO_UNREFINE, for a reason, or was this a bug?
                        {
                          mark |= DO_REFINE;
                        }
                    }
                }
            }
        }
      return mark;
    }


  }

