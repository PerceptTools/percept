// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GenericAlgorithm_update_coordinates_hpp
#define GenericAlgorithm_update_coordinates_hpp

#include <percept/MeshType.hpp>

namespace percept {

template<typename MeshType>
class ReferenceMeshSmootherConjugateGradientImpl;

#if defined(WITH_KOKKOS) //&& defined(KOKKOS_HAVE_OPENMP)
using Array4D = viewType;

template<typename MeshType>
struct A4DMD{
	unsigned blkIndex;
	unsigned numNodes;
	Kokkos::View<typename MeshType::MTNode*, DataLayout, MemSpace> blockNodes;

	A4DMD(){
		blkIndex =0;
		numNodes=0;
	}

	~A4DMD(){}

};
#else
using Array4D = MDArray;

template<typename MeshType>
struct A4DMD{
	//metadata on an A4D
	unsigned blkIndex;
	unsigned numNodes;
	std::vector<typename MeshType::MTNode> blockNodes;

	A4DMD(){

		blockNodes.resize(0);
		blkIndex =0;
		numNodes=0;
	}

	~A4DMD(){}
};
#endif

  template<typename MeshType>
  struct GenericAlgorithm_update_coordinates
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using Double = typename RefMeshSmoother::Double;
    using This = GenericAlgorithm_update_coordinates<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    int spatialDim;
    typename MeshType::MTField *coord_field;
    typename MeshType::MTField *coord_field_current;
    typename MeshType::MTField *coord_field_lagged;

    typename MeshType::MTField *cg_s_field;


    std::vector<typename MeshType::MTNode> nodes;


    Double alpha;

    GenericAlgorithm_update_coordinates(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in);

#if defined (WITH_KOKKOS) //&& defined(KOKKOS_HAVE_OPENMP)
    std::list<A4DMD<MeshType>> blockMetaDatas;
	Kokkos::View<typename MeshType::MTNode*, DataLayout, MemSpace> nodesThisBlock; //will be used to loop over a particular block's nodes.
	Array4D cfc;
	Array4D cfl;
	Array4D cgsf;
#else
    std::list<A4DMD<MeshType>> blockMetaDatas;
    std::vector<typename MeshType::MTNode> * nodesThisBlock;
	Array4D * cfc;
	Array4D * cfl;
	Array4D * cgsf;
#endif

#if defined (WITH_KOKKOS) //&& defined(KOKKOS_HAVE_OPENMP)
    KOKKOS_INLINE_FUNCTION
    void operator()(const int64_t& index) const;
#else
    inline
    void operator()(const int64_t& index) const;
#endif

    void run();

  };

  #if defined (WITH_KOKKOS) //&& defined(KOKKOS_HAVE_OPENMP)

  struct simplified_gatm_1_BSG{ //uses A4D's directly

    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid>;
	  using This = simplified_gatm_1_BSG;
	  typedef long double Double;
	  using Array4D = viewType;

	  RefMeshSmoother * m_rms; //adding smoother pointer a member didn't effect performance
	  PerceptMesh * m_eMesh;

    StructuredGrid::MTSelector on_locally_owned_part;
    StructuredGrid::MTSelector on_globally_shared_part;
    int spatialDim;
    MTSGridField * coord_field;
	  MTSGridField * coord_field_current;
	  MTSGridField * coord_field_lagged;
	  MTSGridField * cg_s_field;


	  std::vector<StructuredGrid::MTNode> nodes;

	  Double alpha;

	  std::list<A4DMD<StructuredGrid>> blockMetaDatas;
	  Kokkos::View<StructuredGrid::MTNode*, DataLayout, MemSpace> nodesThisBlock;

	  Array4D cfc;
	  Array4D cfl;
	  Array4D cgsf;



	  Array4D cfcB;
	  Array4D cflB;
	  Array4D cgsfB;
	  Kokkos::View< StructuredCellIndex*, DataLayout , MemSpace > nodesV;
	  unsigned numNodes;
	  unsigned numBlocks;

	  simplified_gatm_1_BSG(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in);

  	KOKKOS_INLINE_FUNCTION
  	void
  	operator()(int64_t index) const
  	{
  		int i  = nodesV(index)[0];
  		int j  = nodesV(index)[1];
  		int k  = nodesV(index)[2];

  		for (int iDim=0; iDim < spatialDim; iDim++)
  	        {
  			  Double dt = alpha*cgsfB(i,j,k,iDim);
  	  	  	  cfcB(i,j,k,iDim) += dt;
  	        }
  	}

	void run();
  };





        //testing to see if templating degrades performance?

template<typename boople>
struct templeKokkos {
	Kokkos::View<double*, DataLayout , MemSpace > templeDubs;
	unsigned size;

	templeKokkos(unsigned into);

	KOKKOS_INLINE_FUNCTION
	void
	operator()(int64_t index) const;

	void run();
};

#endif

} // namespace percept

#endif
