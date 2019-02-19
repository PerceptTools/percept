/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef RFunit_OneD_KL_Example_h
#define RFunit_OneD_KL_Example_h

#include "percept/rfgen/RFGen_API_KLSolver.h"

namespace RFGen
{

// class to implement simple discrete spatial RF based on OneD domain (interval)

using shards::Array;
using shards::NaturalOrder;

class OneD_KLSolver : public API_KLSolver
{
 public: 
  explicit
  OneD_KLSolver(
    const int num_intervals,
    std::vector<double> &eigenVal_mem,
    std::vector<double> &eigenVec_mem);

  unsigned getSpatialDim() const {return 1;}

  void computeLocalIntgDataSizes(
    int &localNumElem,
    int &localMaxIntgPts);

  void computeLocalIntgData(
    Array<double,NaturalOrder,Cell,Point,Dim> &localIntgPtCoords,
    Array<double,NaturalOrder,Cell,Point> &localVolumeWeights);

  MPI_Comm getParallelComm() const;

  void setKLSolution(
    const int &numTerms,
    const Array<double,NaturalOrder,Eigen> &eigenValues,
    const Array<double,NaturalOrder,Eigen,Cell> &eigenVectors);

 private:
  const int num_intervals_;
  const double h_;

  int intervals_per_proc_, num_local_intervals_;

  int stochasticDim_;
  Array<double,NaturalOrder,Eigen> eigenValues_;
  Array<double,NaturalOrder,Eigen,Cell> eigenVectors_;

  Array<double,NaturalOrder,Cell> rfValuesCell_;
  
  MPI_Comm mpiComm_;

  std::vector<double> &eigenVal_mem_;
  std::vector<double> &eigenVec_mem_;
};


}

#endif
