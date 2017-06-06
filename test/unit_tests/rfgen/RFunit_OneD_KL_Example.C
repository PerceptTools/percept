#include "RFunit_OneD_KL_Example.h"

#include <Teuchos_TestForException.hpp>

#include <gtest/gtest.h>

#include <mpi.h>

namespace RFGen
{

OneD_KLSolver::OneD_KLSolver(
  const int num_intervals,
  std::vector<double> &eigenVal_mem,
  std::vector<double> &eigenVec_mem)
  : 
  API_KLSolver(),
  num_intervals_(num_intervals),
  h_(1.0/(double)num_intervals),
  stochasticDim_(0),
  mpiComm_(MPI_COMM_WORLD),
  eigenVal_mem_(eigenVal_mem),
  eigenVec_mem_(eigenVec_mem)
{}

void 
OneD_KLSolver::computeLocalIntgDataSizes(
  int &localNumElem,
  int &localMaxIntgPts)
{
  int numProcs, myRank;
  EXPECT_TRUE(MPI_SUCCESS == MPI_Comm_size(mpiComm_, &numProcs));
  EXPECT_TRUE(MPI_SUCCESS == MPI_Comm_rank(mpiComm_, &myRank));

  intervals_per_proc_ = num_intervals_/numProcs;

  if (myRank != numProcs-1)
    localNumElem = intervals_per_proc_;
  else
    localNumElem = num_intervals_ - (numProcs-1)*intervals_per_proc_;

  localMaxIntgPts = 1;

  num_local_intervals_ = localNumElem;
}

void
OneD_KLSolver::computeLocalIntgData(
  Array<double,NaturalOrder,Cell,Point,Dim> &localIntgPtCoords,
  Array<double,NaturalOrder,Cell,Point> &localVolumeWeights)
{
  EXPECT_TRUE(num_local_intervals_==localIntgPtCoords.dimension(0));
  EXPECT_TRUE(1==localIntgPtCoords.dimension(1));
  EXPECT_TRUE(1==localIntgPtCoords.dimension(2));

  EXPECT_TRUE(num_local_intervals_==localVolumeWeights.dimension(0));
  EXPECT_TRUE(1==localVolumeWeights.dimension(1));

  int myRank;
  EXPECT_TRUE(MPI_SUCCESS == MPI_Comm_rank(mpiComm_, &myRank));
  
  const int offset = myRank * intervals_per_proc_;

  for (int i=0; i<num_local_intervals_; i++)
  {
    localIntgPtCoords(i,0,0) = (offset+i+0.5)*h_;
    localVolumeWeights(i,0) = h_;
  }
}

MPI_Comm
OneD_KLSolver::getParallelComm() const
{
  return mpiComm_;
}

void 
OneD_KLSolver::setKLSolution(
  const int &stochasticDim,
  const Array<double,NaturalOrder,Eigen> &eigenValues,
  const Array<double,NaturalOrder,Eigen,Cell> &eigenVectors)
{
  stochasticDim_  = stochasticDim;

  // FIXME this would be simpler if there were an "ArrayContainer" class 
  // that could manage its own memory

  eigenVal_mem_.resize(stochasticDim_);
  eigenVec_mem_.resize(stochasticDim_*num_local_intervals_);

  memcpy(&eigenVal_mem_[0], eigenValues.contiguous_data(), stochasticDim_*sizeof(double));
  eigenValues_.assign(&eigenVal_mem_[0], stochasticDim_);

  memcpy(&eigenVec_mem_[0], eigenVectors.contiguous_data(), stochasticDim_*num_local_intervals_*sizeof(double));
  eigenVectors_.assign(&eigenVec_mem_[0], stochasticDim_, num_local_intervals_);
}

}
