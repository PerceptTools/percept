/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <cmath>

#include "Teuchos_Assert.hpp"

#include "percept/rfgen/RFGen_RealizeKL.h"

namespace RFGen
{

TEST(UQ, RealizeKL)
{
  // application-owned data that is passed to RealizeRF
  // this may represent all the cells owned by a processor or it 
  // may be only a "bucket" of cells from a workset algorithm
  const int numCells = 5;
  const int rv_dim = 2;
  std::vector<double> kl_evals_mem(rv_dim);
  std::vector<double> kl_evecs_mem(rv_dim*numCells);
  std::vector<double> rf_values_mem(numCells);

  // wrap the data (flat arrays) in Array wrappers
  shards::Array<double,shards::NaturalOrder,Eigen>      kl_evals(&kl_evals_mem[0],rv_dim);
  shards::Array<double,shards::NaturalOrder,Cell,Eigen> kl_evecs(&kl_evecs_mem[0],numCells,rv_dim);
  shards::Array<double,shards::NaturalOrder,Cell>       rf_vals(&rf_values_mem[0],numCells);

  // for this test, we make up the data for a 5-cell, 2-term KL solution
  kl_evals(0) = 100;
  kl_evals(1) = 30;

  for (int i=0; i<numCells; i++)
  {
    kl_evecs(i,0) = i+1;
    kl_evecs(i,1) = -i*i;
  }

  // global scalars needed for RF realizations
  const double rf_mean = 1.0;
  const double rf_variance = 2.0;

  // build an instance of the KL series realization class
  RealizeKL realizeKL(rv_dim, rf_mean, rf_variance);

  // create random variable coefficients for a particular realization
  std::vector<double> rv_coeffs_mem(rv_dim);
  shards::Array<double,shards::NaturalOrder,Eigen> rv_coeffs(&rv_coeffs_mem[0], rv_dim);

  rv_coeffs(0) = 1.5;
  rv_coeffs(1) = -2.0;

  // build a specific KL realization using the KL data and rv coeffs
  realizeKL.computeRealization(rv_coeffs, kl_evals, kl_evecs, rf_vals);
  
  // verify the resulting solution
  ASSERT_NEAR(rf_vals(0), 1.0+2.0*(std::sqrt(100.0)*(1.5)*(1.0)+std::sqrt(30.0)*(-2.0)*(0.0)), 1e-12);
  ASSERT_NEAR(rf_vals(2), 1.0+2.0*(std::sqrt(100.0)*(1.5)*(3.0)+std::sqrt(30.0)*(-2.0)*(-4.0)), 1e-12);
  ASSERT_NEAR(rf_vals(4), 1.0+2.0*(std::sqrt(100.0)*(1.5)*(5.0)+std::sqrt(30.0)*(-2.0)*(-16.0)), 1e-12);
}

}
