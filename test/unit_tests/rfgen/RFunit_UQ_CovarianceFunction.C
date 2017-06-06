/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "percept/rfgen/RFGen_CovarianceFunction.h"

namespace RFGen
{

TEST(UQ, CovarianceFunction1D)
{
  /* %TRACE[OFF]% */  /* %TRACE% */
  const double Lx=0.3;
  
  const unsigned dim = 1;
  std::vector<double> scales(dim,Lx);
  
  std::vector<double> x1_mem(dim*2);
  std::vector<double> x2_mem(dim*2);
  std::vector<double> res_mem(2*2);
  
  // TODO: can we create an "ArrayContainer" class to do this for us??
  shards::Array<double,shards::NaturalOrder,Point,Dim> x1(&x1_mem[0], 2, dim);
  shards::Array<double,shards::NaturalOrder,Point,Dim> x2(&x2_mem[0], 2, dim);
  shards::Array<double,shards::NaturalOrder,Point,Point> res(&res_mem[0], 2, 2);
  
  x1(0,0) = 0.4;
  x1(1,0) = 0.5;
  
  x2(0,0) = 0.3;
  x2(1,0) = -0.5;
  
  {
    ExpMagL2CovarianceFunction expMagL2CovarianceFunction(dim, scales);
    expMagL2CovarianceFunction.computeValues(x1, x2, res);
    
    ASSERT_NEAR(res(0,0), std::exp(-std::abs(0.4-0.3)/Lx), 1e-12);
    ASSERT_NEAR(res(1,0), std::exp(-std::abs(0.5-0.3)/Lx), 1e-12);
    ASSERT_NEAR(res(0,1), std::exp(-std::abs(0.4+0.5)/Lx), 1e-12);
    ASSERT_NEAR(res(1,1), std::exp(-std::abs(0.5+0.5)/Lx), 1e-12);
  }

  {
    ExpMagL1CovarianceFunction expMagL1CovarianceFunction(dim, scales);    
    
    expMagL1CovarianceFunction.computeValues(x1, x2, res);
    
    ASSERT_NEAR(res(0,0), std::exp(-std::abs(0.4-0.3)/Lx), 1e-12);
    ASSERT_NEAR(res(1,0), std::exp(-std::abs(0.5-0.3)/Lx), 1e-12);
    ASSERT_NEAR(res(0,1), std::exp(-std::abs(0.4+0.5)/Lx), 1e-12);
    ASSERT_NEAR(res(1,1), std::exp(-std::abs(0.5+0.5)/Lx), 1e-12);
  }
}
  
TEST(UQ, CovarianceFunction2D)
{
  /* %TRACE[OFF]% */  /* %TRACE% */
  const double Lx=0.3;
  const double Ly=0.4;

  const unsigned dim = 2;
  std::vector<double> scales(dim);
  scales[0] = Lx;
  scales[1] = Ly;
  
  std::vector<double> x1_mem(dim*1);
  std::vector<double> x2_mem(dim*1);
  std::vector<double> res_mem(1*1);
  
  shards::Array<double,shards::NaturalOrder,Point,Dim> x1(&x1_mem[0], 1, dim);
  shards::Array<double,shards::NaturalOrder,Point,Dim> x2(&x2_mem[0], 1, dim);
  shards::Array<double,shards::NaturalOrder,Point,Point> res(&res_mem[0], 1, 1);
  
  x1(0,0) = 0.2;
  x1(0,1) = 1.7;
  
  x2(0,0) = 1.3;
  x2(0,1) = 1.5;
  
  {
    ExpMagL2CovarianceFunction expMagL2CovarianceFunction(dim, scales);

    expMagL2CovarianceFunction.computeValues(x1, x2, res);
  
    ASSERT_NEAR(res(0,0), std::exp(-std::sqrt( (0.2-1.3)*(0.2-1.3)/(Lx*Lx) + (1.7-1.5)*(1.7-1.5)/(Ly*Ly) ) ), 1e-12);
  }
  
  {
    ExpMagL1CovarianceFunction expMagL1CovarianceFunction(dim, scales);
    
    expMagL1CovarianceFunction.computeValues(x1, x2, res);
    
    ASSERT_NEAR(res(0,0), std::exp(-(std::abs(0.2-1.3)/Lx + std::abs(1.7-1.5)/Ly ) ), 1e-12);
  }
}

TEST(UQ, CovarianceFunction3D)
{
  /* %TRACE[OFF]% */  /* %TRACE% */
  const double Lx=0.3;
  const double Ly=0.4;
  const double Lz=0.5;
  
  const unsigned dim = 3;
  const unsigned numPoints = 1;
  
  std::vector<double> scales(dim);
  scales[0] = Lx;
  scales[1] = Ly;
  scales[2] = Lz;
  
  std::vector<double> x1_mem(dim*numPoints);
  std::vector<double> x2_mem(dim*numPoints);
  std::vector<double> res_mem(numPoints*numPoints);
  
  shards::Array<double,shards::NaturalOrder,Point,Dim> x1(&x1_mem[0], numPoints, dim);
  shards::Array<double,shards::NaturalOrder,Point,Dim> x2(&x2_mem[0], numPoints, dim);
  shards::Array<double,shards::NaturalOrder,Point,Point> res(&res_mem[0], numPoints, numPoints);
  
  x1(0,0) = 0.4;
  x1(0,1) = 0.1;
  x1(0,2) = 3.2;
  
  x2(0,0) = 0.3;
  x2(0,1) = 0.23;
  x2(0,2) = 2.1;

  {
    ExpMagL2CovarianceFunction expMagL2CovarianceFunction(dim, scales);
    
    expMagL2CovarianceFunction.computeValues(x1, x2, res);
    
    ASSERT_NEAR(res(0,0), std::exp(-std::sqrt( (0.4-0.3)*(0.4-0.3)/(Lx*Lx) + (0.1-0.23)*(0.1-0.23)/(Ly*Ly) + (3.2-2.1)*(3.2-2.1)/(Lz*Lz) ) ), 1e-12);
  }

  {
    ExpMagL1CovarianceFunction expMagL1CovarianceFunction(dim, scales);
    
    expMagL1CovarianceFunction.computeValues(x1, x2, res);
    
    ASSERT_NEAR(res(0,0), std::exp(-(std::abs(0.4-0.3)/Lx + std::abs(0.1-0.23)/Ly + std::abs(3.2-2.1)/Lz ) ), 1e-12);
  }
}


} 
