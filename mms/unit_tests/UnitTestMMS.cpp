/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <mms/mms.hpp>

//#include <valgrind/callgrind.h>

namespace mms {
namespace unit_tests {

//=============================================================================
//=============================================================================
//=============================================================================

TEST(derivatives, firstOrderPartials)
{
//   EXCEPTWATCH;

  const double xb=1.234, yb=2.345, zb=2.3, tb = 23.4;

  FAD_Type x(4, 0, xb), y(4, 1, yb), z(4, 2, zb), t(4, 3, tb);

  FAD_Type expr = (1.0+(x-t)*y*z)+(1.0+sin(x+y-t*z))+(2.0+exp(-x*x-t))+(3.0+pow(y-z*t,3));

  const double dx=DX(expr), dy=DY(expr), dz=DZ(expr), dt=DT(expr);

  // test against output from Mathematica
  EXPECT_DOUBLE_EQ(expr.val(), -136504.5806378632);
  EXPECT_DOUBLE_EQ(dx, 6.393200319571155);
  EXPECT_DOUBLE_EQ(dy, 7899.044775319607);
  EXPECT_DOUBLE_EQ(dz, -186082.60113247877);
  EXPECT_DOUBLE_EQ(dt, -18290.45462323511);
}

TEST(exampleMMS, euler3D)
{
  const double xb=1.234, yb=2.345, zb=2.3, tb = 23.4;

  FAD_Type x(4, 0, xb), y(4, 1, yb), z(4, 2, zb), t(4, 3, tb);

  FAD_Type rho = 1.0+(x-t)*y*z;
  FAD_Type u = 1.0+sin(x+y-t*z);
  FAD_Type v = 2.0+exp(-x*x-t);
  FAD_Type w = 3.0+pow(y-z*t,3);
  FAD_Type T = 1.0+cos(z+x-t*y);

  const double Rgas = 287.0;
  FAD_Type p = Rgas*rho*T;

  const double gamma = 1.4;
  const double Cv = Rgas/(gamma-1.0);
  FAD_Type E = rho*(Cv*T + 0.5*(u*u+v*v+w*w));

  const double mass_src =   DT(rho)   + DX(rho*u)   + DY(rho*v)   + DZ(rho*w);
  const double momx_src =   DT(rho*u) + DX(rho*u*u) + DY(rho*u*v) + DZ(rho*u*w) + DX(p);
  const double momy_src =   DT(rho*v) + DX(rho*v*u) + DY(rho*v*v) + DZ(rho*v*w) + DY(p);
  const double momz_src =   DT(rho*w) + DX(rho*w*u) + DY(rho*w*v) + DZ(rho*w*w) + DZ(p);
  const double energy_src = DT(E)     + DX(u*(E+p)) + DY(v*(E+p)) + DZ(w*(E+p));

  // test against output from Mathematica
  EXPECT_DOUBLE_EQ(mass_src, 2.914077175792222e7);
  EXPECT_DOUBLE_EQ(momx_src, -3.4842036495958334e8);
  EXPECT_DOUBLE_EQ(momy_src, 5.895967611187038e7);
  EXPECT_DOUBLE_EQ(momz_src, -6.98207732332361e12);
  EXPECT_DOUBLE_EQ(energy_src, 6.81241027548468e17);
}

TEST(exampleMMS, euler2D_FAD2)
{
  const int ND = 2;
  double nodalSrc[5];
  const int irho_=0, irhoU_=1, irhoE_=4;

  const double xb=1.234, yb=2.345, zb=0.0, tb = 23.4;

  mms::FAD2_Type x(4, 0, xb);
  mms::FAD2_Type y(4, 1, yb);
  mms::FAD2_Type z(4, 2, zb);
  mms::FAD2_Type t(4, 3, tb);

  x.val() = mms::FAD_Type(4, 0, xb);
  y.val() = mms::FAD_Type(4, 1, yb);
  z.val() = mms::FAD_Type(4, 2, zb);
  t.val() = mms::FAD_Type(4, 3, tb);

  FAD2_Type p = 1.0+(x-t)*y;
  FAD2_Type T = 1.0+cos(1.0+x-t*y);

  FAD2_Type velocity[3];
  velocity[0] = 1.0+sin(x+y-t);
  velocity[1] = 2.0+exp(-x*y-t);
  velocity[2] = 0.0;

  const double gamma = 1.4;
  const double Rgas = 287.0;
  const double Cv = Rgas/(gamma-1.0);

  FAD2_Type rho = p/(Rgas*T);
  FAD2_Type E = rho*(Cv*T);
  for (int i = 0; i < ND; ++i ) {
    E += rho*0.5*velocity[i]*velocity[i];
  }

  // 1) time derivative terms
  nodalSrc[irho_]  =     DT2(rho).val();
  nodalSrc[irhoE_] =     DT2(E).val();
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irhoU_+i] = DT2(rho*velocity[i]).val();
  }

  // 2) advection terms + pressure
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irho_]    +=   D2(  rho*velocity[i], i).val();
    nodalSrc[irhoE_]   +=   D2((E+p)*velocity[i], i).val();

    nodalSrc[irhoU_+i] +=   D2(  p, i).val();

    for (int j = 0; j < ND; ++j ) {
      nodalSrc[irhoU_+i] += D2(  rho*velocity[i]*velocity[j], j).val();
    }
  }

  // compare with output from Mathematica
  EXPECT_NEAR(nodalSrc[irho_],    -77.60303215919267,  1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+0], -11.600229605425282, 1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+1], -177.37206431867156, 1e-10);
  EXPECT_NEAR(nodalSrc[irhoE_],   -416.70415397545025, 1e-10);
}

TEST(exampleMMS, euler3D_FAD2)
{
  const int ND = 3;
  double nodalSrc[5];
  const int irho_=0, irhoU_=1, irhoE_=4;

  const double xb=1.234, yb=2.345, zb=2.3, tb = 2.4;

  mms::FAD2_Type x(4, 0, xb);
  mms::FAD2_Type y(4, 1, yb);
  mms::FAD2_Type z(4, 2, zb);
  mms::FAD2_Type t(4, 3, tb);

  x.val() = mms::FAD_Type(4, 0, xb);
  y.val() = mms::FAD_Type(4, 1, yb);
  z.val() = mms::FAD_Type(4, 2, zb);
  t.val() = mms::FAD_Type(4, 3, tb);

  FAD2_Type p = 1.0+(x-t)*y*z;
  FAD2_Type T = 1.0+cos(z+x-t*y);

  FAD2_Type velocity[3];
  velocity[0] = 1.0+sin(x+y-t*z);
  velocity[1] = 2.0+exp(-x*y-t);
  velocity[2] = 3.0+pow(y-z*t,2);

  const double gamma = 1.4;
  const double Rgas = 287.0;
  const double Cv = Rgas/(gamma-1.0);

  FAD2_Type rho = p/(Rgas*T);
  FAD2_Type E = rho*(Cv*T);
  for (int i = 0; i < ND; ++i ) {
    E += rho*0.5*velocity[i]*velocity[i];
  }

  // 1) time derivative terms
  nodalSrc[irho_]  =     DT2(rho).val();
  nodalSrc[irhoE_] =     DT2(E).val();
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irhoU_+i] = DT2(rho*velocity[i]).val();
  }

  // 2) advection terms + pressure
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irho_]    +=   D2(  rho*velocity[i], i).val();
    nodalSrc[irhoE_]   +=   D2((E+p)*velocity[i], i).val();

    nodalSrc[irhoU_+i] +=   D2(  p, i).val();

    for (int j = 0; j < ND; ++j ) {
      nodalSrc[irhoU_+i] += D2(  rho*velocity[i]*velocity[j], j).val();
    }
  }

  // compare with output from Mathematica
  EXPECT_NEAR(nodalSrc[irho_],    -0.4872484486213013, 1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+0], 4.9391274393182325,  1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+1], -3.658072223766688,  1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+2], -16.518922637184943, 1e-10);
  EXPECT_NEAR(nodalSrc[irhoE_],   -571.1339486153646,  1e-10);
}

TEST(exampleMMS, navierstokes2D_SFAD2)
{
  const int ND = 2;
  double nodalSrc[5];
  const int irho_=0, irhoU_=1, irhoE_=4;

  const double xb=1.234, yb=2.345, zb=0.0, tb = 2.4;

  mms::SFAD2_Type x(4, 0, xb);
  mms::SFAD2_Type y(4, 1, yb);
  mms::SFAD2_Type z(4, 2, zb);
  mms::SFAD2_Type t(4, 3, tb);

  x.val() = mms::SFAD_Type(4, 0, xb);
  y.val() = mms::SFAD_Type(4, 1, yb);
  z.val() = mms::SFAD_Type(4, 2, zb);
  t.val() = mms::SFAD_Type(4, 3, tb);

  SFAD2_Type p = 1.0+(x-t)*y;
  SFAD2_Type T = 1.0+cos(1.0+x-t*y);

  SFAD2_Type velocity[3];
  velocity[0] = 1.0+sin(x+y-t);
  velocity[1] = 2.0+exp(-x*y-t);
  velocity[2] = 0.0;

  const double gamma = 1.4;
  const double Rgas = 287.0;
  const double Cv = Rgas/(gamma-1.0);

  SFAD2_Type rho = p/(Rgas*T);
  SFAD2_Type E = rho*(Cv*T);
  for (int i = 0; i < ND; ++i ) {
    E += rho*0.5*velocity[i]*velocity[i];
  }

  // 1) time derivative terms
  nodalSrc[irho_]  =     DT2(rho).val();
  nodalSrc[irhoE_] =     DT2(E).val();
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irhoU_+i] = DT2(rho*velocity[i]).val();
  }

  // 2) advection terms + pressure
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irho_]    +=   D2(  rho*velocity[i], i).val();
    nodalSrc[irhoE_]   +=   D2((E+p)*velocity[i], i).val();

    nodalSrc[irhoU_+i] +=   D2(  p, i).val();

    for (int j = 0; j < ND; ++j ) {
      nodalSrc[irhoU_+i] += D2(  rho*velocity[i]*velocity[j], j).val();
    }
  }

  // 3.0) compute material props: viscosity and thermal_conductivity
  const double suth_c1 = 0.1;
  const double suth_c2 = 100.0;

  SFAD_Type viscosity = ( suth_c1 * T * sqrt(T)/(T + suth_c2) ).val();

  const double Pr = 0.72;
  const double Cp = gamma*Cv;

  SFAD_Type thermal_conductivity = viscosity*Cp/Pr;

  // 3.1) form stress tensor (no pressure)
  SFAD_Type tau[3][3];
  const double oneThird = 1./3.;

  SFAD_Type div_vel = 0.0;
  for (int i = 0; i < ND; ++i ) {
    div_vel += D2(velocity[i], i);
  }

  for (int i = 0; i < ND; ++i ) {
    tau[i][i] = 2.0*viscosity*(D2(velocity[i], i) - oneThird*div_vel);

    for (int j = 0; j < i; ++j ) {
      tau[i][j] = viscosity*(D2(velocity[i], j) + D2(velocity[j], i));
      tau[j][i] = tau[i][j];
    }
  }

  // 3.2) form total thermal heat flux
  SFAD_Type thermal_flux[3];

  for (int i = 0; i < ND; ++i ) {
    thermal_flux[i] = -thermal_conductivity*D2(T, i);
    for (int j = 0; j < ND; ++j ) {
      thermal_flux[i] -= velocity[j].val()*tau[j][i];
    }
  }

  // 3.3) sum terms into src (using divg)
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irhoE_] += D(thermal_flux[i], i);
    for (int j = 0; j < ND; ++j ) {
      nodalSrc[irhoU_+i] -= D(tau[j][i], j);
    }
  }

  // compare with output from Mathematica
  EXPECT_NEAR(nodalSrc[irho_],    7.774957427548612,  1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+0], 17.092439451251725, 1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+1], 14.43068705973624,  1e-10);
  EXPECT_NEAR(nodalSrc[irhoE_],   28.886475272725438, 1e-10);
}

TEST(exampleMMS, navierstokes3D_FAD2)
{
  const int ND = 3;
  double nodalSrc[5];
  const int irho_=0, irhoU_=1, irhoE_=4;

  const double xb=1.234, yb=2.345, zb=2.3, tb = 2.4;

  mms::SFAD2_Type x(4, 0, xb);
  mms::SFAD2_Type y(4, 1, yb);
  mms::SFAD2_Type z(4, 2, zb);
  mms::SFAD2_Type t(4, 3, tb);

  x.val() = mms::SFAD_Type(4, 0, xb);
  y.val() = mms::SFAD_Type(4, 1, yb);
  z.val() = mms::SFAD_Type(4, 2, zb);
  t.val() = mms::SFAD_Type(4, 3, tb);

  SFAD2_Type p = 1.0+(x-t)*y*z;
  SFAD2_Type T = 1.0+cos(z+x-t*y);

  SFAD2_Type velocity[3];
  velocity[0] = 1.0+sin(x+y-t*z);
  velocity[1] = 2.0+exp(-x*y-t);
  //velocity[2] = 3.0+pow(y-z*t,2);
  velocity[2] = 3.0+(y-z*t)*(y-z*t);

  const double gamma = 1.4;
  const double Rgas = 287.0;
  const double Cv = Rgas/(gamma-1.0);

  SFAD2_Type rho = p/(Rgas*T);
  SFAD2_Type E = rho*(Cv*T);
  for (int i = 0; i < ND; ++i ) {
    E += rho*0.5*velocity[i]*velocity[i];
  }

  // 1) time derivative terms
  nodalSrc[irho_]  =     DT2(rho).val();
  nodalSrc[irhoE_] =     DT2(E).val();
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irhoU_+i] = DT2(rho*velocity[i]).val();
  }

  // 2) advection terms + pressure
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irho_]    +=   D2(  rho*velocity[i], i).val();
    nodalSrc[irhoE_]   +=   D2((E+p)*velocity[i], i).val();

    nodalSrc[irhoU_+i] +=   D2(  p, i).val();

    for (int j = 0; j < ND; ++j ) {
      nodalSrc[irhoU_+i] += D2(  rho*velocity[i]*velocity[j], j).val();
    }
  }

  // 3.0) compute material props: viscosity and thermal_conductivity
  const double suth_c1 = 0.1;
  const double suth_c2 = 100.0;

  SFAD_Type viscosity = ( suth_c1 * T * sqrt(T)/(T + suth_c2) ).val();

  const double Pr = 0.72;
  const double Cp = gamma*Cv;

  SFAD_Type thermal_conductivity = viscosity*Cp/Pr;

  // 3.1) form stress tensor (no pressure)
  SFAD_Type tau[3][3];
  const double oneThird = 1./3.;

  SFAD_Type div_vel = 0.0;
  for (int i = 0; i < ND; ++i ) {
    div_vel += D2(velocity[i], i);
  }

  for (int i = 0; i < ND; ++i ) {
    tau[i][i] = 2.0*viscosity*(D2(velocity[i], i) - oneThird*div_vel);

    for (int j = 0; j < i; ++j ) {
      tau[i][j] = viscosity*(D2(velocity[i], j) + D2(velocity[j], i));
      tau[j][i] = tau[i][j];
    }
  }

  // 3.2) form total thermal heat flux
  SFAD_Type thermal_flux[3];

  for (int i = 0; i < ND; ++i ) {
    thermal_flux[i] = -thermal_conductivity*D2(T, i);
    for (int j = 0; j < ND; ++j ) {
      thermal_flux[i] -= velocity[j].val()*tau[j][i];
    }
  }

  // 3.3) sum terms into src (using divg)
  for (int i = 0; i < ND; ++i ) {
    nodalSrc[irhoE_] += D(thermal_flux[i], i);
    for (int j = 0; j < ND; ++j ) {
      nodalSrc[irhoU_+i] -= D(tau[j][i], j);
    }
  }

  // compare with output from Mathematica
  EXPECT_NEAR(nodalSrc[irho_],    -0.4872484486213013, 1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+0], 4.944557100416018,   1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+1], -3.673218932569313,  1e-10);
  EXPECT_NEAR(nodalSrc[irhoU_+2], -16.558200406094826, 1e-10);
  EXPECT_NEAR(nodalSrc[irhoE_],   -582.2566990445401,  1e-10);
}

TEST(exampleMMS, solidMechanics3D)
{
  // NOTE: some of the code below should be moved into library calls

  // point in ref config = unit cube
  const double Xb=0.8124, Yb=0.763, Zb=0.659, tb = 1.0;

  FAD2_Type X(4,0,FAD_Type(4,0,Xb));
  FAD2_Type Y(4,1,FAD_Type(4,0,Yb));
  FAD2_Type Z(4,2,FAD_Type(4,0,Zb));
  FAD2_Type t(4,3,FAD_Type(4,0,tb));

  // displacement function - user specified
  FAD2_Type phi[3];
  phi[0] = X + 0.2*sin(0.5*M_PI*t)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
  phi[1] = Y + 0.3*sin(0.5*M_PI*t)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
  phi[2] = Z + 0.4*sin(0.5*M_PI*t)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);

  // deformation gradient - library call
  FAD_Type F[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      F[i][j] = phi[i].dx(j);
    }
  }

  // determinant - library call
  FAD_Type J =
     F[0][0]*(F[1][1]*F[2][2]-F[1][2]*F[2][1])
    -F[0][1]*(F[1][0]*F[2][2]-F[1][2]*F[2][0])
    +F[0][2]*(F[1][0]*F[2][1]-F[1][1]*F[2][0]);

  // left Cauchy-Green tensor (b = F*F^t) - library call
  FAD_Type b[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b[i][j] = 0;
      for (int k=0; k<3; k++) {
	b[i][j] += F[i][k]*F[j][k];
      }
    }
  }

  const double mu  = 1.0;
  const double lambda = 1.0;
  const double bulk = lambda + mu*(2.0/3.0);

  // Cauchy stress tensor - material dependent.
  // Here we use neo-Hookean
  FAD_Type sigma[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      sigma[i][j] = 0;
    }
  }

  FAD_Type pressure = 0.5*bulk*(J-1.0/J);

  FAD_Type factor = pow(J, -5.0/3.0);
  FAD_Type trace_b = b[0][0]+b[1][1]+b[2][2];

  for (int i=0; i<3; i++) {
    sigma[i][i] += pressure - mu*factor*trace_b/3.0;
    for (int j=0; j<3; j++) {
      sigma[i][j] += mu*factor*b[i][j];
    }
  }

  // inverse deformation gradient - library call
  FAD_Type inv_F[3][3];

  inv_F[0][0] =  (F[2][2]*F[1][1]-F[2][1]*F[1][2])/J;
  inv_F[1][0] = -(F[2][2]*F[1][0]-F[2][0]*F[1][2])/J;
  inv_F[2][0] =  (F[2][1]*F[1][0]-F[2][0]*F[1][1])/J;

  inv_F[0][1] = -(F[2][2]*F[0][1]-F[2][1]*F[0][2])/J;
  inv_F[1][1] =  (F[2][2]*F[0][0]-F[2][0]*F[0][2])/J;
  inv_F[2][1] = -(F[2][1]*F[0][0]-F[2][0]*F[0][1])/J;

  inv_F[0][2] =  (F[1][2]*F[0][1]-F[1][1]*F[0][2])/J;
  inv_F[1][2] = -(F[1][2]*F[0][0]-F[1][0]*F[0][2])/J;
  inv_F[2][2] =  (F[1][1]*F[0][0]-F[1][0]*F[0][1])/J;

  // source term = -div(sigma) + rho*ddot(u) ...
  double mom_src[3];

  const double rho = 1.0;

  // compute div(sigma) in current coordinates using chain rule

  // this is because we compute everything based on a point
  // in the original configuration, but need the divergence
  // w.r.t. the current configuration

  for (int i=0; i<3; i++) {
    mom_src[i] = 0;
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	// mom_i = - D_x_j sigma_ij = - D_X_k sigma_ij D_X_k/D_x_j
	mom_src[i] -= (sigma[i][j]).dx(k)*inv_F[k][j].val();
      }
    }
    mom_src[i] += rho*phi[i].dx(_T).dx(_T);
  }

  // TODO: test against output from Mathematica

  // DEBUG
  /*
  std::cout << "Time: " << t.val().val() << std::endl;
  std::cout << "Initial point: "
	    << X.val().val() << " "
	    << Y.val().val() << " "
	    << Z.val().val() << std::endl;
  std::cout << "Deformed point: "
	    << phi[0].val().val() << " "
	    << phi[1].val().val() << " "
	    << phi[2].val().val() << std::endl;
  std::cout << "Deformation gradient: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << F[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Inverse deformation gradient: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << inv_F[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Determinant: " << J.val() << std::endl;
  std::cout << "Left Cauchy-Green tensor: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << b[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Pressure: " << pressure.val() << std::endl;
  std::cout << "Cauchy stress tensor: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << sigma[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "momentum source: " << std::endl;
  for (int i=0; i<3; i++) {
    std::cout << mom_src[i] << " ";
  }
  std::cout << std::endl;
  */
}

TEST(exampleMMS, viscousNavierStokes2D)
{
  const int ND = 2;

  const double gamma = 1.4;
  const double Rgas = 287.097384766765;

  const double suth_c1 = 0.01458;
  const double suth_c2 = 110.4;
  const double Pr = 0.72;

  const double Cv = Rgas/(gamma-1.0);

  const double M0 = 2.5;
  const double p0 = 35651.28116;
  const double T0 = 236.215;

  const double rho0 = p0/(Rgas*T0);
  const double u0 = M0*sqrt(gamma*p0/rho0);

  const double eps = 0.05;
  const double Cp = gamma*Cv;

  const double xb=1.34, yb=0.356, tb = 1.0;
  SFAD2_Type x(4, 0, xb), y(4, 1, yb), t(4, 3, tb);

  x.val() = SFAD_Type(4, 0, xb);
  y.val() = SFAD_Type(4, 1, yb);
  t.val() = SFAD_Type(4, 3, tb);

  SFAD2_Type u = u0*(1.0+eps*pow(sin(0.5*M_PI*x),2)*pow(sin(M_PI*y),2));
  SFAD2_Type v = u0*eps*pow(sin(0.5*M_PI*x),2)*pow(sin(M_PI*y),2);
  SFAD2_Type p = p0*(1.0 + eps*sin(0.5*M_PI*x)*cos(M_PI*y));
  SFAD2_Type T = T0*(1.0 + eps*sin(0.25*M_PI*x)*cos(M_PI*y));

  SFAD2_Type velocity[3];
  velocity[0] = u;
  velocity[1] = v;
  velocity[2] = 0.0;

  SFAD_Type viscosity = ( suth_c1 * T * sqrt(T)/(T + suth_c2) ).val();

  SFAD_Type thermal_conductivity = viscosity*Cp/Pr;

  // form stress tensor (no pressure)
  SFAD_Type tau[3][3];
  const double oneThird = 1./3.;

  SFAD_Type div_vel = 0.0;
  for (int i = 0; i < ND; ++i ) {
    div_vel += D2(velocity[i], i);
  }

  for (int i = 0; i < ND; ++i ) {
    tau[i][i] = 2.0*viscosity*(D2(velocity[i], i) - oneThird*div_vel);

    for (int j = 0; j < i; ++j ) {
      tau[i][j] = viscosity*(D2(velocity[i], j) + D2(velocity[j], i));
      tau[j][i] = tau[i][j];
    }
  }

  // form total thermal heat flux
  SFAD_Type thermal_flux[3];

  for (int i = 0; i < ND; ++i ) {
    thermal_flux[i] = -thermal_conductivity*D2(T, i);
  }

  for (int i = 0; i < ND; ++i ) {
    for (int j = 0; j < ND; ++j ) {
      thermal_flux[i] -= velocity[j].val()*tau[j][i];
    }
  }

  double mom_src[3];
  for (int i = 0; i < ND; ++i ) {
    mom_src[i] = 0.0;
  }

  double energy_src = 0.0;

  // sum terms into src (using divg)
  for (int i = 0; i < ND; ++i ) {
    energy_src += D(thermal_flux[i], i);
    for (int j = 0; j < ND; ++j ) {
      mom_src[i] -= D(tau[j][i], j);
    }
  }

  // test against output from Mathematica
  EXPECT_NEAR(mom_src[0], 76.58413133799213, 1e-10);
  EXPECT_NEAR(mom_src[1], 92.05081838663132, 1e-10);
  EXPECT_NEAR(energy_src, 70295.4770236427,  1e-10);
}

TEST(exampleMMS, FADbug)
{
  const double xb=1.234, yb=2.345, zb=2.3, tb = 2.4;

  mms::FAD2_Type x(4, 0, xb);
  mms::FAD2_Type y(4, 1, yb);
  mms::FAD2_Type z(4, 2, zb);
  mms::FAD2_Type t(4, 3, tb);

  x.val() = mms::FAD_Type(4, 0, xb);
  y.val() = mms::FAD_Type(4, 1, yb);
  z.val() = mms::FAD_Type(4, 2, zb);
  t.val() = mms::FAD_Type(4, 3, tb);

  mms::FAD2_Type w = pow(y-z*t,2);

  std::cout << "w = "
	    << w.val() << " "
	    << w.dx(0) << " "
	    << w.dx(1) << " "
	    << w.dx(2) << " "
	    << w.dx(3) << std::endl;

}

TEST(performanceMMS, linearElasticity)
{
  //CALLGRIND_START_INSTRUMENTATION;

  const int num_points = 100; //1e6;

  const double tol = 1e-10;
    
  const bool use_mms = false;

  const double rho = 1349.0;  
  const double mu  = 1123.0;
  const double lam = 2342.0;
      
  // source term = -div(sigma) + rho*ddot(u) ...
  double mom_src[3];
      
  //CALLGRIND_TOGGLE_COLLECT;

  for (int ip=0; ip<num_points; ip++) {

    // point in ref config = unit cube
    const double x=((double)rand()/(double)RAND_MAX);
    const double y=((double)rand()/(double)RAND_MAX);
    const double z=((double)rand()/(double)RAND_MAX);
    const double t=((double)rand()/(double)RAND_MAX);

    if (use_mms) {

      FAD2_Type X(4,0,FAD_Type(4,0,x));
      FAD2_Type Y(4,1,FAD_Type(4,1,y));
      FAD2_Type Z(4,2,FAD_Type(4,2,z));
      FAD2_Type T(4,3,FAD_Type(4,3,t));
      
      // displacement function - user specified
      FAD2_Type phi[3];
      phi[0] = X + 0.2*sin(0.5*M_PI*T)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
      phi[1] = Y + 0.3*sin(0.5*M_PI*T)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
      phi[2] = Z + 0.4*sin(0.5*M_PI*T)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
      
      // strain tensor
      FAD_Type strain[3][3];
      FAD_Type trace_strain = 0;
      for (int i=0; i<3; i++) {
	trace_strain += phi[i].dx(i);
	for (int j=0; j<3; j++) {
	  strain[i][j] = 0.5*(phi[i].dx(j) + phi[j].dx(i));
	}
      }
      
      // Cauchy stress tensor - material dependent.
      // Here we use neo-Hookean
      FAD_Type sigma[3][3];
      for (int i=0; i<3; i++) {
	sigma[i][i] += lam * trace_strain;
	for (int j=0; j<3; j++) {
	  sigma[i][j] += 2 * mu * strain[i][j];
	}
      }
      
      for (int i=0; i<3; i++) {
	mom_src[i] = rho*phi[i].dx(_T).dx(_T);
	for (int j=0; j<3; j++) {      
	  mom_src[i] += -sigma[i][j].dx(j);
	}
      }
    
    // DEBUG
    /*
    std::cout << "Time: " << T.val().val() << std::endl;
    std::cout << "Initial point: "
	      << X.val().val() << " "
	      << Y.val().val() << " "
	      << Z.val().val() << std::endl;
    std::cout << "Deformed point: "
	      << phi[0].val().val() << " "
	      << phi[1].val().val() << " "
	      << phi[2].val().val() << std::endl;
    std::cout << "Cauchy stress tensor: " << std::endl;
    for (int i=0; i<3; i++) {
      for (int j=0; j<3; j++) {
	std::cout << sigma[i][j].val() << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "momentum source: " << std::endl;
    for (int i=0; i<3; i++) {
      std::cout << mom_src[i] << " ";
    }
    std::cout << std::endl;
    */

    // compare with SymPy

    EXPECT_NEAR(mom_src[0], -0.05*pow(M_PI, 2)*rho*sin(0.5*M_PI*t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) - pow(M_PI, 2)*(-0.2*lam*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*lam*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z) + 0.3*lam*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y) - 0.8*mu*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*mu*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z) + 0.3*mu*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y))*sin(0.5*M_PI*t), tol);

    EXPECT_NEAR(mom_src[1], -0.075*pow(M_PI, 2)*rho*sin(0.5*M_PI*t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) - pow(M_PI, 2)*(-0.3*lam*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*lam*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*lam*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y) - 1.2*mu*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*mu*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*mu*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y))*sin(0.5*M_PI*t), tol);

    EXPECT_NEAR(mom_src[2], -0.1*pow(M_PI, 2)*rho*sin(0.5*M_PI*t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) - pow(M_PI, 2)*(-0.4*lam*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.3*lam*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*lam*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z) - 1.6*mu*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.3*mu*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*mu*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z))*sin(0.5*M_PI*t), tol);
    
    }
    else { // use exact src expression

      mom_src[0] = -0.05*pow(M_PI, 2)*rho*sin(0.5*M_PI*t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) - pow(M_PI, 2)*(-0.2*lam*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*lam*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z) + 0.3*lam*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y) - 0.8*mu*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*mu*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z) + 0.3*mu*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y))*sin(0.5*M_PI*t);
      
      mom_src[1] = -0.075*pow(M_PI, 2)*rho*sin(0.5*M_PI*t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) - pow(M_PI, 2)*(-0.3*lam*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*lam*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*lam*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y) - 1.2*mu*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.4*mu*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*mu*sin(M_PI*z)*cos(M_PI*x)*cos(M_PI*y))*sin(0.5*M_PI*t);
      
      mom_src[2] = -0.1*pow(M_PI, 2)*rho*sin(0.5*M_PI*t)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) - pow(M_PI, 2)*(-0.4*lam*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.3*lam*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*lam*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z) - 1.6*mu*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z) + 0.3*mu*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z) + 0.2*mu*sin(M_PI*y)*cos(M_PI*x)*cos(M_PI*z))*sin(0.5*M_PI*t);
    }

  }

  //CALLGRIND_TOGGLE_COLLECT;

  //CALLGRIND_STOP_INSTRUMENTATION;
}

const double gamma = 1.4;
const double Rgas = 287;
const double eps(0.05);
const double Cv(Rgas/(gamma-1.0));
const double M0(2.5);
const double p0(35651.28116);
const double T0(236.215);
const double rho0(p0/(Rgas*T0));
const double u0(M0*sqrt(gamma*p0/rho0));

const double suth_c1 = 1.458e-1;
const double suth_c2 = 110.4;
const double Pr = 0.72;
const double Cp = gamma*Cv;

template <typename ScalarT>
ScalarT
pressure(const ScalarT & x,const ScalarT & y,const ScalarT & z,const ScalarT & t)
  {return p0*(1.0 - eps*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z/1.5)*sin(M_PI*(t*(t/eps-x-y-z))));}

template <typename ScalarT>
ScalarT
temperature(const ScalarT & x,const ScalarT & y,const ScalarT & z,const ScalarT & t)
  {return T0*(1.0 + eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z/1.5)*sin(M_PI*(t*(t/eps-x-y-z))));}

template <typename ScalarT>
ScalarT
xvelocity(const ScalarT & x,const ScalarT & y,const ScalarT & z,const ScalarT & t)
  {return u0*(1.0 - eps*sin(M_PI*x)*cos(M_PI*y)*cos(M_PI*z/1.5)*sin(M_PI*(t*(t/eps-x-y-z))));}

template <typename ScalarT>
ScalarT
yvelocity(const ScalarT & x,const ScalarT & y,const ScalarT & z,const ScalarT & t)
  {return u0*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z/1.5)*sin(M_PI*(t*(t/eps-x-y-z)));}

template <typename ScalarT>
ScalarT
zvelocity(const ScalarT & x,const ScalarT & y,const ScalarT & z,const ScalarT & t)
  {return u0*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z/1.5)*sin(M_PI*(t*(t/eps-x-y-z)));}

double mass_src(const double & x,const double & y,const double & z,const double & t)
{return -2*M_PI*eps*p0*t*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.666666666666667*M_PI*eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - eps*p0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - eps*p0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 2*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*u0*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2));}

double xmom_src(const double & x,const double & y,const double & z,const double & t)
{return -2*M_PI*eps*p0*t*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.666666666666667*M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 2*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 1.33333333333333*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - 2*eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 3*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2));}

double ymom_src(const double & x,const double & y,const double & z,const double & t)
{return -4*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.33333333333333*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*sin(0.666666666666667*M_PI*z)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - M_PI*eps*p0*t*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 2*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0));}

double zmom_src(const double & x,const double & y,const double & z,const double & t)
{return -4*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.33333333333333*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*sin(0.666666666666667*M_PI*z)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - M_PI*eps*p0*t*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 2*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0));}

double energy_src(const double & x,const double & y,const double & z,const double & t)
{return -Cv*eps*p0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/Rgas - 2*M_PI*eps*t*u0*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)))*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)))*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) + eps*u0*(Cv*p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))/Rgas - 2.0*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2.0*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 3*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + eps*u0*(Cv*p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/Rgas - 2.0*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.33333333333333*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*sin(0.666666666666667*M_PI*z)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 2.0*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(Cv*p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/Rgas - 2.0*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2.0*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*x)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 3*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2))) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))) - 1.0*pow(eps, 3)*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(sin(M_PI*x), 3)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - 1.0*pow(eps, 3)*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 3)*pow(sin(M_PI*y), 3)*pow(sin(0.666666666666667*M_PI*z), 3)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 2.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - 1.5*eps*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - 0.5*eps*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2));}

TEST(performanceMMS, Euler3DTransient)
{
  const int num_points = 100; //1e6;

  const double tol = 3e-6;

  const bool use_mms = true;

  const int ND = 3;
  double nodalSrc[5];
  const int irho_=0, irhoU_=1, irhoE_=4;

  for (int ip=0; ip<num_points; ip++) {

    // point in ref config = unit cube
    const double x=((double)rand()/(double)RAND_MAX);
    const double y=((double)rand()/(double)RAND_MAX);
    const double z=((double)rand()/(double)RAND_MAX);
    const double t=((double)rand()/(double)RAND_MAX);

    if (use_mms) {
    	mms::SFAD_Type X(4, 0, x);
    	mms::SFAD_Type Y(4, 1, y);
    	mms::SFAD_Type Z(4, 2, z);
    	mms::SFAD_Type T(4, 3, t);

    	SFAD_Type p = pressure<SFAD_Type>(X,Y,Z,T);
    	SFAD_Type Temp = temperature<SFAD_Type>(X,Y,Z,T);

    	SFAD_Type velocity[3];
    	velocity[0] = xvelocity<SFAD_Type>(X,Y,Z,T);
    	velocity[1] = yvelocity<SFAD_Type>(X,Y,Z,T);
    	velocity[2] = zvelocity<SFAD_Type>(X,Y,Z,T);

        SFAD_Type rho = p/(Rgas*Temp);
        SFAD_Type E = rho*(Cv*Temp);
        for (int i = 0; i < ND; ++i ) {
          E += rho*0.5*velocity[i]*velocity[i];
        }

        // eval source terms common to Euler and NavierStokes
        //   NOTE last .val() is to reduce Fad_Type to scalar

        // 1) time derivative terms
        nodalSrc[irho_]  =     rho.dx(3);
        nodalSrc[irhoE_] =     E.dx(3);
        for (int i = 0; i < ND; ++i ) {
          nodalSrc[irhoU_+i] = (rho*velocity[i]).dx(3);
        }

        // 2) advection terms + pressure
        for (int i = 0; i < ND; ++i ) {
          nodalSrc[irho_]    +=   (  rho*velocity[i]).dx(i);
          nodalSrc[irhoE_]   +=   ((E+p)*velocity[i]).dx(i);

          nodalSrc[irhoU_+i] +=   p.dx(i);

          for (int j = 0; j < ND; ++j ) {
            nodalSrc[irhoU_+i] += (  rho*velocity[i]*velocity[j]).dx(j);
          }
        }

	// compare with output from SymPy

	EXPECT_NEAR(nodalSrc[irho_],    mass_src(x,y,z,t), tol);
	EXPECT_NEAR(nodalSrc[irhoU_+0], xmom_src(x,y,z,t), tol);
	EXPECT_NEAR(nodalSrc[irhoU_+1], ymom_src(x,y,z,t), tol);
	EXPECT_NEAR(nodalSrc[irhoU_+2], zmom_src(x,y,z,t), tol);
	EXPECT_NEAR(nodalSrc[irhoE_], energy_src(x,y,z,t), tol);

    }
    else { // use exact src expression

      nodalSrc[irho_] =    mass_src(x,y,z,t);
      nodalSrc[irhoU_+0] = xmom_src(x,y,z,t);
      nodalSrc[irhoU_+1] = ymom_src(x,y,z,t);
      nodalSrc[irhoU_+2] = zmom_src(x,y,z,t);
      nodalSrc[irhoE_] = energy_src(x,y,z,t);

    }

  }

  //CALLGRIND_TOGGLE_COLLECT;

  //CALLGRIND_STOP_INSTRUMENTATION;
}

double ns_mass_src(const double & x,const double & y,const double & z,const double & t)
{return -2*M_PI*eps*p0*t*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.666666666666667*M_PI*eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - eps*p0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - eps*p0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + p0*u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 2*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*u0*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2));}

double ns_xmom_src(const double & x,const double & y,const double & z,const double & t)
{return -2*M_PI*eps*p0*t*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.666666666666667*M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + 2.0*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(0.666666666666667*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) + 0.666666666666667*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) - T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 1.33333333333333*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*u0*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 2*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(0.666666666666667*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.222222222222222*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.333333333333333*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*pow(M_PI, 2)*eps*u0*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z) - 0.333333333333333*pow(M_PI, 2)*eps*u0*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y) + 0.666666666666667*u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 2*pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3.0*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(0.666666666666667*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) + 0.666666666666667*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 2*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 1.33333333333333*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - 2*eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 3*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2));}

double ns_ymom_src(const double & x,const double & y,const double & z,const double & t)
{return -4*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.33333333333333*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*sin(0.666666666666667*M_PI*z)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - M_PI*eps*p0*t*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + 2.0*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) - T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2*pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-2*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2.0*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.444444444444444*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.222222222222222*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 1.0*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.333333333333333*u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3.0*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 2*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0));}

double ns_zmom_src(const double & x,const double & y,const double & z,const double & t)
{return -4*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.33333333333333*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*sin(0.666666666666667*M_PI*z)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - M_PI*eps*p0*t*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + M_PI*eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2.0*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) - T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2*pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.666666666666667*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-2*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 3*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.333333333333333*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.296296296296296*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.333333333333333*u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.666666666666667*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3.0*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3*T0*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - pow(eps, 2)*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(2*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 2*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + eps*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + eps*p0*u0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0));}

double ns_energy_src(const double & x,const double & y,const double & z,const double & t)
{return Cp*pow(T0, 3)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*pow(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z), 2)/(Pr*pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2)) + Cp*pow(T0, 3)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*pow(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y), 2)/(Pr*pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2)) + Cp*pow(T0, 3)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*pow(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x), 2)/(Pr*pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2)) - Cp*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 1.33333333333333*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.444444444444444*pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)))/(Pr*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - Cp*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)))/(Pr*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - Cp*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2*pow(M_PI, 2)*eps*t*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)))/(Pr*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3*Cp*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*pow(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z), 2)/(2*Pr*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3*Cp*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*pow(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y), 2)/(2*Pr*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3*Cp*pow(T0, 2)*suth_c1*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*pow(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x), 2)/(2*Pr*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - Cv*eps*p0*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/Rgas + 2*M_PI*T0*eps*suth_c1*t*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) + M_PI*T0*eps*suth_c1*t*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) + M_PI*T0*eps*suth_c1*t*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) + 2.0*M_PI*T0*eps*suth_c1*t*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) + 2.0*M_PI*T0*eps*suth_c1*t*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 0.666666666666667*M_PI*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - M_PI*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - M_PI*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - M_PI*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*M_PI*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 1.33333333333333*M_PI*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2*M_PI*eps*t*u0*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)))*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)))*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) + pow(T0, 2)*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + 2.0*pow(T0, 2)*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + 2.0*pow(T0, 2)*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + pow(T0, 2)*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) + 2.0*pow(T0, 2)*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(0.666666666666667*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) + 0.666666666666667*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/pow(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2, 2) - T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2*pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2*pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.666666666666667*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-2*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 2.0*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.444444444444444*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-2*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 3*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.333333333333333*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.296296296296296*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.333333333333333*u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.666666666666667*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-0.333333333333333*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.222222222222222*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 1.0*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*pow(M_PI, 2)*eps*u0*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 0.333333333333333*u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3.0*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-2*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3.0*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-0.333333333333333*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) + 0.666666666666667*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) - 0.333333333333333*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3*T0*eps*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*u0*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 1.33333333333333*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.444444444444444*pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) - pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*u0*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y) + u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) - 2*pow(M_PI, 2)*eps*t*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(0.666666666666667*pow(M_PI, 2)*eps*pow(t, 2)*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 0.222222222222222*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.333333333333333*pow(M_PI, 2)*eps*t*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*cos(M_PI*y)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*pow(M_PI, 2)*eps*t*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*x)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*pow(M_PI, 2)*eps*u0*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(0.666666666666667*M_PI*z) - 0.333333333333333*pow(M_PI, 2)*eps*u0*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y) + 0.666666666666667*u0*(pow(M_PI, 2)*eps*pow(t, 2)*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 2*pow(M_PI, 2)*eps*t*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + pow(M_PI, 2)*eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 2.0*T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*(0.666666666666667*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) + 0.666666666666667*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) - 3*T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3*T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*(-M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*u0*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)))/(2*(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2)) - 3.0*T0*suth_c1*u0*sqrt(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(-M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*(0.666666666666667*M_PI*eps*t*u0*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.222222222222222*M_PI*eps*u0*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z) - 0.333333333333333*M_PI*eps*u0*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y) + 0.666666666666667*u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)))/(T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0) + suth_c2) + eps*u0*(Cv*p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))/Rgas - 2.0*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2.0*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*sin(M_PI*y)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 3*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + eps*u0*(Cv*p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/Rgas - 2.0*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.33333333333333*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*sin(0.666666666666667*M_PI*z)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(0.666666666666667*M_PI*z)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) + 2.0*M_PI*eps*sin(M_PI*x)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 0.666666666666667*M_PI*eps*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*t*(-x - y - z + t/eps))*cos(0.666666666666667*M_PI*z))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)))*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + u0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(Cv*p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/Rgas - 2.0*M_PI*pow(eps, 2)*p0*t*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 2.0*M_PI*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*sin(M_PI*x)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*x)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + p0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*(3*M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - 3*M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*(M_PI*eps*t*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2))) + u0*(M_PI*eps*t*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps)) - M_PI*eps*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z))*(Cv*p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)/Rgas + p0*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0) + 1.0*pow(eps, 2)*p0*pow(u0, 2)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) + 0.5*p0*pow(u0, 2)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0))) - 1.0*pow(eps, 3)*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(sin(M_PI*x), 3)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - 1.0*pow(eps, 3)*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 3)*pow(sin(M_PI*y), 3)*pow(sin(0.666666666666667*M_PI*z), 3)*pow(sin(M_PI*t*(-x - y - z + t/eps)), 2)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2)) + 2.0*pow(eps, 2)*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0)*pow(sin(M_PI*x), 2)*pow(sin(M_PI*y), 2)*pow(sin(0.666666666666667*M_PI*z), 2)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - 1.5*eps*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 2)*sin(M_PI*x)*cos(M_PI*y)*cos(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0)) - 0.5*eps*p0*pow(u0, 2)*(M_PI*(-x - y - z + t/eps) + M_PI*t/eps)*pow(-eps*sin(M_PI*x)*sin(M_PI*t*(-x - y - z + t/eps))*cos(M_PI*y)*cos(0.666666666666667*M_PI*z) + 1.0, 3)*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*cos(M_PI*t*(-x - y - z + t/eps))/(Rgas*T0*pow(eps*sin(M_PI*x)*sin(M_PI*y)*sin(0.666666666666667*M_PI*z)*sin(M_PI*t*(-x - y - z + t/eps)) + 1.0, 2));}

TEST(performanceMMS, NavierStokes3DTransient)
{
  const int num_points = 100; //1e7;

  const double tol = 2e-6;

  const bool use_mms = true;

  const int ND = 3;
  double nodalSrc[5];
  const int irho_=0, irhoU_=1, irhoE_=4;

  for (int ip=0; ip<num_points; ip++) {

    // point in ref config = unit cube
    const double x=((double)rand()/(double)RAND_MAX);
    const double y=((double)rand()/(double)RAND_MAX);
    const double z=((double)rand()/(double)RAND_MAX);
    const double t=((double)rand()/(double)RAND_MAX);

    if (use_mms) {
      mms::SFAD2_Type X(4, 0, x);
      mms::SFAD2_Type Y(4, 1, y);
      mms::SFAD2_Type Z(4, 2, z);
      mms::SFAD2_Type T(4, 3, t);

      X.val() = mms::SFAD_Type(4, 0, x);
      Y.val() = mms::SFAD_Type(4, 1, y);
      Z.val() = mms::SFAD_Type(4, 2, z);
      T.val() = mms::SFAD_Type(4, 3, t);

      SFAD2_Type p = pressure<SFAD2_Type>(X,Y,Z,T);
      SFAD2_Type Temp = temperature<SFAD2_Type>(X,Y,Z,T);
      
      SFAD2_Type velocity[3];
      velocity[0] = xvelocity<SFAD2_Type>(X,Y,Z,T);
      velocity[1] = yvelocity<SFAD2_Type>(X,Y,Z,T);
      velocity[2] = zvelocity<SFAD2_Type>(X,Y,Z,T);
      
      SFAD2_Type rho = p/(Rgas*Temp);
      SFAD2_Type E = rho*(Cv*Temp);
      for (int i = 0; i < ND; ++i ) {
	E += rho*0.5*velocity[i]*velocity[i];
      }
      
      // eval source terms common to Euler and NavierStokes
      //   NOTE last .val() is to reduce Fad_Type to scalar
      
      // 1) time derivative terms
      nodalSrc[irho_]  =     rho.dx(3).val();
      nodalSrc[irhoE_] =     E.dx(3).val();
      for (int i = 0; i < ND; ++i ) {
	nodalSrc[irhoU_+i] = (rho*velocity[i]).dx(3).val();
      }
      
      // 2) advection terms + pressure
      for (int i = 0; i < ND; ++i ) {
	nodalSrc[irho_]    +=   (  rho*velocity[i]).dx(i).val();
	nodalSrc[irhoE_]   +=   ((E+p)*velocity[i]).dx(i).val();
	
	nodalSrc[irhoU_+i] +=   p.dx(i).val();
	
	for (int j = 0; j < ND; ++j ) {
	  nodalSrc[irhoU_+i] += (  rho*velocity[i]*velocity[j]).dx(j).val();
	}
      }
      
      // 3) compute the viscous flux terms
      
      // 3.0) compute material props: viscosity and thermal_conductivity
      SFAD_Type viscosity = ( suth_c1 * Temp * sqrt(Temp)/(Temp + suth_c2) ).val();
      
      SFAD_Type thermal_conductivity = viscosity*Cp/Pr;
      
      // 3.1) form stress tensor (no pressure)
      SFAD_Type tau[3][3];
      const double oneThird = 1./3.;
      
      SFAD_Type div_vel = 0.0;
      for (int i = 0; i < ND; ++i ) {
	div_vel += velocity[i].dx(i);
      }
      
      for (int i = 0; i < ND; ++i ) {
	tau[i][i] = 2.0*viscosity*(velocity[i].dx(i) - oneThird*div_vel);
	
	for (int j = 0; j < i; ++j ) {
	  tau[i][j] = viscosity*(velocity[i].dx(j) + velocity[j].dx(i));
	  tau[j][i] = tau[i][j];
	}
      }
      
      // 3.2) form total thermal heat flux
      SFAD_Type thermal_flux[3];
      
      for (int i = 0; i < ND; ++i ) {
	thermal_flux[i] = -thermal_conductivity*Temp.dx(i);
	for (int j = 0; j < ND; ++j ) {
	  thermal_flux[i] -= velocity[j].val()*tau[j][i];
	}
      }
      
      // 3.3) sum terms into src (using divg)
      for (int i = 0; i < ND; ++i ) {
	nodalSrc[irhoE_] += thermal_flux[i].dx(i);
	for (int j = 0; j < ND; ++j ) {
	  nodalSrc[irhoU_+i] -= tau[j][i].dx(j);
	}
      }
      
      // compare with output from SymPy

      EXPECT_NEAR(nodalSrc[irho_],    ns_mass_src(x,y,z,t), tol);
      EXPECT_NEAR(nodalSrc[irhoU_+0], ns_xmom_src(x,y,z,t), tol);
      EXPECT_NEAR(nodalSrc[irhoU_+1], ns_ymom_src(x,y,z,t), tol);
      EXPECT_NEAR(nodalSrc[irhoU_+2], ns_zmom_src(x,y,z,t), tol);
      EXPECT_NEAR(nodalSrc[irhoE_], ns_energy_src(x,y,z,t), tol);

    }
    else { // use exact src expression

      nodalSrc[irho_] =    ns_mass_src(x,y,z,t);
      nodalSrc[irhoU_+0] = ns_xmom_src(x,y,z,t);
      nodalSrc[irhoU_+1] = ns_ymom_src(x,y,z,t);
      nodalSrc[irhoU_+2] = ns_zmom_src(x,y,z,t);
      nodalSrc[irhoE_] = ns_energy_src(x,y,z,t);

    }

  }

  //CALLGRIND_TOGGLE_COLLECT;

  //CALLGRIND_STOP_INSTRUMENTATION;
}

} // namespace unit_tests
} // namespace mms

