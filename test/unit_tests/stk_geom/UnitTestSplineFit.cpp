// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/Util.hpp>
#include <percept/mesh/geometry/stk_geom/SplineFit.hpp>
#include <percept/mesh/geometry/stk_geom/LocalCubicSplineFit.hpp>
#include <gtest/gtest.h>

#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

  namespace geom
  {
    namespace unit_tests
    {


      //=============================================================================
      //=============================================================================
      //=============================================================================

#if !defined(NO_GEOM_SUPPORT)
      TEST(unit_stk_geom, test_1)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        /// fit a quadratic
        LocalCubicSplineFit cf;
        int n = 11;
        Vectors2D Q(n);
        for (int i=0; i < n; i++)
          {
            double x = double(i)/double(n-1);
            double y = x*x;
            Q[i] = Vector2D(x,y);
          }
        if (debug_print) std::cout << "Fitting quadratic with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        double t0=0,t1=1;
        curve->GetDomain( &t0, &t1);
        DPRINTLN2(t0,t1);
        DPRINTLN(curve->SpanCount());
        std::vector<double> spanKnots(curve->SpanCount()+1);
        curve->GetSpanVector(&spanKnots[0]);
        DPRINTLN(spanKnots);
        DPRINTLN(curve->HasNurbForm());
        DPRINTLN(curve->PointAtStart());
        DPRINTLN(curve->PointAtEnd());
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->PointAt(cf.U()[i]));
          }

        double tol = 1.e-5;
        int kk=0;
        // can't expect this type of curve to intepolate the knots - see below - must use
        // GetClosestPoint
        if (0)
          {
            for (size_t i=1; i < cf.U().size()-1; i += 2)
              {
                Point3D pt = curve->PointAt(cf.U()[i]);
                Point3D Qkk = Q[kk++];
                DPRINTLN2(pt, Qkk);
                double dist = pt.DistanceTo(Qkk);
                EXPECT_NEAR(dist, 0, tol);
              }
          }
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->DerivativeAt(cf.U()[i]));
          }
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->CurvatureAt(cf.U()[i]));
          }

        for (size_t i=0; i < cf.U().size(); i++)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(cf.U()[i], point, first_derivative, second_derivative);
            if (debug_print)
              std::cout << PR(i) << PR(cf.U()[i]) << PR(point)
                        << PR(first_derivative) << PR(second_derivative) << std::endl;
          }

        ((ON_NurbsCurve *)curve)->SetProjectionTolerance(1.e-8);
        for (size_t i = 0; i < Q.size(); i++)
        {
          Point3D p = Q[i];
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN(closest_point);
          double dist = closest_point.DistanceTo(p);
          DPRINTLN(dist);
          EXPECT_NEAR(dist, 0.0, tol);
        }

        {
          Point3D p(0.1,.01,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          EXPECT_NEAR(dist, 0.0, tol);
          EXPECT_NEAR(closest_point[0], 0.1, tol);
          EXPECT_NEAR(closest_point[1], 0.01, tol);
        }

        {
          Point3D p(0.5,0.25,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          EXPECT_NEAR(dist, 0, tol);
        }

        {
          Point3D p(1,1,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          EXPECT_NEAR(dist, 0,tol);
        }

        delete curve;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(unit_stk_geom, test_five_point)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        /// fit a quadratic
        LocalCubicSplineFit cf(BSplineFit::FivePoint);
        int n = 14;
        Vectors2D Q(3);
        Q[0] = Vector2D(-2,0);
        Q[1] = Vector2D(-1,0);
        Q[2] = Vector2D(-.5,0);
        for (int i=0; i < n; i++)
          {
            double x = double(i)/double(n-1);
            double y = x*x;
            Q.push_back( Vector2D(x,y) );
          }
        if (debug_print) std::cout << "Fitting quadratic with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        delete curve;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(unit_stk_geom, test_five_point_1)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        /// fit a quadratic
        LocalCubicSplineFit cf (BSplineFit::FivePoint);
        int n = 7;
        Vectors2D Q(n);
        Q[0] = Vector2D(-2,0);
        Q[1] = Vector2D(-1,0);
        Q[2] = Vector2D(-.5,0);
        Q[3] = Vector2D(0,0);
        Q[4] = Vector2D(1,1);
        Q[5] = Vector2D(2,2);
        Q[6] = Vector2D(2.5,2.5);
        if (debug_print) std::cout << "Fitting quadratic with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        delete curve;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(unit_stk_geom, test_five_point_2)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        /// fit a quadratic
        LocalCubicSplineFit cf (BSplineFit::ThreePoint);
        int n = 7;
        Vectors2D Q(n);
        Q[0] = Vector2D(-3,0);
        Q[1] = Vector2D(-2,0.3);
        Q[2] = Vector2D(-1,0.3);
        Q[3] = Vector2D(0,0);
        Q[4] = Vector2D(1,1);
        Q[5] = Vector2D(2,2);
        Q[6] = Vector2D(2.5,2.5);
        if (debug_print) std::cout << "Fitting quadratic with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        std::vector<int> isCorner(n, 0);
        isCorner[3] = 1;
        cf.setIsCorner(isCorner);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        delete curve;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(unit_stk_geom, test_corner)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        /// fit a quadratic
        LocalCubicSplineFit cf (BSplineFit::ThreePoint);
        int n = 8;
        Vectors2D Q(n);
        Q[0] = Vector2D(-3.2,3);
        Q[1] = Vector2D(-2.2,2);
        Q[2] = Vector2D(-1.2,1);
        Q[3] = Vector2D(-.2,0);
        Q[4] = Vector2D(.2,0);
        Q[5] = Vector2D(1.2,1);
        Q[6] = Vector2D(2.2,2);
        Q[7] = Vector2D(3.2,3);

        int offset=0;
        if (1)
          {
            Q.insert(Q.begin()+4, Vector2D(0,0));
            ++n;
            offset=1;
          }
        if (debug_print) std::cout << "Fitting quadratic with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        std::vector<int> isCorner(n, 0);
        isCorner[3] = 1;
        isCorner[4+offset] = 1;
        cf.setIsCorner(isCorner);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        delete curve;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(unit_stk_geom, test_periodic)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        /// fit a circle
        LocalCubicSplineFit cf;
        cf.setIsPeriodic(true);
        int n = 13;
        Vectors2D Q(n);
        const double r = 1.0;
        for (int i=0; i < n; i++)
          {
            double theta = 2.0*M_PI*double(i)/double(n-1);
            double x = r*std::cos(theta);
            double y = r*std::sin(theta);
            Q[i] = Vector2D(x,y);
          }
        if (debug_print) std::cout << "Fitting circle with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        double t0=0,t1=1;
        curve->GetDomain( &t0, &t1);
        DPRINTLN2(t0,t1);
        DPRINTLN(curve->SpanCount());
        std::vector<double> spanKnots(curve->SpanCount()+1);
        curve->GetSpanVector(&spanKnots[0]);
        DPRINTLN(spanKnots);
        DPRINTLN(curve->HasNurbForm());
        DPRINTLN(curve->PointAtStart());
        DPRINTLN(curve->PointAtEnd());
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->PointAt(cf.U()[i]));
          }

        double tol = 1.e-5;
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->DerivativeAt(cf.U()[i]));
          }
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->CurvatureAt(cf.U()[i]));
          }

        for (size_t i=0; i < cf.U().size(); i++)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(cf.U()[i], point, first_derivative, second_derivative);
            if (debug_print)
              std::cout << PR(i) << PR(cf.U()[i]) << PR(point)
                        << PR(first_derivative) << PR(second_derivative) << std::endl;
          }

        ((ON_NurbsCurve *)curve)->SetProjectionTolerance(1.e-8);
        for (size_t i = 0; i < Q.size(); i++)
        {
          Point3D p = Q[i];
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN(closest_point);
          double dist = closest_point.DistanceTo(p);
          DPRINTLN(dist);
          EXPECT_NEAR(dist, 0.0, tol);
        }

        {
          Point3D p(1.01,0.0,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          EXPECT_NEAR(dist, 0.01, tol);
          EXPECT_NEAR(closest_point[0], 1.0, tol);
          EXPECT_NEAR(closest_point[1], 0.0, tol);
        }

        {
          Point3D p(0,1.1,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          EXPECT_NEAR(dist, 0.1, tol);
        }

        delete curve;
      }

      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(unit_stk_geom, test_closed)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        double naca_0012_21pts[][2] =
          {{1,  0},
           {0.97552826,  0.00341331},
           {0.9045085, 0.01277464},
           {0.79389263,  0.02590486},
           {0.6545085, 0.04009273},
           {0.5, 0.05231025},
           {0.3454915, 0.0591394},
           {0.20610737,  0.05745444},
           {0.0954915, 0.04592861},
           {0.02447174,  0.02586248},
           {0, 0},
           {0.02447174,  -0.02586248},
           {0.0954915, -0.04592861},
           {0.20610737,  -0.05745443},
           {0.3454915, -0.0591394},
           {0.5, -0.05231025},
           {0.6545085, -0.04009273},
           {0.79389263,  -0.02590486},
           {0.9045085, -0.01277464},
           {0.97552826,  -0.00341331},
           {1, 0}};

        /// fit a circle
        LocalCubicSplineFit cf;
        cf.setIsPeriodic(false);
        int n = 21;
        Vectors2D Q(n);
        for (int i=0; i < n; i++)
          {
            Q[i] = Vector2D(naca_0012_21pts[i][0], naca_0012_21pts[i][1]);
          }
        if (debug_print) std::cout << "Fitting naca_0012_21pts with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        double t0=0,t1=1;
        curve->GetDomain( &t0, &t1);
        DPRINTLN2(t0,t1);
        DPRINTLN(curve->SpanCount());
        std::vector<double> spanKnots(curve->SpanCount()+1);
        curve->GetSpanVector(&spanKnots[0]);
        DPRINTLN(spanKnots);
        DPRINTLN(curve->HasNurbForm());
        DPRINTLN(curve->PointAtStart());
        DPRINTLN(curve->PointAtEnd());
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->PointAt(cf.U()[i]));
          }

        double tol = 1.e-4;
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->DerivativeAt(cf.U()[i]));
          }
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->CurvatureAt(cf.U()[i]));
          }

        for (size_t i=0; i < cf.U().size(); i++)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(cf.U()[i], point, first_derivative, second_derivative);
            if (debug_print)
              std::cout << PR(i) << PR(cf.U()[i]) << PR(point)
                        << PR(first_derivative) << PR(second_derivative) << std::endl;
          }

        ((ON_NurbsCurve *)curve)->SetProjectionTolerance(1.e-8);
        for (size_t i = 0; i < Q.size(); i++)
        {
          Point3D p = Q[i];
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN(closest_point);
          double dist = closest_point.DistanceTo(p);
          DPRINTLN(dist);
          EXPECT_NEAR(dist, 0.0, tol);
        }

        {
          Point3D p(1.01,0.0,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          EXPECT_NEAR(dist, 0.01, tol);
          EXPECT_NEAR(closest_point[0], 1.0, tol);
          EXPECT_NEAR(closest_point[1], 0.0, tol);
        }

        {
          Point3D p(-0.01,0.,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          EXPECT_NEAR(dist, 0.01, tol);
        }

        delete curve;
      }

      void test_it(double points[][2], int npts, std::vector<int>& corners, const bool debug_print = false)
      {

        SplineFit::s_debug_print = debug_print;

        LocalCubicSplineFit cf;
        cf.setIsPeriodic(false);
        int n = npts;
        Vectors2D Q(n);
        for (int i=0; i < n; i++)
          {
            Q[i] = Vector2D(points[i][0], points[i][1]);
          }
        if (debug_print) std::cout << "Fitting airfoil with cubic spline - points = " << std::endl;
        DPRINTLN(Q);

        cf.setIsCorner(corners);

        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        double t0=0,t1=1;
        curve->GetDomain( &t0, &t1);
        int ndiv = 10;
        Vectors2D Qref;
        std::vector<double> tv;
        for (size_t i=0; i < cf.U().size()-1; i++)
          {
            tv.push_back(cf.U()[i]);
            for (int id=1; id <= ndiv-1; ++id)
              {
                double t = cf.U()[i] + double(id)/double(ndiv)*(cf.U()[i+1] - cf.U()[i]);
                tv.push_back(t);
              }
          }
        tv.push_back(t1);
        percept::Util::makeUnique(tv);
        if (debug_print) std::cout << "tv= " << tv << std::endl;
        for (unsigned ii=0; ii < tv.size(); ++ii)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(tv[ii], point, first_derivative, second_derivative);
            Vector2D v(point.x, point.y);
            Qref.push_back(v);
          }
        if (debug_print) std::cout << "Fitting airfoil with cubic spline - refined = " << std::endl;
        DPRINTLN(Qref);
#if 1
            std::ofstream ff("/ascldap/users/srkenno/xfer/scraps.txt");
            ff  << "Q= " << Q << std::endl;
            ff  << "Qref2= " << Qref << std::endl;
#endif
        if (debug_print)
          {
            std::cout  << "Qref2= " << Qref << std::endl;
          }

        DPRINTLN2(t0,t1);
        DPRINTLN(curve->SpanCount());
        std::vector<double> spanKnots(curve->SpanCount()+1);
        curve->GetSpanVector(&spanKnots[0]);
        DPRINTLN(spanKnots);
        DPRINTLN(curve->HasNurbForm());
        DPRINTLN(curve->PointAtStart());
        DPRINTLN(curve->PointAtEnd());
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->PointAt(cf.U()[i]));
          }

        double tol = 1.e-4;
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->DerivativeAt(cf.U()[i]));
          }
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->CurvatureAt(cf.U()[i]));
          }

        for (size_t i=0; i < cf.U().size(); i++)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(cf.U()[i], point, first_derivative, second_derivative);
            if (debug_print)
              std::cout << PR(i) << PR(cf.U()[i]) << PR(point)
                        << PR(first_derivative) << PR(second_derivative) << std::endl;
          }

        ((ON_NurbsCurve *)curve)->SetProjectionTolerance(1.e-8);
        for (size_t i = 0; i < Q.size(); i++)
          {
            Point3D p = Q[i];
            double u=0;
            bool success = curve->GetClosestPoint(p, &u);
            EXPECT_TRUE(success);
            Point3D closest_point = curve->PointAt(u);
            DPRINTLN(closest_point);
            double dist = closest_point.DistanceTo(p);
            DPRINTLN(dist);
            EXPECT_NEAR(dist, 0.0, tol);
          }

        if (0)
          {
            Point3D topTE(1, 0.00025, 0);
            Point3D topTEm(0.990190797590092, 0.00219410540613287, 0);
            Point3D p = 0.5*(topTE+topTEm);
            std::cout << "point being projected = " << p << std::endl;
            double u=0;
            bool success = curve->GetClosestPoint(p, &u);
            EXPECT_TRUE(success);
            Point3D closest_point = curve->PointAt(u);
            DPRINTLN2(u,closest_point);
            double dist = closest_point.DistanceTo(p);
            DPRINTLN(dist);
          }

        delete curve;

      }

      TEST(unit_stk_geom, test_closed2)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        double points[][2] =
          {  {-8.6688924, -5.0029264},   {-8.660254, -5},   {-8.6534005, -4.9939822},   {-8.6488179, -4.9860635},   {-8.6460052, -4.9773397},
             {-8.6446617, -4.9703024},   {-8.6444374, -4.9631208},   {-8.6446886, -4.9559395},
             {-8.645216, -4.9487721},   {-8.6459287, -4.9416206},   {-8.6467765, -4.9344837},   {-8.6483853, -4.9259673},   {-8.6501983, -4.9174934},
             {-8.6522198, -4.909066},   {-8.6544529, -4.9006923},   {-8.6569007, -4.8923789},
             {-8.6595656, -4.8841325},   {-8.665133, -4.868257},   {-8.6710449, -4.8525063},   {-8.6773016, -4.8368908},   {-8.683904, -4.8214174},
             {-8.6908518, -4.8060952},   {-8.6981427, -4.7909349},   {-8.7055575, -4.7759944},
             {-8.7131503, -4.7611436},   {-8.7209204, -4.7463848},   {-8.7288668, -4.7317202},   {-8.7369886, -4.717152},   {-8.7452848, -4.7026824},
             {-8.7623129, -4.6739792},   {-8.7797119, -4.6454994},   {-8.7974787, -4.6172475},
             {-8.8156097, -4.5892278},   {-8.8341012, -4.5614448},   {-8.8529491, -4.5339023},   {-8.8720808, -4.5063929},   {-8.8914044, -4.4790178},
             {-8.9109184, -4.4517782},   {-8.9306217, -4.424675},   {-8.9505127, -4.3977094},
             {-8.9705901, -4.3708822},   {-8.9908884, -4.3440402},   {-9.0113346, -4.3173105},   {-9.0319277, -4.2906938},   {-9.0526668, -4.2641907},
             {-9.0735509, -4.2378016},   {-9.094579, -4.2115271},   {-9.1051948, -4.1983867},
             {-9.1158516, -4.1852795},   {-9.1265493, -4.1722057},   {-9.1372877, -4.1591651},   {-9.1480665, -4.1461581},   {-9.1588857, -4.1331846},
             {-9.1616224, -4.1347646},   {-9.1557966, -4.150621},   {-9.1499216, -4.1664593},
             {-9.1439973, -4.1822793},   {-9.1380239, -4.1980807},   {-9.1320011, -4.2138634},   {-9.1259291, -4.2296271},   {-9.1136888, -4.2609752},
             {-9.1012772, -4.2922559},   {-9.0886943, -4.323468},   {-9.0759402, -4.3546106},
             {-9.0630147, -4.3856823},   {-9.049918, -4.4166822},   {-9.0367237, -4.4474833},   {-9.0233163, -4.4781923},   {-9.0096959, -4.5088074},
             {-8.9958627, -4.5393269},   {-8.9818169, -4.5697491},   {-8.9675589, -4.6000723},
             {-8.9531304, -4.6301663},   {-8.9383153, -4.6600719},   {-8.9231151, -4.6897837},   {-8.9075316, -4.7192961},   {-8.8915669, -4.7486041},
             {-8.8752233, -4.7777024},   {-8.8668403, -4.7921219},   {-8.8582848, -4.8064397},
             {-8.8495581, -4.8206538},   {-8.8406616, -4.8347623},   {-8.8315968, -4.8487633},   {-8.8223654, -4.8626549},   {-8.8128816, -4.8765492},
             {-8.8030861, -4.8902272},   {-8.792987, -4.9036818},   {-8.7825919, -4.916908},
             {-8.7719073, -4.9299032},   {-8.7609425, -4.9426625},   {-8.7551333, -4.9490935},   {-8.7491576, -4.9553701},   {-8.7430223, -4.9614909},
             {-8.7367348, -4.9674552},   {-8.7303026, -4.9732623},   {-8.7237315, -4.9789137},
             {-8.7179747, -4.9832164},   {-8.7121377, -4.9874094},   {-8.7061943, -4.9914499},   {-8.7001006, -4.995258},   {-8.6937691, -4.9986546},
             {-8.6870029, -5.0010097},   {-8.6780415, -5.0029358},   {-8.6688924, -5.0029264}
          };

        int n = sizeof(points)/(2*sizeof(double));

        std::vector<int> corners(n,0);
        corners[52] = 1;
        corners[53] = 1;
        test_it(points, n, corners, debug_print);
      }

      TEST(unit_stk_geom, test_closed1)
      {
        const bool debug_print = false;
        SplineFit::s_debug_print = debug_print;

        double naca64418[][2] =
{  {0.99002628828708583164, 0.00047397265058935026403},   {1, -0.0002500000000000000052},   {1, -0.00021000000000000000871},   {1, -0.00012717143142914233596},   {1, 0},   {1, 0.00012717143142914233596},   {1, 0.00021000000000000000871},   {1, 0.0002500000000000000052},
  {0.99019079759009220876, 0.0021941054061328656329},   {0.97939672348457318396, 0.0043456828769943721055},   {0.96754799652658007858, 0.006738140551001466777},   {0.95457740200951124443, 0.0094148844216804614821},   {0.94042220338285342773, 0.012429663166003478841},   {0.92502003395103205019, 0.015817681683521832986},   {0.90830984484028376436, 0.01959191961047875713},   {0.89023626787365739421, 0.023753451853701704893},
  {0.87075460505293444946, 0.02830161685569335811},   {0.84983396894195406546, 0.033234802411729415239},   {0.82745831813032222346, 0.038542107530149104233},   {0.80363057603858389921, 0.044208562324423494871},   {0.77837265041734315929, 0.050203999904281197686},   {0.75172087002855070281, 0.056453948516654445244},   {0.72373225546400277164, 0.062857167391083906827},   {0.69449509754184657329, 0.069325707963627816732},
  {0.66412361369131334143, 0.075761518303150560127},   {0.6327543626068404814, 0.082042775691104977143},   {0.60054920175778159219, 0.088045607110938214901},   {0.56769284071806480618, 0.09364737841273459984},   {0.53438703852640079273, 0.098715779356176425186},   {0.5008456368711438067, 0.10309988954840600128},   {0.46729422943978454752, 0.10665655576324731268},   {0.43397194987540377298, 0.10933287316996238847},
  {0.40110407853824148194, 0.11087708945996455068},   {0.36892533482516109977, 0.11095875686748987721},   {0.33767022574360006093, 0.10991989995920561918},   {0.3075221641753079771, 0.10808763563184944911},   {0.27863797169662024578, 0.10551679761701227012},   {0.25114270211111250353, 0.10230043857842227584},   {0.22512691322579930775, 0.098561404261393162352},   {0.20065073074612110871, 0.094412461356085150554},
  {0.17774665516336701776, 0.089954064851589854435},   {0.15642853707453546752, 0.085247405797907363501},   {0.13668950157889939168, 0.08034743959251201828},   {0.11850438567566352888, 0.075305807121968007523},   {0.10183232476853765203, 0.070171277722418878842},   {0.086617038650746991379, 0.064997167293847066261},   {0.072798753498379831228, 0.059813935349415842113},   {0.060314165470808228653, 0.054641543390762013777},
  {0.049102264687040486635, 0.049485014440722036033},   {0.039096363517174340108, 0.044359183588986307589},   {0.030236433423895866385, 0.039267147284332462598},   {0.022458938379117655143, 0.034222952210548600316},   {0.015707084232797447615, 0.029236942065808409208},   {0.010019635813366760055, 0.024217173961720041969},   {0.005431695309709636682, 0.019136007698792076304},   {0.0020489450882541786376, 0.013984042174987025217},
  {0.00020410562666559633918, 0.0087534448761868235978},   {-0.00035033511722699597251, 0.0037944450508796629187},   {0.00025537761657368971942, -0.0011558285011150260142},   {0.0022877633538304483971, -0.0063092248219418691629},   {0.0057895193014080286592, -0.011365340644720967403},   {0.01072332838221700009, -0.016071223344123291427},   {0.016967754358253032965, -0.020318052163419461631},   {0.024339606193296229458, -0.02423550027702111101},
  {0.03277017241382997742, -0.027947367318039997852},   {0.042254174541994933556, -0.031546730447824189447},   {0.052843416550667987597, -0.035033235712240373694},   {0.064586014642605479863, -0.038434124382554238353},   {0.077540000101318395931, -0.041756746999623367556},   {0.09176294717196496753, -0.045011355006699879655},   {0.10731466659203625635, -0.048186887993079074999},   {0.12425068037641283369, -0.051263190485404348806},
  {0.14261827452153957863, -0.054215638720243236603},   {0.16245145914992878411, -0.057031129744904179857},   {0.18377167556204154764, -0.059687102101209338345},   {0.20658629441019629724, -0.062135616702981809334},   {0.2308823038757764623, -0.064326275993137352338},   {0.25662423684751717312, -0.066201828967897308198},   {0.28375207116483019965, -0.06769416328162937424},   {0.31218001097248571085, -0.068701775916454793647},
  {0.34179248829714869995, -0.069148527608371451736},   {0.37244604063915648373, -0.068998069262887540276},   {0.40396431936271415264, -0.067948207091934914592},   {0.43612842972917137407, -0.065698082345921879344},   {0.46873781798448782565, -0.062593365483713553354},   {0.50159300286180530382, -0.058887398385269650036},   {0.5344760369259542454, -0.05464349647854067904},   {0.56717863440844307199, -0.049996226620562016218},
  {0.5995018446399557055, -0.045093793935311494991},   {0.63125613583719863975, -0.040067491105955144182},   {0.66226103687797865938, -0.035007022126730548417},   {0.69235692013529193112, -0.030012728031143430518},   {0.72140845019166544017, -0.025197638094096464628},   {0.74929776419191285175, -0.02063500784971983304},   {0.77593045418062112084, -0.016388829465100260979},   {0.80124300331585196489, -0.012562565185579830512},
  {0.82518508874050966462, -0.0091930936548549838788},   {0.84771037950456895604, -0.0061791789085143822291},   {0.86880158240480642728, -0.00346332249124713153},   {0.88848382018411453664, -0.0011748376029786177475},   {0.90679174274713503223, 0.00057067652435943455776},   {0.92376005610720912653, 0.0017232011544815475056},   {0.93942514503353724553, 0.0022826905047788406677},   {0.95382757744262558486, 0.0022866684086450535067},
  {0.96701657450649303183, 0.0018611136536592579743},   {0.97905764462256894554, 0.0012085019055270856698},   {0.99002628828708583164, 0.00047397265058935026403}};

        LocalCubicSplineFit cf;
        cf.setIsPeriodic(false);
        int n = 107;
        int n1 = sizeof(naca64418)/(2*sizeof(double));
        VERIFY_OP_ON(n, ==, n1, "bad");
        Vectors2D Q(n);
        for (int i=0; i < n; i++)
          {
            Q[i] = Vector2D(naca64418[i][0], naca64418[i][1]);
          }
        if (debug_print) std::cout << "Fitting airfoil with cubic spline - points = " << std::endl;
        DPRINTLN(Q);

        std::vector<int> corners(n,0);
        corners[1] = 1;
        corners[7] = 1;
        cf.setIsCorner(corners);

        ON_Curve *curve = cf.fit(Q);
        if (debug_print)
          cf.print();

        double t0=0,t1=1;
        curve->GetDomain( &t0, &t1);
        int ndiv = 10;
        Vectors2D Qref;
        std::vector<double> tv;
        for (size_t i=0; i < cf.U().size()-1; i++)
          {
            tv.push_back(cf.U()[i]);
            for (int id=1; id <= ndiv-1; ++id)
              {
                double t = cf.U()[i] + double(id)/double(ndiv)*(cf.U()[i+1] - cf.U()[i]);
                tv.push_back(t);
              }
          }
        tv.push_back(t1);
        percept::Util::makeUnique(tv);
        if (debug_print) std::cout << "tv= " << tv << std::endl;
        for (unsigned ii=0; ii < tv.size(); ++ii)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(tv[ii], point, first_derivative, second_derivative);
            Vector2D v(point.x, point.y);
            Qref.push_back(v);
          }
        if (debug_print) std::cout << "Fitting airfoil with cubic spline - refined = " << std::endl;
        DPRINTLN(Qref);
#if 1
            std::ofstream ff("/ascldap/users/srkenno/xfer/scraps.txt");
            ff  << "Q= " << Q << std::endl;
            ff  << "Qref2= " << Qref << std::endl;
#endif
        if (debug_print)
          {
            std::cout  << "Qref2= " << Qref << std::endl;
          }

        DPRINTLN2(t0,t1);
        DPRINTLN(curve->SpanCount());
        std::vector<double> spanKnots(curve->SpanCount()+1);
        curve->GetSpanVector(&spanKnots[0]);
        DPRINTLN(spanKnots);
        DPRINTLN(curve->HasNurbForm());
        DPRINTLN(curve->PointAtStart());
        DPRINTLN(curve->PointAtEnd());
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->PointAt(cf.U()[i]));
          }

        double tol = 1.e-4;
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->DerivativeAt(cf.U()[i]));
          }
        for (size_t i=0; i < cf.U().size(); i++)
          {
            DPRINTLN2(i, curve->CurvatureAt(cf.U()[i]));
          }

        for (size_t i=0; i < cf.U().size(); i++)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(cf.U()[i], point, first_derivative, second_derivative);
            if (debug_print)
              std::cout << PR(i) << PR(cf.U()[i]) << PR(point)
                        << PR(first_derivative) << PR(second_derivative) << std::endl;
          }

        ((ON_NurbsCurve *)curve)->SetProjectionTolerance(1.e-8);
        for (size_t i = 0; i < Q.size(); i++)
          {
            Point3D p = Q[i];
            double u=0;
            bool success = curve->GetClosestPoint(p, &u);
            EXPECT_TRUE(success);
            Point3D closest_point = curve->PointAt(u);
            DPRINTLN(closest_point);
            double dist = closest_point.DistanceTo(p);
            DPRINTLN(dist);
            EXPECT_NEAR(dist, 0.0, tol);
          }

        if (1)
          {
            Point3D topTE(1, 0.00025, 0);
            Point3D topTEm(0.990190797590092, 0.00219410540613287, 0);
            Point3D p = 0.5*(topTE+topTEm);
            //std::cout << "point being projected = " << p << std::endl;
            double u=0;
            bool success = curve->GetClosestPoint(p, &u);
            EXPECT_TRUE(success);
            Point3D closest_point = curve->PointAt(u);
            DPRINTLN2(u,closest_point);
            double dist = closest_point.DistanceTo(p);
            DPRINTLN(dist);
            //EXPECT_NEAR(dist, 0.01, tol);
            //EXPECT_NEAR(closest_point[0], 1.01, tol);
            //EXPECT_NEAR(closest_point[1], -2.50000e-04, tol);
          }


        if (0)
          {
            Point3D p(1.01000e+00, -2.50000e-04, 0);
            double u=0;
            bool success = curve->GetClosestPoint(p, &u);
            EXPECT_TRUE(success);
            Point3D closest_point = curve->PointAt(u);
            DPRINTLN2(u,closest_point);
            double dist = closest_point.DistanceTo(p);

            EXPECT_NEAR(dist, 0.01, tol);
            EXPECT_NEAR(closest_point[0], 1.01, tol);
            EXPECT_NEAR(closest_point[1], -2.50000e-04, tol);
          }

        if (0)
          {
            Point3D p(-0.01,0.,0);
            double u=0;
            bool success = curve->GetClosestPoint(p, &u);
            EXPECT_TRUE(success);
            Point3D closest_point = curve->PointAt(u);
            DPRINTLN2(u,closest_point);
            double dist = closest_point.DistanceTo(p);
            EXPECT_NEAR(dist, 0.01, tol);
          }

        delete curve;

      }


      //=============================================================================
      //=============================================================================
      //=============================================================================

      TEST(unit_stk_geom, test_2)
      {
        const bool debug_print = false;
        if (1)
          {
            // write a wiggly cubic curve on the "green NURBS wiggle" layer
            ON_NurbsCurve* wiggle = new ON_NurbsCurve(
                                                      3, // dimension
                                                      false, // true if rational
                                                      4,     // order = degree+1
                                                      6      // number of control vertices
                                                      );
            int i;
            Points3D pts;
            for ( i = 0; i < wiggle->CVCount(); i++ ) {
              Point3D pt( 2*i, -i, (i-3)*(i-3) ); // pt = some 3d point
              DPRINTLN2(i,pt);
              wiggle->SetCV( i, pt );
              pts.push_back(pt);
            }

            // ON_NurbsCurve's have order+cv_count-2 knots.
            wiggle->SetKnot(0, 0.0);
            wiggle->SetKnot(1, 0.0);
            wiggle->SetKnot(2, 0.0);
            wiggle->SetKnot(3, 1.5);
            wiggle->SetKnot(4, 2.3);
            wiggle->SetKnot(5, 4.0);
            wiggle->SetKnot(6, 4.0);
            wiggle->SetKnot(7, 4.0);

            DPRINTLN(wiggle->IsValid());
            DPRINTLN(wiggle->IsRational());
            DPRINTLN(wiggle->CVCount());
            DPRINTLN(wiggle->CVSize());
            DPRINTLN(wiggle->Order());
            DPRINTLN(wiggle->KnotCount());

            DPRINTLN(wiggle->HasNurbForm());
            DPRINTLN(wiggle->PointAtStart());
            DPRINTLN(wiggle->PointAtEnd());

            int kn=wiggle->KnotCount();
            std::vector<double> knots;
            for (int k = 0; k < kn; k++)
              {
                DPRINTLN2(k, wiggle->PointAt(wiggle->Knot(k)));
                knots.push_back(wiggle->Knot(k));
              }
            if (debug_print)
              {
                ON_TextLog log;
                wiggle->Dump(log);
              }

            DPRINTLN(pts);
            DPRINTLN(knots);

            delete wiggle;
          }

        if (1)
          {
            double pts[][2] = {{0, 1},     {0.2, 1.3}, {1.5, 2},
                               {1.75, 2},  {2, 2}, {2.1, 0.5},
                               {2.5, .5},  {2.9, .5}, {3.3, .9},
                               {4, 2}};
            //  Piegl & Tiller convention - start/end with order number of repeated knots
            //double[] knots =  {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7};
            // OpenNURBS convention - start/end with (order-1) repeated knots
            //double knots[] =  {0, 0,   0, 1, 2, 3, 4, 5, 6, 7,   7, 7 };
            double knots[] =  {0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4 };

            // n = 3, Q[0]..Q[n] --> n = 3;
            int np=10;
            //int nk=12;

            ON_NurbsCurve* wiggle = new ON_NurbsCurve(
                                                      2, // dimension
                                                      false, // true if rational
                                                      4,     // order = degree+1
                                                      np      // number of control vertices
                                                      );
            int i;
            for ( i = 0; i < wiggle->CVCount(); i++ ) {
              Point3D pt( pts[i][0], pts[i][1], 0);
              DPRINTLN2(i,pt);
              wiggle->SetCV( i, pt );
            }

            // ON_NurbsCurve's have order+cv_count-2 knots.
            int kn=wiggle->KnotCount();
            for (int k = 0; k < kn; k++)
              {
                wiggle->SetKnot(k, knots[k]);
              }

            DPRINTLN(wiggle->IsValid());
            DPRINTLN(wiggle->IsRational());
            DPRINTLN(wiggle->CVCount());
            DPRINTLN(wiggle->CVSize());
            DPRINTLN(wiggle->Order());
            DPRINTLN(wiggle->KnotCount());

            DPRINTLN(wiggle->HasNurbForm());
            DPRINTLN(wiggle->PointAtStart());
            DPRINTLN(wiggle->PointAtEnd());

            for (int k = 0; k < kn; k++)
              {
                DPRINTLN2(k, wiggle->PointAt(wiggle->Knot(k)));
              }
            delete wiggle;
          }

      }
#endif

    }
  }
