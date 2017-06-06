// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_Percept_hpp
#define percept_Percept_hpp

//#define HAVE_INTREPID_DEBUG 0

#define DO_TIMING 1
#define DO_MEMORY 1
#define ENABLE_SMOOTHER3 0

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- PGI Compiler bug workaround (switch statements in template code)
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

#define USE_PGI_7_1_COMPILER_BUG_WORKAROUND

#if !defined(PGI_INSTANTIATION_FILE) && defined(__PGI) && defined(USE_PGI_7_1_COMPILER_BUG_WORKAROUND)
// ensure bases get included only once
#define INTREPID_HGRAD_HEX_C1_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_C2_FEMDEF_HPP

#define INTREPID_HGRAD_LINE_C1_FEMDEF_HPP

#define INTREPID_HGRAD_TRI_C1_FEMDEF_HPP
#define INTREPID_HGRAD_TRI_C2_FEMDEF_HPP

#define INTREPID_HGRAD_QUAD_C1_FEMDEF_HPP
#define INTREPID_HGRAD_QUAD_C2_FEMDEF_HPP

#define INTREPID_HGRAD_HEX_C1_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_C2_FEMDEF_HPP

#define INTREPID_HGRAD_TET_C1_FEMDEF_HPP
#define INTREPID_HGRAD_TET_C2_FEMDEF_HPP

#define INTREPID_HGRAD_WEDGE_C1_FEMDEF_HPP

#define INTREPID_HGRAD_QUAD_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_HEX_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_WEDGE_C2_SERENDIPITY_FEMDEF_HPP
#define INTREPID_HGRAD_QUAD_C2_SERENDIPITY_FEMDEF_HPP

      // Shells
#define INTREPID_HGRAD_TRI_C1_FEMDEF_HPP
#define INTREPID_HGRAD_TRI_C2_FEMDEF_HPP

#define INTREPID_HGRAD_QUAD_C1_FEMDEF_HPP


#define INTREPID_CELLTOOLSDEF_HPP

#endif


//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- PGI Compiler old-ness problem - unsupported on boost::ublas
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

/* The following is to avoid the following error from boost::ublas libraries:

"/scratch/srkenno/code/TPLs_src/boost/boost/numeric/ublas/detail/config.hpp", line 170: catastrophic error:
          #error directive: Your compiler and/or configuration is unsupported
          by this verions of uBLAS. Define BOOST_UBLAS_UNSUPPORTED_COMPILER=0
          to override this message. Boost 1.32.0 includes uBLAS with support
          for many older compilers.
  #error Your compiler and/or configuration is unsupported by this verions of uBLAS. Define BOOST_UBLAS_UNSUPPORTED_COMPILER=0 to
  override this message. Boost 1.32.0 includes uBLAS with support for many older compilers.
*/

#if defined(__PGI)
#define BOOST_UBLAS_UNSUPPORTED_COMPILER 0
#endif


//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- Geometry configuration
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

#if defined(STK_BUILT_IN_SIERRA) && !defined(STK_PERCEPT_HAS_GEOMETRY)
#define STK_PERCEPT_HAS_GEOMETRY
#define STK_PERCEPT_USE_INTREPID
#endif

#if !defined(STK_BUILT_IN_SIERRA) && defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 0
#  if !defined(STK_PERCEPT_HAS_GEOMETRY)
#    define STK_PERCEPT_HAS_GEOMETRY
#  endif
#  if !defined(STK_PERCEPT_USE_INTREPID)
#    define STK_PERCEPT_USE_INTREPID
#  endif
#endif

#if defined(NO_GEOM_SUPPORT) && defined(STK_PERCEPT_HAS_GEOMETRY)
#undef STK_PERCEPT_HAS_GEOMETRY
#endif

#if defined(STK_PERCEPT_LITE) && STK_PERCEPT_LITE == 1
#if defined(STK_PERCEPT_HAS_GEOMETRY)
#undef STK_PERCEPT_HAS_GEOMETRY
#endif
#if defined(STK_PERCEPT_USE_INTREPID)
#undef STK_PERCEPT_USE_INTREPID
#endif
#endif


//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- YAML configuration
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------


#define STK_ADAPT_ALLOW_YAML_CPP 1
#ifndef STK_ADAPT_HAVE_YAML_CPP
#define STK_ADAPT_HAVE_YAML_CPP (STK_ADAPT_ALLOW_YAML_CPP && defined(STK_BUILT_IN_SIERRA))
#endif

//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//----- Kokkos configuration
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

//#undef WITH_KOKKOS
#define WITH_KOKKOS

#include <percept/ExceptionWatch.hpp>


#endif
