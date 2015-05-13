#ifndef _2D_MORPH_H
#define _2D_MORPH_H
/**
 *  Thin Plate Spline 2D point morpher example.
 *
 *    Takes in sample point coordinates and their displacements,
 *    fits a TPS into them and allowes morphing new points
 *    accordint to the approximated deformation.
 *
 *    Supports TPS approximation 3 as suggested in paper
 *    Gianluca Donato and Serge Belongie, 2002: "Approximation
 *    Methods for Thin Plate Spline Mappings and Principal Warps"
 *
 *  This should be considered more as an example than a ready made module!
 *  The code has been extracted from a working "deformed shape matching"
 *  application and has been optimized for that particular case.
 *  I don't even know if this compiles or not.
 *
 *  Copyright (C) 2003-2005 by Jarno Elonen
 *
 *  This is Free Software / Open Source with a very permissive
 *  license:
 *
 *  Permission to use, copy, modify, distribute and sell this software
 *  and its documentation for any purpose is hereby granted without fee,
 *  provided that the above copyright notice appear in all copies and
 *  that both that copyright notice and this permission notice appear
 *  in supporting documentation.  The authors make no representations
 *  about the suitability of this software for any purpose.
 *  It is provided "as is" without express or implied warranty.
 */

#include "ludecomposition.h"

#include <vector>
#include <utility>

#include <cassert>
#include <cmath>

#include <stdio.h>

using namespace boost::numeric;
typedef ublas::matrix<double> Matrix;
typedef ublas::matrix_row<Matrix> Matrix_Row;
typedef ublas::matrix_column<Matrix> Matrix_Col;

inline double SQR( double x ) { return x*x; }

/// A 2D point
struct Point
{
  double x, y;
};

/// A displacement of one 2D point: x,y is the original location
/// and dx,dy is the offset.
struct Coord_Diff
{
  Coord_Diff()
  {}

  Coord_Diff( double x, double y, double dx, double dy )
    : x(x), y(y), dx(dx), dy(dy)
  {}

  double x, y, dx, dy;
};

/// 2D point morph interpolator using TPS (Thin Plate Spline).
class TPS_Morpher
{
public:

  /// Calculate the morph weights from sample points.
  /// Builds a matrix of (approximately) size N(p_samples) x N(p_samples*subsampling_factor),
  /// and inverts it using LU decomposition (O(N^3), so be careful not to input too large sample sets.
  ///
  /// For performance reasons, this function assumes ownership of the sample vector so don't
  /// change it after passing it here. Once you have done using the morpher, you can claim
  /// it back by calling grab_samples().
  ///
  /// The function assumes that the average distance between the points to be about 0.5 to
  /// save some computation. If this is not the case in your application, you need to calculate
  /// (or approximate) this value by yourself. See the comments concerning "double a" in the code.
  ///
  /// @param p_samples Displacement samples to be interpolated
  /// @param regularization Amount of "relaxation", 0.0 = exact interpolation
  /// @param subsampling_factor 1.0 = use all points, 0.5 = use 1/2 of the points etc.
  TPS_Morpher(
    std::auto_ptr<std::vector<Coord_Diff> > p_samples,
    double regularization,
    double subsampling_factor = 1.0 )
      : samples( p_samples )
  {
    assert( samples->size() >= 3 );
    unsigned p = samples->size();

    unsigned m = (unsigned)(p * subsampling_factor);
    if ( m < 3 ) m = 3;
    if ( m > p ) m = p;

    if ( m < p )
    {
      // Randomize the input if subsampling is used
      for ( unsigned i=0; i<samples->size()-1; ++i )
      {
        int j = i + ((unsigned)rand()) % (samples->size()-i);
        Coord_Diff tmp = (*samples)[j];
        (*samples)[i] = (*samples)[j];
        (*samples)[j] = tmp;
      }
    }

    // Allocate the matrix and vector
    mtx_l.resize(p+3, m+3);
    mtx_v.resize(p+3, 2);
    mtx_orig_k.resize(p, m);

    // Fill K (p x m, upper left of L)
    for ( unsigned i=0; i<p; ++i )
    {
      const Coord_Diff& pt_i = (*samples)[i];
      for ( unsigned j=0; j<m; ++j )
      {
        const Coord_Diff& pt_j = (*samples)[j];
        double elen2 = SQR(pt_i.x-pt_j.x) + SQR(pt_i.y-pt_j.y);
        mtx_l(i,j) = mtx_orig_k(i,j) = base_func(elen2);
      }
    }

    // Empiric value for avg. distance between points
    //
    // This variable is normally calculated to make regularization
    // scale independent, but since our shapes in this application are always
    // normalized to maxspect [-.5,.5]x[-.5,.5], this approximation is pretty
    // safe and saves us p*p square roots
    const double a = 0.5;

    // Fill the rest of L
    for ( unsigned i=0; i<p; ++i )
    {
      const Coord_Diff pt_i = (*samples)[i];

      // P (p x 3, upper right)
      mtx_l(i, m+0) = 1.0;
      mtx_l(i, m+1) = pt_i.x;
      mtx_l(i, m+2) = pt_i.y;

      if ( i<m )
      {
        // diagonal: reqularization parameters (lambda * a^2)
        mtx_l(i,i) = mtx_orig_k(i,i) =
          regularization * (a*a);

        // P transposed (3 x p, bottom left)
        mtx_l(p+0, i) = 1.0;
        mtx_l(p+1, i) = pt_i.x;
        mtx_l(p+2, i) = pt_i.y;
      }
    }

    // O (3 x 3, lower right)
    for ( unsigned i=p; i<p+3; ++i )
      for ( unsigned j=m; j<m+3; ++j )
        mtx_l(i,j) = 0.0;

    // Fill the right hand matrix V
    for ( unsigned i=0; i<p; ++i )
    {
      const Coord_Diff& pt_i = (*samples)[i];
      mtx_v(i,0) = pt_i.dx;
      mtx_v(i,1) = pt_i.dy;
    }

    mtx_v(p+0, 0) = mtx_v(p+1, 0) = mtx_v(p+2, 0) = 0.0;
    mtx_v(p+0, 1) = mtx_v(p+1, 1) = mtx_v(p+2, 1) = 0.0;

    // Solve the linear system "inplace"
    int sret = LU_Solve(mtx_l, mtx_v);
    assert( sret != 2 );
    if (sret == 1)
    {
      puts( "Singular matrix! Aborting." );
      exit(1);
    }
  }

  /// Morph given points according to the TPS
  /// @param pts The points to morph
  void morph( std::vector<Point>& pts )
  {
    assert( samples.get() && "Morpher no longer owns 'samples'");

    const unsigned m = mtx_orig_k.size2();
    for ( std::vector<Point>::iterator ite=pts.begin(), end=pts.end();
          ite != end;
          ++ite )
    {
      double x=ite->x, y=ite->y;
      double dx = mtx_v(m+0, 0) + mtx_v(m+1, 0)*x + mtx_v(m+2, 0)*y;
      double dy = mtx_v(m+0, 1) + mtx_v(m+1, 1)*x + mtx_v(m+2, 1)*y;

      std::vector<Coord_Diff>::const_iterator diff_ite = samples->begin();
      Matrix_Col cv0(mtx_v,0), cv1(mtx_v,1);
      Matrix_Col::const_iterator cv0_ite(cv0.begin()), cv1_ite(cv1.begin());
      for ( unsigned i=0; i<m; ++i, ++diff_ite, ++cv0_ite, ++cv1_ite )
      {
        double d = base_func( SQR(diff_ite->x - x) + SQR(diff_ite->y - y) );
        dx += (*cv0_ite) * d;
        dy += (*cv1_ite) * d;
      }

      ite->x += dx;
      ite->y += dy;
    }
  }

  /// Calculate bending energy, or if subsampling_factor
  /// for constructor was < 1.0, approximate it.
  double calc_bending_energy()
  {
    assert( samples.get() && "Morpher no longer owns 'samples'");

    // bending energy = trace( W^T * A * W ),
    // where A = upper left m x m block of mtx_orig_k
    const unsigned m = mtx_orig_k.size2();
    ublas::matrix_range<Matrix> mtx_w(mtx_v,
      ublas::range(0, m),
      ublas::range(0, 2));
    ublas::matrix_range<Matrix> mtx_a(mtx_orig_k,
      ublas::range(0, m),
      ublas::range(0, m));
    Matrix bm = prod( Matrix(prod(trans(mtx_w), mtx_a)), mtx_w);
    assert( bm.size1() == bm.size2() && bm.size1() == 2 );
    return bm(0,0) + bm(1,1);
  }

  /// Takes away the ownership of 'samples' from the morpher.
  /// After this, calling other functions becomes illegal.
  std::auto_ptr<std::vector<Coord_Diff> > grab_samples()
  {
    return samples;
  }

private:

  inline double base_func(double r2)
  {
    // same as r*r * log(r), but for r^2:
    return ( r2==0 )
      ? 0.0 // function limit at 0
      : r2 * log(r2) * 0.217147241; // = 1/(2*log(10))
  }

  std::auto_ptr<std::vector<Coord_Diff> > samples;
  Matrix mtx_l, mtx_v, mtx_orig_k;
};

#endif
