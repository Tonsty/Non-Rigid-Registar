#ifndef __GAUSSELIM_H
#define __GAUSSELIM_H

/*
 *  Simple Gauss linear system solver for Boost uBlas matrices
 *
 *  Copyright (C) 2003 by Jarno Elonen
 *
 *  Permission to use, copy, modify, distribute and sell this software
 *  and its documentation for any purpose is hereby granted without fee,
 *  provided that the above copyright notice appear in all copies and
 *  that both that copyright notice and this permission notice appear
 *  in supporting documentation.  The authors make no representations
 *  about the suitability of this software for any purpose.
 *  It is provided "as is" without express or implied warranty.
 */

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

 /*
  *  Solves a linear system A x X = B using Gauss elimination,
  *  given by Boost uBlas matrices 'a' and 'b'.
  *
  *  If the elimination succeeds, the function returns true,
  *  A is inverted and B contains values for X.
  *
  *  In case A is singular (linear system is unsolvable), the
  *  function returns false and leaves A and B in scrambled state.
  *
  *  TODO: make further use of uBlas views instead of direct
  *  element access (for better optimizations)
  */
template <class T> bool gauss_solve(
  boost::numeric::ublas::matrix<T>& a,
  boost::numeric::ublas::matrix<T>& b )
{
  int icol, irow;
  int n = a.size1();
  int m = b.size2();

  int* indxc = new int[n];
  int* indxr = new int[n];
  int* ipiv = new int[n];

  typedef boost::numeric::ublas::matrix<T> GJ_Mtx;
  typedef boost::numeric::ublas::matrix_row<GJ_Mtx> GJ_Mtx_Row;
  typedef boost::numeric::ublas::matrix_column<GJ_Mtx> GJ_Mtx_Col;

  for (int j=0; j<n; ++j)
    ipiv[j]=0;

  for (int i=0; i<n; ++i)
  {
    T big=0.0;
    for (int j=0; j<n; j++)
      if (ipiv[j] != 1)
        for (int k=0; k<n; k++)
        {
          if (ipiv[k] == 0)
          {
            T cmpa = a(j,k);
            if ( cmpa < 0)
              cmpa = -cmpa;
            if (cmpa >= big)
            {
              big = cmpa;
              irow=j;
              icol=k;
            }
          }
          else if (ipiv[k] > 1)
            return false;
        }

    ++(ipiv[icol]);
    if (irow != icol)
    {
      GJ_Mtx_Row ar1(a, irow), ar2 (a, icol);
      ar1.swap(ar2);

      GJ_Mtx_Row br1(b, irow), br2(b, icol);
      br1.swap(br2);
    }

    indxr[i] = irow;
    indxc[i] = icol;
    if (a(icol, icol) == 0.0)
      return false;

    T pivinv = 1.0 / a(icol, icol);
    a(icol, icol) = 1.0;
    GJ_Mtx_Row(a, icol) *= pivinv;
    GJ_Mtx_Row(b, icol) *= pivinv;

    for (int ll=0; ll<n; ll++)
      if (ll != icol)
      {
        T dum = a(ll, icol);
        a(ll, icol) = 0.0;
        for (int l=0; l<n; l++)
          a(ll, l) -= a(icol, l) * dum;
        for (int l=0; l<m; l++)
          b(ll, l) -= b(icol, l) * dum;
      }
  }

  // Unscramble A's columns
  for (int l=n-1; l>=0; --l)
    if (indxr[l] != indxc[l])
    {
      GJ_Mtx_Col ac1(a, indxr[l]), ac2(a, indxc[l]);
      ac1.swap(ac2);
    }

  delete[] ipiv;
  delete[] indxr;
  delete[] indxc;

  return true;
}

#endif
