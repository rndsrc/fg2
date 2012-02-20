/* Copyright (C) 2011,2012 Chi-kwan Chan
   Copyright (C) 2011,2012 NORDITA

   This file is part of fd2.

   Fd2 is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Fd2 is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   License for more details.

   You should have received a copy of the GNU General Public License
   along with fd2.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef DERIV_H
#define DERIV_H

__device__ __constant__ R Delta1 = 0.0;
__device__ __constant__ R Delta2 = 0.0;

#ifndef THIS
#define THIS u
#endif

#ifndef STRIDE
#define STRIDE s
#endif

#define _GET_(x, i, j) ( (THIS)[(i) * (STRIDE) + (j)].x  )
#define _SUM_(x, i, j) ( _GET_(x, i, j) + _GET_(x,-i,-j) )
#define _DIF_(x, i, j) ( _GET_(x, i, j) - _GET_(x,-i,-j) )
#define _MIX_(x, i, j) ( _GET_(x, i, j) - _GET_(x,-i, j) \
                       + _GET_(x,-i,-j) - _GET_(x, i,-j) )

#if ORDER == 6 ///////////////////////////////////////////////////////////////

// 1st derivative in one direction: 9 FLOP
#define  D1(x) (( (R)(  45./ 60.) * _DIF_(x, 1, 0) \
                + (R)(  -9./ 60.) * _DIF_(x, 2, 0) \
                + (R)(   1./ 60.) * _DIF_(x, 3, 0) ) / Delta1)
#define  D2(x) (( (R)(  45./ 60.) * _DIF_(x, 0, 1) \
                + (R)(  -9./ 60.) * _DIF_(x, 0, 2) \
                + (R)(   1./ 60.) * _DIF_(x, 0, 3) ) / Delta2)

// 2nd derivative in one direction: 12 FLOP
#define D11(x) (( (R)(-490./180.) * _GET_(x, 0, 0) \
                + (R)( 270./180.) * _SUM_(x, 1, 0) \
                + (R)( -27./180.) * _SUM_(x, 2, 0) \
                + (R)(   2./180.) * _SUM_(x, 3, 0) ) / (Delta1 * Delta1))
#define D22(x) (( (R)(-490./180.) * _GET_(x, 0, 0) \
                + (R)( 270./180.) * _SUM_(x, 0, 1) \
                + (R)( -27./180.) * _SUM_(x, 0, 2) \
                + (R)(   2./180.) * _SUM_(x, 0, 3) ) / (Delta2 * Delta2))

// 2nd derivative in mixed directions: 16 FLOP
#define D12(x) (( (R)( 270./720.) * _MIX_(x, 1, 1) \
                + (R)( -27./720.) * _MIX_(x, 2, 2) \
                + (R)(   2./720.) * _MIX_(x, 3, 3) ) / (Delta1 * Delta2))
#define D21(x) D12(x)

#else ////////////////////////////////////////////////////////////////////////

#error !!! Requested derivatives are not implemented !!!

#endif ///////////////////////////////////////////////////////////////////////

#endif // DERIV_H
