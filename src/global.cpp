/* Copyright (C) 2011,2012 Chi-kwan Chan
   Copyright (C) 2011,2012 NORDITA

   This file is part of fg2.

   Fg2 is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Fg2 is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   License for more details.

   You should have received a copy of the GNU General Public License
   along with fg2.  If not, see <http://www.gnu.org/licenses/>. */

#include <cstdlib>
#include "fg2.h"

namespace global {
  // Values needed in setup()
  Z n1, n2, s;
  Z g1, g2, b1, b2, sz;
  R *u, *v, *host = NULL;
  // Values needed in config()
  R l1 = 1.0, l2 = 1.0, c = 0.5;
  int p1 = 1, p2 = 1;
  double flops = 0.0, bps = 0.0;
}
