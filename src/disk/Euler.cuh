/* Copyright (C) 2012 Chi-kwan Chan
   Copyright (C) 2012 NORDITA

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

static __device__ void first_order(S *u, const R dt, const Z i, const Z j)
{
  const R den    = exp(u->lnd );
  const R eng    = exp(u->lne );
  const R rad    = exp(u->lner);
  const R gamma1 = para_gamma - K(1.0);

  R eq_eng, dt_eng;
  {
    // To reach radiation equilibrium we need to solve the system of
    // equations
    //
    //   er + e == etot
    //   er == ar T^4 / den == ar (gamma - 1)^4 e^4 / den
    //
    // Combining the two equations, we have
    //
    //   [ar (gamma - 1)^4 / den] e^4 + e - etot == 0
    //
    // because ar is a large number, the above equation is not
    // necessary well behave.  Instead, we multiple the whole
    // equation by [ar (gamma - 1)^4 / den]^(1/3) and solve
    //
    //    x^4 + x - x0 == 0
    //
    // where x0 = etot [ar (gamma - 1)^4 / den]^(1/3).

    const R gamma1_2 = gamma1 * gamma1;
    const R scale    = pow(para_ar * gamma1_2 * gamma1_2 / den,
                           K(1.0) / K(3.0));

    const R x0 = (eng + rad) * scale;
    const R a  = pow(sqrt(K(0.00390625) + x0 * x0 * x0 / K(27.0)) + K(0.0625),
                     K(1.0) / K(3.0));
    const R b  = a - x0 / (K(3.0) * a);
    const R c  = sqrt(K(2.0) * b);

    eq_eng = (sqrt(K(0.5) / c - K(0.5) * b) - K(0.5) * c) / scale;
    dt_eng = (eq_eng - eng) / dt;
  }

  const R temp    = gamma1 * eng;
  const R temp2   = temp * temp;
  const R kappa   = para_ff * den * den * pow(temp, K(-3.5));
  const R cooling = kappa * (para_ar * temp2 * temp2 / den - rad);

  if(cooling * cooling > dt_eng * dt_eng) {
    u->lne  = log(            eq_eng);
    u->lner = log(eng + rad - eq_eng);
  }
}
