// -*- mode: c++; c-basic-offset: 4 -*-
// (C) Nicholas Polson, James Scott, Jesse Windle, 2012-2019

// This file is part of BayesLogit.

// BayesLogit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// BayesLogit.  If not, see <https://www.gnu.org/licenses/>.



// MAKE SURE YOU USE GetRNGState() and PutRNGState() !!!

#include "polyagamma_wrapper.h"
#include "PolyaGamma.h"
#include "PolyaGammaApproxAlt.h"
#include "PolyaGammaApproxSP.h"
#include "simple_RNG_wrapper.h"


double rpg_hybrid(double b_, double z_)
{
    PolyaGamma          dv(1000);
    //PolyaGammaApproxAlt al;
    PolyaGammaApproxSP  sp;

    double x;

    double b = (double) b_;
    double z = (double) z_;

    if (b > 170) {
	double m = dv.pg_m1(b,z);
	double v = dv.pg_m2(b,z) - m*m;
	x = (double) norm(m, sqrt(v));
    }
    else if (b > 13) {
	sp.draw(x, b, z);
    }
    else if (b==1 || b==2) {
	x = dv.draw((int)b, z);
    }
    else if (b > 0) {
	x = dv.draw_sum_of_gammas(b, z);
    }
    else {
	x = 0.0;
    }

    return (double) x;
}
