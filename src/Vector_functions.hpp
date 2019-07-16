#ifndef _Vector_functions_hpp_
#define _Vector_functions_hpp_

//@HEADER
// ************************************************************************
//
// MiniFE: Simple Finite Element Assembly and Solve
// Copyright (2006-2013) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// ************************************************************************
//@HEADER

#include <sstream>
#include <fstream>

#include "Vector.hpp"
#include "mytimer.hpp"

namespace miniFE {

	//------------------------------------------------------------
	//Compute the update of a vector with the sum of two scaled vectors where:
	//
	// w = alpha*x + beta*y
	//
	// x,y - input vectors
	//
	// alpha,beta - scalars applied to x and y respectively
	//
	// w - output vector
	//

	inline void waxpby_task(double alpha, const Vector *x,
		double beta,  const Vector *y,
		Vector *w)
	{

		#pragma oss task					\
			in(x[0])						\
			in(x->coefs[0; x->local_size])			\
			in(y[0])						\
			in(y->coefs[0; y->local_size])			\
			in(w[0])						\
			out(w->coefs[0; w->local_size])
		{
			assert(x->local_size <= y->local_size);
			assert(x->local_size <= w->local_size);
			const int n = x->local_size;

			for (int i = 0; i < n; ++i)
				w->coefs[i] = alpha * x->coefs[i] + beta * y->coefs[i];
		}

	}

	//Like waxpby above, except operates on two sets of arguments.
	//In other words, performs two waxpby operations in one loop.

	inline void fused_waxpby(double alpha, const Vector &x,
	                         double beta, const Vector &y,
	                         Vector& w,
	                         double alpha2, const Vector &x2,
	                         double beta2, const Vector &y2,
	                         Vector &w2)
	{

		assert(x.local_size <= y.local_size);
		assert(x.local_size <= w.local_size);
		const int n = x.local_size;

		for(int i = 0; i < n; ++i) {
			w.coefs[i] = alpha * x.coefs[i] + beta * y.coefs[i];
			w2.coefs[i] = alpha2 * x2.coefs[i] + beta2 * y2.coefs[i];
		}
	}

}//namespace miniFE

#endif

