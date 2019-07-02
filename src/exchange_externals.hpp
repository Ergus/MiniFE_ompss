#ifndef EXCHANGE_EXTERNALS_HPP
#define EXCHANGE_EXTERNALS_HPP

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

#include <cstdlib>
#include <iostream>

#include <outstream.hpp>
#include "Vector.hpp"

namespace miniFE {

	template<typename MatrixType>
	void exchange_externals_all(MatrixType *A_array, Vector *x_array, size_t numboxes)
	{

		#ifdef MINIFE_DEBUG
		std::ostream& os = outstream();
		os << "entering exchange_externals\n";
		#endif

		if (numboxes < 2)
			return;

		for (size_t id = 0; id < numboxes; ++id) {

			MatrixType &A = A_array[id];
			Vector &x = x_array[id];

			// TODO: Task here (will create a nested)
			// in A
			// in x
			{
				// in A.elements_to_send[0;A.nelements_to_send];
				// in x.coefs[]
				// out A.send_buffer[0;A.nelements_to_send];
				{
					for (int i = 0; i < A.nelements_to_send; ++i)
						A.send_buffer[i] = x.coefs[A.elements_to_send[i]];
				}
			}
		}

		// TODO:
		for (size_t id = 0; id < numboxes; ++id) {
			MatrixType &A = A_array[id];
			Vector &x = x_array[id];

			// TODO: Task here to copy locally, but then there is
			// the problem inside.  I don't know the whole
			// dependencies I'll need in advance.  in A (full) out
			// x.coefs[A.nrows; A.nexternals]
			{
				double *x_external = &(x.coefs[A.nrows]);

				int index = 0;
				// This can be paralelized if needed easily
				for (int i = 0; i < A.nrecv_neighbors; ++i) {
					for (int j = 0; j < A.recv_length[i]; ++j) {
						x_external[index++] = A.recv_ptr[i][j];
					}
				}

				// Assert that we copied all the elements
				assert(index == A.nexternals);
			}
		}
	}


}//namespace miniFE

#endif

