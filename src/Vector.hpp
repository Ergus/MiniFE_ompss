#ifndef _Vector_hpp_
#define _Vector_hpp_

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

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include "ompss_utils.hpp"

namespace miniFE
{
	struct Vector {

		Vector() : startIndex(-1), local_size(-1), coefs(nullptr)
		{
		}

		~Vector()
		{
			rrd_free(coefs, local_size * sizeof(double));
		}

		// Arrays allocations (boxes can go in local memory)
		static void* operator new[](std::size_t sz)
		{
			void * const tmp = rrl_malloc(sz);
			dbprintf("Calling: %s, size: %lu\n", __PRETTY_FUNCTION__, sz);
			return tmp;
		}

		static void operator delete[](void* ptr, std::size_t sz)
		{
			printf("Calling: %s, address %p size: %lu\n", __PRETTY_FUNCTION__, ptr, sz);
			return rrl_free(ptr, sz);
		}

		void init(const int startIdx, const int local_sz)
		{
			startIndex = startIdx;
			local_size = local_sz;
			coefs = (double *) rrd_malloc(local_sz * sizeof(double));

			// TODO: task HERE
			for (int i = 0; i < local_sz; ++i)
				coefs[i] = 0.;
		}

		int startIndex, local_size;
		double *coefs;

		void sum_into_vector(size_t num_indices,
		                     const int *indices, const double *coefs)
		{
			int first = startIndex;
			int last = first + local_size - 1;

			for (size_t i = 0; i < num_indices; ++i) {
				if (indices[i] < first || indices[i] > last)
					continue;
				const size_t idx = indices[i] - first;
				this->coefs[idx] += coefs[i];
			}
		}

		void write(const std::string& filename, int id, int numboxes) const
		{
			std::ostringstream osstr;
			osstr << filename << "." << numboxes << "." << id;
			std::string full_name = osstr.str();
			std::ofstream ofs(full_name.c_str());

			for (size_t i = 0; i< local_size; ++i)
				ofs << startIndex + i << " " << coefs[i] << std::endl;
		}

		double dot(const Vector& y) const
		{
			const int n = local_size;

			assert(y.local_size >= n);

			double result = 0.0;
			for (int i = 0; i < n; ++i)
				result += coefs[i] * y.coefs[i];

			return result;
		}

		double dot2() const
		{
			double result = 0.0;
			for (int i = 0; i < local_size; ++i)
				result += coefs[i] * coefs[i];

			return result;
		}


	};


}//namespace miniFE

#endif

