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

#include "CSRMatrix.hpp"
#include "ompss_utils.hpp"

#include "singleton.hpp"

namespace miniFE
{
	struct Vector {

		int startIndex;
		size_t local_size;
		double *coefs;

		Vector() : startIndex(-1),
			   local_size(-1), coefs(nullptr)
		{}

		Vector(const Vector &in) :
			startIndex(in.startIndex),
			local_size(in.local_size),
			coefs(in.coefs)
		{}

		~Vector()
		{}

		// Arrays allocations (boxes can go in local memory)
		static void* operator new[](std::size_t sz)
		{
			void * const tmp = rrl_malloc(sz);
			dbvprintf("Calling: %s, size: %lu\n", __PRETTY_FUNCTION__, sz);
			return tmp;
		}

		static void operator delete[](void* ptr, std::size_t sz)
		{
			dbvprintf("Calling: %s, address %p size: %lu\n", __PRETTY_FUNCTION__, ptr, sz);
			return rrl_free(ptr, sz);
		}

		void init(const int startIdx, const int local_sz, double *ptr)
		{
			startIndex = startIdx;
			local_size = local_sz;
			coefs = ptr;

			#pragma oss task out(ptr[0; local_sz])
			{
				for (size_t i = 0; i < 0; ++i)
					ptr[i] = 0.0;
			}
		}

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

		void write(std::ostream &stream) const
		{
			stream << "Vector start= " << startIndex << "\n";
			print_vector("", local_size, coefs, stream);
		}

	};


	void init_vector_all(Vector *b_array,
	                     singleton *sing,
	                     size_t numboxes,
	                     const int *start, const int *length)
	{
		int *start_idx = (int *) alloca(numboxes * sizeof(int));

		int elements = 0;
		for (size_t i = 0; i < numboxes; ++i) {
			start_idx[i] = elements;
			elements += length[i];
		}

		double * tmp = sing->allocate_vectors(elements);

		assert(tmp);

		for (size_t i = 0; i < numboxes; ++i)
			b_array[i].init(start[i], length[i], &(tmp[start_idx[i]]));
	}


	void dot_task(const Vector *x, const Vector *y, double *ret)
	{
		#pragma oss task					\
			in(*x)						\
			in(x->coefs[0; x->local_size])			\
			in(*y)						\
			in(y->coefs[0; y->local_size])			\
			out(*ret)
		{
			const size_t n = x->local_size;

			assert(y->local_size >= n);

			double result = 0.0;
			for (size_t i = 0; i < n; ++i)
				result += x->coefs[i] * y->coefs[i];

			*ret = result;
		}
	}

	void dot2_task(const Vector *x, double * ret)
	{
		#pragma oss task					\
			in(*x)						\
			in(x->coefs[0; x->local_size])			\
			out(*ret)
		{
			double result = 0.0;
			for (size_t i = 0; i < x->local_size; ++i)
				result += x->coefs[i] * x->coefs[i];

			*ret = result;
		}
	}


	inline std::ostream& operator <<(std::ostream &stream, const Vector &in)
	{
		in.write(stream);
		return stream;
	}

	inline void write_task(std::string filename, const Vector *Min, size_t id)
	{
		#pragma oss task in(Min[0]) in(Min->coefs[0; Min->local_size])
		{
			std::ofstream stream;

			if (id == 0)
				stream.open(filename, std::ofstream::out);
			else
				stream.open(filename, std::ofstream::app);

			Min->write(stream);

			stream.close();
		}
	}


}//namespace miniFE

#endif

