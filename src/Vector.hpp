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
			dbvprint_vector("coefs", local_size, coefs, stream);
			stream << std::endl;
		}

	};


	void init_vector_all(Vector *b_array,
	                     singleton *sing,
	                     size_t numboxes,
	                     const int *start, const int *length)
	{

		int elements = 0;
		int *Voffsets = (int *) alloca(numboxes * sizeof(int));

		for (size_t i = 0; i < numboxes; ++i) {
			Voffsets[i] = elements;
			elements += length[i];
		}

		double * tmp = sing->allocate_vectors(elements);

		assert(tmp);

		for (size_t i = 0; i < numboxes; ++i)
			b_array[i].init(start[i], length[i], &(tmp[Voffsets[i]]));
	}


	void dot_task(const Vector *x, const Vector *y, double *ret)
	{

		double *xcoefs = x->coefs;
		size_t xlocal_size = x->local_size;
		double *ycoefs = y->coefs;
		size_t ylocal_size = y->local_size;

		#pragma oss task					\
			in(xcoefs[0; xlocal_size])			\
			in(ycoefs[0; ylocal_size])			\
			out(ret[0])
		{
			const size_t n = xlocal_size;

			assert(ylocal_size >= n);

			double result = 0.0;
			for (size_t i = 0; i < n; ++i)
				result += xcoefs[i] * ycoefs[i];

			ret[0] = result;
			dbv2printf("dot_task: %lg\n", ret[0]);
		}
	}

	void dot2_task(const Vector *x, double *ret)
	{
		double *xcoefs = x->coefs;
		size_t xlocal_size = x->local_size;

		#pragma oss task					\
			in(xcoefs[0; xlocal_size])			\
			out(ret[0; 1])
		{
			double result = 0.0;
			for (size_t i = 0; i < xlocal_size; ++i)
				result += xcoefs[i] * xcoefs[i];

			ret[0] = result;
			dbv2printf("dot2_task: %lg\n", ret[0]);
		}
	}


	inline void waxpby_task(double alpha, const Vector *x,
		double beta,  const Vector *y,
		Vector *w)
	{

		double *xcoefs = x->coefs;
		size_t xlocal_size = x->local_size;
		double *ycoefs = y->coefs;
		size_t ylocal_size = y->local_size;
		double *wcoefs = w->coefs;
		size_t wlocal_size = w->local_size;


		#pragma oss task					\
			in(xcoefs[0; xlocal_size])			\
			in(ycoefs[0; ylocal_size])			\
			out(wcoefs[0; wlocal_size])
		{
			assert(xlocal_size <= ylocal_size);
			assert(xlocal_size <= wlocal_size);
			const int n = xlocal_size;

			for (int i = 0; i < n; ++i)
				wcoefs[i] = alpha * xcoefs[i] + beta * ycoefs[i];
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

