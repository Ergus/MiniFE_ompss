#ifndef _compute_matrix_stats_hpp_
#define _compute_matrix_stats_hpp_

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

#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <outstream.hpp>
#include <utils.hpp>
#include <YAML_Doc.hpp>

namespace miniFE {

	#pragma oss task						\
		in(A_array[0; numboxes])				\
		out(global_nnz)
	void compute_matrix_stats(const CSRMatrix *A_array,
		int numboxes,
		YAML_Doc& ydoc,
		size_t &global_nnz)
{

	int min_nrows = 0, max_nrows = 0, global_nrows = 0;
	int local_nrows[numboxes];

	double dmin_nnz = 0, dmax_nnz = 0, dglobal_nnz = 0;
	double dlocal_nnz[numboxes];

	int min_box = 0, max_box = 0;

	for (int i = 0; i < numboxes; ++i) {

		// TODO: Task here
		// in A_array[i].nrows
		// in A_array[i].num_nonzeros
		// out local_nrows[i]
		local_nrows[i] = A_array[i].nrows;
		dlocal_nnz[i] = A_array[i].nnz;
	}

	// TODO taskwait here
	get_global_min_max(global_nrows,
	                   min_nrows, min_box,
	                   max_nrows, max_box,
	                   local_nrows, numboxes);

	get_global_min_max(dglobal_nnz,
	                   dmin_nnz, min_box,
	                   dmax_nnz, max_box,
	                   dlocal_nnz, numboxes);

	double avg_nrows = (double) global_nrows / numboxes;
	double avg_nnz = dglobal_nnz / numboxes;

	global_nnz = static_cast<size_t>(std::ceil(dglobal_nnz));
	size_t min_nnz = static_cast<size_t>(std::ceil(dmin_nnz));
	size_t max_nnz = static_cast<size_t>(std::ceil(dmax_nnz));
	size_t global_num_rows = global_nrows;

	ydoc.add("Matrix attributes","");
	ydoc.get("Matrix attributes")->add("Global Nrows",global_num_rows);
	ydoc.get("Matrix attributes")->add("Global NNZ",global_nnz);

	//compute how much memory the matrix occupies:
	//num-bytes = sizeof(int)*global_nrows   for A.rows
	//          + sizeof(LocalOrdinal)*global_nrows    for A.rows_offsets
	//          + sizeof(int)*global_nnz     for A.packed_cols
	//          + sizeof(Scalar)*global_nnz            for A.packed_coefs

	double invGB = 1.0 / (1024*1024*1024);
	double memGB = invGB * global_nrows * sizeof(int);
	memGB += invGB * global_nrows * sizeof(int);
	memGB += invGB * global_nnz * sizeof(int);
	memGB += invGB * global_nnz * sizeof(Scalar);
	ydoc.get("Matrix attributes")->add("Global Memory (GB)", memGB);

	size_t min_num_rows = min_nrows;
	size_t max_num_rows = max_nrows;
	ydoc.get("Matrix attributes")->add("Rows per proc MIN", min_num_rows);
	ydoc.get("Matrix attributes")->add("Rows per proc MAX", max_num_rows);
	ydoc.get("Matrix attributes")->add("Rows per proc AVG", avg_nrows);
	ydoc.get("Matrix attributes")->add("NNZ per proc MIN", min_nnz);
	ydoc.get("Matrix attributes")->add("NNZ per proc MAX", max_nnz);
	ydoc.get("Matrix attributes")->add("NNZ per proc AVG", avg_nnz);

	global_nnz;
}

}//namespace miniFE

#endif

