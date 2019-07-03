#ifndef _SparseMatrix_functions_hpp_
#define _SparseMatrix_functions_hpp_

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
#include <vector>
#include <set>
#include <algorithm>
#include <sstream>
#include <fstream>

#include <Vector.hpp>
#include <Vector_functions.hpp>
#include <ElemData.hpp>
#include <exchange_externals.hpp>
#include <mytimer.hpp>

#ifdef MINIFE_HAVE_TBB
#include <LockingMatrix.hpp>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace miniFE
{

	template<typename T, typename U>
	void sort_with_companions(ptrdiff_t len, T* array, U* companions)
	{
		ptrdiff_t i, j, index;
		U companion;

		for (i = 1; i < len; i++) {
			index = array[i];
			companion = companions[i];
			j = i;
			while ((j > 0) && (array[j-1] > index)) {
				array[j] = array[j-1];
				companions[j] = companions[j-1];
				j = j - 1;
			}
			array[j] = index;
			companions[j] = companion;
		}
	}

	inline void sum_into_row(int row_len, int *row_indices, double *row_coefs,
	                  size_t num_inputs, const int *input_indices, const double *input_coefs)
	{
		for (size_t i = 0; i < num_inputs; ++i) {
			int *loc = std::lower_bound(row_indices, row_indices + row_len, input_indices[i]);

			if (loc - row_indices < row_len && *loc == input_indices[i])
				row_coefs[loc - row_indices] += input_coefs[i];
		}
	}

	template<typename MatrixType>
	void sum_in_symm_elem_matrix(size_t num, int *indices, double *coefs, MatrixType &mat)
	{
		//indices is length num (which should be nodes-per-elem)
		//coefs is the upper triangle of the element diffusion matrix
		//which should be length num*(num+1)/2
		//std::cout<<std::endl;

		int row_offset = 0;

		for(size_t i = 0; i < num; ++i) {
			int row = indices[i];

			const double *row_coefs = &coefs[row_offset];
			const int *row_col_inds = &indices[i];
			size_t row_len = num - i;
			row_offset += row_len;

			size_t mat_row_len = 0;
			int *mat_row_cols = NULL;
			double *mat_row_coefs = NULL;

			mat.get_row_pointers(row, mat_row_len, mat_row_cols, mat_row_coefs);
			if (mat_row_len == 0)
				continue;

			sum_into_row(mat_row_len, mat_row_cols, mat_row_coefs,
			             row_len, row_col_inds, row_coefs);

			int offset = i;
			for(size_t j = 0; j < i; ++j) {
				double coef = coefs[offset];
				sum_into_row(mat_row_len, mat_row_cols, mat_row_coefs,
				             1, &indices[j], &coef);
				offset += num - (j+1);
			}
		}
	}

	template<typename MatrixType>
	void sum_in_elem_matrix(size_t num, const int *indices, const double *coefs,
	                        MatrixType& mat)
	{
		size_t offset = 0;

		for(size_t i = 0; i < num; ++i) {
			sum_into_row(indices[i], num,
			             &indices[0], &coefs[offset], mat);
			offset += num;
		}
	}

	template<typename MatrixType>
	void sum_into_global_linear_system(ElemData &elem_data, MatrixType &A, Vector &b)
	{
		sum_in_symm_elem_matrix(elem_data.nodes_per_elem,
		                        elem_data.elem_node_ids,
		                        elem_data.elem_diffusion_matrix, A);
		b.sum_into_vector(elem_data.nodes_per_elem,
		                  elem_data.elem_node_ids,
		                  elem_data.elem_source_vector);
	}

	#ifdef MINIFE_HAVE_TBB
	template<typename MatrixType>
	void sum_in_elem_matrix(size_t num, const int *indices, const double* coefs,
	                        LockingMatrix<MatrixType>& mat)
	{
		size_t offset = 0;

		for(size_t i=0; i<num; ++i) {
			mat.sum_in(indices[i], num, &indices[0], &coefs[offset]);
			offset += num;
		}
	}

	template<typename GlobalOrdinal, typename Scalar,
	         typename MatrixType, typename VectorType>
	void
	sum_into_global_linear_system(ElemData<GlobalOrdinal,Scalar>& elem_data,
	                              LockingMatrix<MatrixType>& A, LockingVector<VectorType>& b)
	{
		sum_in_elem_matrix(elem_data.nodes_per_elem, elem_data.elem_node_ids,
		                   elem_data.elem_diffusion_matrix, A);
		sum_into_vector(elem_data.nodes_per_elem, elem_data.elem_node_ids,
		                elem_data.elem_source_vector, b);
	}
	#endif

	template<typename MatrixType>
	void
	add_to_diagonal(typename MatrixType::ScalarType value, MatrixType& mat)
	{
		for(size_t i=0; i<mat.rows.size(); ++i)
			sum_into_row(mat.rows[i], 1, &mat.rows[i], &value, mat);
	}


	template<typename MatrixType>
	void rearrange_matrix_local_external(MatrixType& A)
	{
		//This function will rearrange A so that local entries are contiguous at the front
		//of A's memory, and external entries are contiguous at the back of A's memory.
		//
		//A.row_offsets will describe where the local entries occur, and
		//A.row_offsets_external will describe where the external entries occur.

		typedef typename MatrixType::GlobalOrdinalType GlobalOrdinal;
		typedef typename MatrixType::LocalOrdinalType LocalOrdinal;
		typedef typename MatrixType::ScalarType Scalar;

		size_t nrows = A.rows.size();
		std::vector<LocalOrdinal> tmp_row_offsets(nrows*2);
		std::vector<LocalOrdinal> tmp_row_offsets_external(nrows*2);

		LocalOrdinal num_local_nz = 0;
		LocalOrdinal num_extern_nz = 0;

		//First sort within each row of A, so that local entries come
		//before external entries within each row.
		//tmp_row_offsets describe the locations of the local entries, and
		//tmp_row_offsets_external describe the locations of the external entries.
		//
		for(size_t i=0; i<nrows; ++i) {
			GlobalOrdinal* row_begin = &A.packed_cols[A.row_offsets[i]];
			GlobalOrdinal* row_end = &A.packed_cols[A.row_offsets[i+1]];

			Scalar* coef_row_begin = &A.packed_coefs[A.row_offsets[i]];

			tmp_row_offsets[i*2] = A.row_offsets[i];
			tmp_row_offsets[i*2+1] = A.row_offsets[i+1];
			tmp_row_offsets_external[i*2] = A.row_offsets[i+1];
			tmp_row_offsets_external[i*2+1] = A.row_offsets[i+1];

			ptrdiff_t row_len = row_end - row_begin;

			sort_with_companions(row_len, row_begin, coef_row_begin);

			GlobalOrdinal* row_iter = std::lower_bound(row_begin, row_end, nrows);

			LocalOrdinal offset = A.row_offsets[i] + row_iter-row_begin;
			tmp_row_offsets[i*2+1] = offset;
			tmp_row_offsets_external[i*2] = offset;

			num_local_nz += tmp_row_offsets[i*2+1]-tmp_row_offsets[i*2];
			num_extern_nz += tmp_row_offsets_external[i*2+1]-tmp_row_offsets_external[i*2];
		}

		//Next, copy the external entries into separate arrays.

		std::vector<GlobalOrdinal> ext_cols(num_extern_nz);
		std::vector<Scalar> ext_coefs(num_extern_nz);
		std::vector<LocalOrdinal> ext_offsets(nrows+1);
		LocalOrdinal offset = 0;
		for(size_t i=0; i<nrows; ++i) {
			ext_offsets[i] = offset;
			for(LocalOrdinal j=tmp_row_offsets_external[i*2];
			    j<tmp_row_offsets_external[i*2+1]; ++j) {
				ext_cols[offset] = A.packed_cols[j];
				ext_coefs[offset++] = A.packed_coefs[j];
			}
		}
		ext_offsets[nrows] = offset;

		//Now slide all local entries down to the beginning of A's packed arrays

		A.row_offsets.resize(nrows+1);
		offset = 0;
		for(size_t i=0; i<nrows; ++i) {
			A.row_offsets[i] = offset;
			for(LocalOrdinal j=tmp_row_offsets[i*2]; j<tmp_row_offsets[i*2+1]; ++j) {
				A.packed_cols[offset] = A.packed_cols[j];
				A.packed_coefs[offset++] = A.packed_coefs[j];
			}
		}
		A.row_offsets[nrows] = offset;

		//Finally, copy the external entries back into A.packed_cols and
		//A.packed_coefs, starting at the end of the local entries.

		for(LocalOrdinal i=offset; i<offset+ext_cols.size(); ++i) {
			A.packed_cols[i] = ext_cols[i-offset];
			A.packed_coefs[i] = ext_coefs[i-offset];
		}

		A.row_offsets_external.resize(nrows+1);
		for(size_t i = 0; i <= nrows; ++i)
			A.row_offsets_external[i] = ext_offsets[i] + offset;
	}

	//------------------------------------------------------------------------
	template<typename MatrixType>
	void zero_row_and_put_1_on_diagonal(MatrixType *A, const int row)
	{
		size_t row_len = 0;
		int* cols = NULL;
		double *coefs = NULL;
		A->get_row_pointers(row, row_len, cols, coefs);

		for(size_t i = 0; i < row_len; ++i)
			coefs[i] = (cols[i] == row);
	}

	//------------------------------------------------------------------------
	template<typename MatrixType>
	void impose_dirichlet(double prescribed_value,
	                      MatrixType *A, Vector *b,
	                      int global_nx, int global_ny, int global_nz,
	                      const int *bc_rows_array, size_t bc_rows_size)
	{
		const int first_local_row = A->nrows > 0 ? A->rows[0] : 0;
		const int last_local_row  = A->nrows > 0 ? A->rows[A->nrows - 1] : -1;

		for (size_t i = 0; i < bc_rows_size; ++i) {
			const int row = bc_rows_array[i];

			if (row >= first_local_row && row <= last_local_row) {
				const size_t local_row = row - first_local_row;
				b->coefs[local_row] = prescribed_value;
				zero_row_and_put_1_on_diagonal(A, row);
			}
		}

		for (size_t i = 0; i < A->nrows; ++i) {
			const int row = A->rows[i];

			if (std::binary_search(bc_rows_array,
			                       &bc_rows_array[bc_rows_size],
			                       row))
				continue;

			size_t row_length = 0;
			int *cols = NULL;
			double *coefs = NULL;
			A->get_row_pointers(row, row_length, cols, coefs);

			double sum = 0.0;
			for(size_t j = 0; j < row_length; ++j) {
				if (std::binary_search(bc_rows_array,
				                       &bc_rows_array[bc_rows_size],
				                       cols[j])) {
					sum += coefs[j];
					coefs[j] = 0.0;
				}
			}

			b->coefs[i] -= sum * prescribed_value;
		}
	}

}//namespace miniFE

#endif

