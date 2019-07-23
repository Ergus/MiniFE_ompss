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

#include "Vector.hpp"
#include "Vector_functions.hpp"
#include "ElemData.hpp"
#include "mytimer.hpp"

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
			int *mat_row_cols = nullptr;
			double *mat_row_coefs = nullptr;

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

	// TODO: make this a task 
	void init_row(const int *row_coords,
	              const int *global_nodes,
	              int global_nrows,
	              const simple_mesh_description *mesh,
	              int *packed_cols,
	              double *packed_coefs)
	{
		const int ix = row_coords[0];
		const int iy = row_coords[1];
		const int iz = row_coords[2];
		int idx = 0;
		for(int sz = -1; sz <= 1; ++sz) {
			for(int sy = -1; sy <= 1; ++sy) {
				for(int sx = -1; sx <= 1; ++sx) {
					const int col_id = get_id(global_nodes[0],
					                          global_nodes[1],
					                          global_nodes[2],
					                          ix + sx, iy + sy, iz + sz);
					if (col_id >= 0 && col_id < global_nrows) {
						const int col = mesh->map_id_to_row(col_id);
						assert (col < global_nrows);
						packed_cols[idx] = col;
						packed_coefs[idx] = 0;
						++idx;
					}
				}
			}
		}

		sort_if_needed(packed_cols, idx);
	}

	// TODO: make all this weak when things works
	#pragma oss task						\
		in(row_coords[0; 3 * nrows])				\
		in(mesh[0])						\
		in(mesh_ompss2_ids_to_rows[0; mesh_ids_to_rows_size])	\
		in(row_offsets[0; nrows + 1])				\
		out(packed_cols[0; nnz])				\
		out(packed_coefs[0; nnz])
	void init_matrix_task(const int *row_coords,
	                      int global_nrows,
	                      const simple_mesh_description *mesh,
	                      const std::pair<int,int> *mesh_ompss2_ids_to_rows,
	                      size_t mesh_ids_to_rows_size,
	                      const int *row_offsets,
	                      int *packed_cols,
	                      double *packed_coefs,
	                      int nrows, int nnz)
	{
		int global_nodes[3] = {
			mesh->global_box[0][1] + 1,
			mesh->global_box[1][1] + 1,
			mesh->global_box[2][1] + 1 };

		for(int i = 0; i < nrows; ++i) {

			init_row(&row_coords[3 * i],
			         global_nodes,
			         global_nrows,
			         mesh,
			         &packed_cols[row_offsets[i]],
			         &packed_coefs[row_offsets[i]]);
		}
	}


	#pragma oss task						\
		out(row_coords[0; nrows * 3])				\
		out(rows[0; nrows])					\
		out(row_offsets[0; nrows + 1])				\
		in(mesh[0])						\
		in(mesh_ompss2_ids_to_rows[0; mesh_ids_to_rows_size])	\
		out(nnz[0])						\
		out(first_row[0])					\
		out(stop_row[0])
	void init_offsets_task(int *row_coords,
		int *rows,
		int *row_offsets,
		const simple_mesh_description *mesh,
		const std::pair<int,int> *mesh_ompss2_ids_to_rows,
		size_t mesh_ids_to_rows_size,
		size_t global_nrows,
		size_t nrows,
		size_t *nnz,
		int *first_row,
		int *stop_row)
	{
		const Box &box = mesh->extended_box;
		size_t tnnz = 0;
		size_t roffset = 0;

		int global_nodes[3] = {
			mesh->global_box[0][1] + 1,
			mesh->global_box[1][1] + 1,
			mesh->global_box[2][1] + 1 };


		for(int iz = box[2][0]; iz < box[2][1]; ++iz) {
			for(int iy = box[1][0]; iy < box[1][1]; ++iy) {
				for(int ix = box[0][0]; ix < box[0][1]; ++ix) {
					const int row_id =
						get_id(global_nodes[0],
						       global_nodes[1],
						       global_nodes[2],
						       ix, iy, iz);

					rows[roffset] = mesh->map_id_to_row(row_id);
					row_coords[roffset * 3] = ix;
					row_coords[roffset * 3 + 1] = iy;
					row_coords[roffset * 3 + 2] = iz;
					row_offsets[roffset++] = tnnz;

					for(int sz = -1; sz <= 1; ++sz) {
						for(int sy = -1; sy <= 1; ++sy) {
							for(int sx = -1; sx <= 1; ++sx) {
								const size_t col_id =
									get_id(global_nodes[0],
									       global_nodes[1],
									       global_nodes[2],
									       ix + sx,
									       iy + sy,
									       iz + sz);
								if (col_id >= 0 &&
								    col_id < global_nrows)
									++tnnz;
							}
						}
					}
				}
			}
		}
		row_offsets[nrows] = tnnz;
		*nnz = tnnz;
		*first_row = rows[0];
		*stop_row = rows[nrows - 1];
		assert(roffset == nrows);
	}



}//namespace miniFE

#endif

