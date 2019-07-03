#ifndef _CSRMatrix_hpp_
#define _CSRMatrix_hpp_

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
#include <algorithm>

#include "SparseMatrix_functions.hpp"
#include "simple_mesh_description.hpp"

#include "utils.hpp"
#include "ompss_utils.hpp"

namespace miniFE
{

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

	// TODO: make this weak
	void init_matrix_task(const int *row_coords,
		const int *global_nodes,
		int global_nrows,
		const simple_mesh_description *mesh,
		const int *row_offsets,
		int *packed_cols,
		double *packed_coefs,
		int nrows, int nnz)
	{
		#pragma oss task			\
			in(row_coords[0; 3 * nrows])	\
			in(global_nodes[0; 3])		\
			in(global_nrows[0; nrows])	\
			in(*mesh)			\
			in(row_offsets[0; nrows + 1])	\
			out(packed_cols[0; nnz])		\
			out(packed_coefs[0; nnz])
		{
			for(int i = 0; i < nrows; ++i) {

				const int offset = row_offsets[i];
				const int next_offset = row_offsets[i + 1];

				#pragma oss task			\
					in(row_coords[3 * i; 3])	\
					in(global_nodes[0; 3]		\
					in(*mesh)			\
					in(mesh->ompss2_ids_to_rows[0; mesh->ids_to_rows_size] \
					out(packed_cols[offset: next_offset]) \
					out(packed_coefs[offset: next_offset])
				{
					init_row(&row_coords[3 * i],
					         global_nodes,
					         global_nrows,
					         mesh,
					         &packed_cols[offset],
					         &packed_coefs[offset]);
				}
			}
		}
	}


	void init_offsets_task(int *row_coords,
	                       int *rows,
	                       int *row_offsets,
	                       int *global_nodes,
	                       const simple_mesh_description *mesh,
	                       size_t global_nrows,
	                       size_t nrows,
	                       size_t &nnz)
	{
		#pragma oss task \
			out(row_coords[0; nrows * 3])			\
			out(rows[0; nrows])				\
			out(row_offsets[0; tnrows + 1])			\
			in(*box)					\
			in(global_nodes[0; 3])				\
			in(mesh[0])					\
			in(mesh_i[0].ompss2_bc_rows_0[0; mesh_i[0].bc_rows_0_size]) \
			in(mesh_i[0].ompss2_bc_rows_1[0; mesh_i[0].bc_rows_1_size]) \
			in(mesh_i[0].ompss2_ids_to_rows[0; mesh_i[0].ids_to_rows_size]) \
			in(tnrows)					\
			out(nnz)

		{
			const Box &box = mesh->extended_box;
			size_t tnnz = 0;
			size_t roffset = 0;

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
			nnz = tnnz;
			assert(roffset == nrows);

		}
		#pragma oss taskwait
	}

	class CSRMatrix {
	public:
		bool has_local_indices;
		int *rows;
		int *row_offsets;
		//int *row_offsets_external;
		size_t global_nrows, nrows, nnz;

		int *packed_cols;
		double *packed_coefs;

		int num_cols;

		int nrecv_neighbors;
		int nexternals;
		int *recv_neighbors;
		double **recv_ptr;   // List of remote pointers
		int *recv_length;
		int *external_index;

		int nsend_neighbors;
		int *send_neighbors;
		int *send_length;
		int nelements_to_send;
		int *elements_to_send;
		double *send_buffer;

		CSRMatrix()
			: has_local_indices(false),
			  rows(), row_offsets(nullptr), //row_offsets_external(),
			  global_nrows(-1), nrows(-1),
			  nnz(0), packed_cols(nullptr), packed_coefs(nullptr),
			  num_cols(0),

			  nrecv_neighbors(0), nexternals(0), recv_neighbors(nullptr),
			  recv_ptr(nullptr), recv_length(nullptr), external_index(nullptr),

			  nsend_neighbors(0), send_neighbors(nullptr), send_length(nullptr),
			  nelements_to_send(0), elements_to_send(nullptr), send_buffer(nullptr)

		{
		}

		~CSRMatrix()
		{
			// generate_matrix_structure
			rrd_free(rows, nrows * sizeof(int));
			rrd_free(row_offsets, (nrows + 1) * sizeof(int));
			//rrd_free(row_offsets_external, (nrows + 1) * sizeof(int));
			rrd_free(packed_cols, nnz * sizeof(int));
			rrd_free(packed_coefs, nnz * sizeof(double));

			// make_local_matrix
			// Receive
			rrd_free(recv_neighbors, nrecv_neighbors * sizeof(int));
			rrd_free(recv_ptr, nrecv_neighbors * sizeof(double *));

			rrd_free(recv_length, nrecv_neighbors * sizeof(int));
			rrd_free(external_index, nexternals * sizeof(int));

			// Send
			rrd_free(send_neighbors, nsend_neighbors * sizeof(int));
			rrd_free(send_length, nsend_neighbors * sizeof(int));
			rrd_free(elements_to_send, nelements_to_send * sizeof(int));
			rrd_free(send_buffer, nelements_to_send * sizeof(double));
		}


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



		void get_row_pointers(int row, size_t &row_length, int *&cols, double *&coefs) const
		{
			ptrdiff_t local_row = -1;
			//first see if we can get the local-row index using fast direct lookup:
			if (nrows > 0) {
				int idx = row - rows[0];
				if (idx < nrows && rows[idx] == row)
					local_row = idx;
			}

			//if we didn't get the local-row index using direct lookup, try a
			//more expensive binary-search:
			if (local_row == -1) {
				auto row_iter = std::lower_bound(rows, &rows[nrows], row);

				//if we still haven't found row, it's not local so jump out:
				if (row_iter == &rows[nrows] || *row_iter != row) {
					row_length = 0;
					return;
				}

				local_row = row_iter - rows;
			}

			const int offset = row_offsets[local_row];
			row_length = row_offsets[local_row + 1] - offset;
			cols = &packed_cols[offset];
			coefs = &packed_coefs[offset];
		}

		void generate_matrix_structure(const simple_mesh_description *mesh)
		{

			int global_nodes[3] = {
				mesh->global_box[0][1] + 1,
				mesh->global_box[1][1] + 1,
				mesh->global_box[2][1] + 1 };

			nrows = mesh->extended_box.get_num_ids();

			//num-owned-nodes in each dimension is num-elems+1
			//only if num-elems > 0 in that dimension *and*
			//we are at the high end of the global range in that dimension:
			global_nrows = global_nodes[0] * global_nodes[1] * global_nodes[2];

			rows = (int *) rrd_malloc(nrows * sizeof(int));
			row_offsets = (int *) rrd_malloc((nrows + 1) * sizeof(int));

			int *row_coords = (int *) rrl_malloc(nrows * 3 * sizeof(int));

			init_offsets_task(row_coords,
			                  rows,
			                  row_offsets,
			                  global_nodes,
			                  mesh,
			                  global_nrows,
			                  nrows,
			                  nnz);

			packed_cols = (int *) rrd_malloc(nnz * sizeof(int));
			packed_coefs = (double *) rrd_malloc(nnz * sizeof(double));

			 init_matrix_task(row_coords,
			                 global_nodes,
			                 global_nrows,
			                 mesh,
			                 row_offsets,
			                 packed_cols,
			                 packed_coefs,
			                 nrows,
			                 nnz);

			// TODO: Taskwait here
			rrl_free(row_coords, nrows * 3 * sizeof(int));

		}


		void sum_into_row(int row, size_t num_indices, int *col_inds, double *coefs)
		{
			size_t row_len = 0;
			int *mat_row_cols = NULL;
			double *mat_row_coefs = NULL;

			get_row_pointers(row, row_len, mat_row_cols, mat_row_coefs);
			if (row_len == 0)
				return;

			miniFE::sum_into_row(row_len, mat_row_cols, mat_row_coefs, num_indices, col_inds, coefs);
		}


		void write(const std::string &filename, int id, int numboxes) const
		{
			std::ostringstream osstr;

			std::string full_name = osstr.str();
			std::ofstream ofs(full_name.c_str());

			if (id == 0)
				ofs << nrows << " " << nnz << std::endl;

			for (size_t i = 0; i < nrows; ++i) {
				size_t row_len = 0;
				int *cols = NULL;
				double *coefs = NULL;
				get_row_pointers(rows[i], row_len, cols, coefs);

				for (size_t j = 0; j < row_len; ++j)
					ofs << rows[i] << " " << cols[j] << " " << coefs[j] << std::endl;

			}
		}

		void matvec(const Vector &x, Vector &y)
		{
			double beta = 0.0;

			for (size_t row = 0; row < nrows; ++row) {
				double sum = beta * y.coefs[row];

				for(int i = row_offsets[row]; i < row_offsets[row + 1]; ++i) {
					const int col = packed_cols[i];
					sum += packed_coefs[i] * x.coefs[col];
				}

				//std::cout << "row[" << row << "] = " << sum << std::endl;
				y.coefs[row] = sum;
			}
		}
	};

}//namespace miniFE

#endif

