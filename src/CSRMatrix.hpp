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

	class CSRMatrix {
	private:
		void init_row(int i,
		              const int *row_coords_vec,
		              int global_nx, int global_ny, int global_nz,
		              int global_nrows,
		              const simple_mesh_description &mesh)
		{
			const int offset = row_offsets[i];
			int ix = row_coords_vec[i * 3];
			int iy = row_coords_vec[i * 3 + 1];
			int iz = row_coords_vec[i * 3 + 2];
			int nnz = 0;
			for(int sz = -1; sz <= 1; ++sz) {
				for(int sy = -1; sy <= 1; ++sy) {
					for(int sx = -1; sx <= 1; ++sx) {
						const int col_id = get_id(global_nx, global_ny, global_nz,
						                          ix + sx, iy + sy, iz + sz);
						if (col_id >= 0 && col_id < global_nrows) {
							const int col = mesh.map_id_to_row(col_id);
							assert (col < global_nrows);
							packed_cols[offset + nnz] = col;
							packed_coefs[offset + nnz] = 0;
							++nnz;
						}
					}
				}
			}

			sort_if_needed(&packed_cols[offset], nnz);
		}

		void init_matrix(const int *row_coords_vec,
		                 int global_nx, int global_ny, int global_nz,
		                 int global_nrows,
		                 const simple_mesh_description &mesh)
		{
			printf("nrows %d\n",nrows);
			for (int i = 0; i < nrows; ++i)
				init_row(i, row_coords_vec,
				         global_nx, global_ny, global_nz,
				         global_nrows, mesh);
		}

	public:
		bool has_local_indices;
		int *rows;
		int *row_offsets;
		int *row_offsets_external;
		int global_nrows, nrows;

		int nnz;
		int *packed_cols;
		double *packed_coefs;

		int num_cols;
		int num_nonzeros;

		int nrecv_neighbors;
		int nexternals;
		int *recv_neighbors;
		double **recv_ptr;   // List of remote pointers
		int *recv_length;
		int *external_index;

		int nsend_neighbors;
		int nelements_to_send;
		int *send_neighbors;
		int *send_length;
		int *elements_to_send;
		double *send_buffer;

		CSRMatrix()
			: has_local_indices(false),
			  rows(), row_offsets(nullptr), row_offsets_external(),
			  nnz(0), packed_cols(nullptr), packed_coefs(nullptr),
			  num_cols(0), global_nrows(-1), nrows(-1), num_nonzeros(0),

			  nrecv_neighbors(0), recv_neighbors(nullptr), recv_length(nullptr),
			  nexternals(0), external_index(nullptr),

			  nsend_neighbors(0), send_neighbors(nullptr), send_length(nullptr),
			  nelements_to_send(0), elements_to_send(nullptr), send_buffer(nullptr)

		{
		}

		~CSRMatrix()
		{
			// generate_matrix_structure
			rrd_free(rows, nrows * sizeof(int));
			rrd_free(row_offsets, (nrows + 1) * sizeof(int));
			rrd_free(row_offsets_external, (nrows + 1) * sizeof(int));
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
			dbprintf("Calling: %s, size: %lu\n", __PRETTY_FUNCTION__, sz);
			return tmp;
		}

		static void operator delete[](void* ptr, std::size_t sz)
		{
			printf("Calling: %s, address %p size: %lu\n", __PRETTY_FUNCTION__, ptr, sz);
			return rrl_free(ptr, sz);
		}



		void get_row_pointers(int row, size_t &row_length, int *&cols, double *&coefs) const
		{
			ptrdiff_t local_row = -1;
			//first see if we can get the local-row index using fast direct lookup:
			if (nrows > 0) {
				ptrdiff_t idx = row - rows[0];
				if (idx < nrows && rows[idx] == row) {
					local_row = idx;
				}
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

		void generate_matrix_structure(const simple_mesh_description &mesh)
		{
			int threw_exc = 0;
			try {

				const int global_nodes_x = mesh.global_box[0][1] + 1;
				const int global_nodes_y = mesh.global_box[1][1] + 1;
				const int global_nodes_z = mesh.global_box[2][1] + 1;

				Box box(mesh.local_box);

				//num-owned-nodes in each dimension is num-elems+1
				//only if num-elems > 0 in that dimension *and*
				//we are at the high end of the global range in that dimension:
				if (box[0][1] > box[0][0] && box[0][1] == mesh.global_box[0][1])
					++box[0][1];
				if (box[1][1] > box[1][0] && box[1][1] == mesh.global_box[1][1])
					++box[1][1];
				if (box[2][1] > box[2][0] && box[2][1] == mesh.global_box[2][1])
					++box[2][1];

				global_nrows = global_nodes_x * global_nodes_y * global_nodes_z;
				nrows = box.get_num_ids();

				try {
					// This substituted the reserve function
					rows = (int *) rrd_malloc(nrows * sizeof(int));
					row_offsets = (int *) rrd_malloc((nrows + 1) * sizeof(int));
					row_offsets_external = (int *) rrd_malloc((nrows + 1) * sizeof(int));
				} catch(std::exception& exc) {
					std::ostringstream osstr;
					osstr << "One of rows.resize, row_offsets.resize, "
					      << "packed_cols.reserve or packed_coefs.reserve: nrows = "
					      << nrows <<": ";
					osstr << exc.what();
					std::string str1 = osstr.str();
					throw std::runtime_error(str1);
				}

				int *row_coords = (int *) rrl_malloc(nrows * 3 * sizeof(int));

				unsigned roffset = 0;
				nnz = 0;

				// TODO: Task here to touch the array row_offsets
				// out row_coords[0;]
				{
					for(int iz = box[2][0]; iz < box[2][1]; ++iz) {
						for(int iy = box[1][0]; iy < box[1][1]; ++iy) {
							for(int ix = box[0][0]; ix < box[0][1]; ++ix) {
								int row_id =
									get_id(global_nodes_x, global_nodes_y, global_nodes_z,
									       ix, iy, iz);

								rows[roffset] = mesh.map_id_to_row(row_id);
								row_coords[roffset * 3] = ix;
								row_coords[roffset * 3 + 1] = iy;
								row_coords[roffset * 3 + 2] = iz;
								row_offsets[roffset++] = nnz;

								for(int sz = -1; sz <= 1; ++sz) {
									for(int sy = -1; sy <= 1; ++sy) {
										for(int sx = -1; sx <= 1; ++sx) {
											const int col_id = get_id(global_nodes_x,
											                          global_nodes_y,
											                          global_nodes_z,
											                          ix + sx, iy + sy, iz + sz);
											if (col_id >= 0 && col_id < global_nrows)
												++nnz;
										}
									}
								}
							}
						}
					}
					row_offsets[nrows] = nnz;
				}

				assert(roffset == nrows);
				num_nonzeros = nnz;

				packed_cols = (int *) rrd_malloc(nnz * sizeof(int));
				packed_coefs = (double *) rrd_malloc(nnz * sizeof(double));

				// TODO: Task Here
				init_matrix(row_coords,
				            global_nodes_x, global_nodes_y, global_nodes_z,
				            global_nrows, mesh);

				// TODO: Taskwait here
				rrl_free(row_coords, nrows * 3 * sizeof(int));

			} catch(...) {
				std::cout << " threw an exception in generate_matrix_structure,"
				          << " probably due to running out of memory."
				          << std::endl;
				threw_exc = 1;
			}
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
				ofs << nrows << " " << num_nonzeros << std::endl;

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

			for (int row = 0; row < nrows; ++row) {
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

