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
#include "perform_element_loop.hpp"

#include "utils.hpp"
#include "ompss_utils.hpp"

#include <fstream>      // std::ofstream

namespace miniFE
{

	class CSRMatrix {
	public:
		int id;
		bool has_local_indices, is_copy;
		int *rows;
		int *row_offsets;
		//int *row_offsets_external;
		size_t global_nrows, nrows, nnz;

		int *packed_cols;
		double *packed_coefs;

		int first_row;
		int num_cols;

		// The send elements are not released here because they are in
		// the singleton
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
			: has_local_indices(false), is_copy(false),
			  rows(nullptr), row_offsets(nullptr), //row_offsets_external(),
			  global_nrows(-1), nrows(-1),
			  nnz(0), packed_cols(nullptr), packed_coefs(nullptr),
			  first_row(-1), num_cols(0),

			  nrecv_neighbors(0), nexternals(0), recv_neighbors(nullptr),
			  recv_ptr(nullptr), recv_length(nullptr), external_index(nullptr),

			  nsend_neighbors(0), send_neighbors(nullptr), send_length(nullptr),
			  nelements_to_send(0), elements_to_send(nullptr), send_buffer(nullptr)

		{}

		CSRMatrix(const CSRMatrix &in)
			: has_local_indices(in.has_local_indices),
			  is_copy(true),
			  rows(in.rows),
			  row_offsets(in.row_offsets),
			  global_nrows(in.global_nrows),
			  nrows(in.nrows),
			  nnz(in.nnz),
			  packed_cols(in.packed_cols),
			  packed_coefs(in.packed_coefs),
			  first_row(in.first_row),
			  num_cols(in.num_cols),

			  nrecv_neighbors(in.nrecv_neighbors),
			  nexternals(in.nexternals),
			  recv_neighbors(in.recv_neighbors),
			  recv_ptr(in.recv_ptr),
			  recv_length(in.recv_length),
			  external_index(in.external_index),

			  nsend_neighbors(in.nsend_neighbors),
			  send_neighbors(in.send_neighbors),
			  send_length(in.send_length),
			  nelements_to_send(in.nelements_to_send),
			  elements_to_send(in.elements_to_send),
			  send_buffer(in.send_buffer)

		{}

		~CSRMatrix()
		{
			if (!is_copy) {
				// generate_matrix_structure
				assert(nrows > 0);
				rrd_free(rows, nrows * sizeof(int));
				rrd_free(row_offsets, (nrows + 1) * sizeof(int));

				//rrd_free(row_offsets_external, (nrows + 1) * sizeof(int));
				assert(nnz > 0);
				rrd_free(packed_cols, nnz * sizeof(int));
				rrd_free(packed_coefs, nnz * sizeof(double));
			}
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
				if (idx < (int)nrows && rows[idx] == row)
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

		void write(std::ostream &stream) const
		{
			stream << "CSRMatrix: " << id << "\n";
			stream << "has_local_indices" << "="<< has_local_indices << "\n";
			stream << "global_nrows" << "="<< global_nrows << "\n";
			stream << "nrows" << "="<< nrows << "\n";
			stream << "nnz" << "="<< nnz << "\n";
			stream << "first_row" << "="<< first_row << "\n";
			stream << "num_cols" << "="<< num_cols << "\n";
			stream << "nrecv_neighbors" << "="<< nrecv_neighbors << "\n";
			stream << "nexternals" << "="<< nexternals << "\n";
			stream << "nsend_neighbors" << "="<< nsend_neighbors << "\n";
			stream << "nelements_to_send" << "="<< nelements_to_send << "\n";

			#ifdef VERBOSE

			for (size_t i = 0; i < nrows; ++i) {
				int *Acols = NULL;
				double *Acoefs = NULL;
				size_t row_len = 0;
				get_row_pointers(rows[i], row_len, Acols, Acoefs);

				stream << i << ":" << rows[i] << ":" << row_offsets[i] << "={";

				for (size_t j = 0; j < row_len; ++j) {
					if (j > 0)
						stream << "; ";
					stream << "<" << Acols[j] << ";" << Acoefs[j] << ">";
				}
				stream << "}\n";
			}

			print_vector("recv_neighbors", nrecv_neighbors, recv_neighbors, stream);
			print_vector("recv_length", nrecv_neighbors, recv_length, stream);
			print_vector("externals", nexternals, external_index, stream);

			print_vector("send_neighbors", nsend_neighbors, send_neighbors, stream);
			print_vector("send_length", nsend_neighbors, send_length, stream);
			print_vector("nelements_to_send", nelements_to_send, elements_to_send, stream);

			#endif

		}


	}; // class CSRMatrix


	#pragma oss task						\
		inout(A[0])						\
		in(mesh[0])						\
		in(mesh_ompss2_ids_to_rows[0; mesh_ids_to_rows_size])       // weak
	void generate_matrix_structure_task(CSRMatrix *A,
	                                    const simple_mesh_description *mesh,
	                                    std::pair<int,int> *mesh_ompss2_ids_to_rows,
	                                    size_t mesh_ids_to_rows_size,
	                                    size_t _id)
	{

		#ifdef VERBOSE
		std::string filename = "VERB_mesh_gen_mat_" + std::to_string(_id) + ".verb";
		std::ofstream stream(filename);
		mesh->write(stream);
		stream.close();
		#endif

		int global_nodes[3] = {
			mesh->global_box[0][1] + 1,
			mesh->global_box[1][1] + 1,
			mesh->global_box[2][1] + 1 };

		A->id = _id;
		A->nrows = mesh->extended_box.get_num_ids();

		//num-owned-nodes in each dimension is num-elems+1
		//only if num-elems > 0 in that dimension *and*
		//we are at the high end of the global range in that dimension:
		A->global_nrows = global_nodes[0] * global_nodes[1] * global_nodes[2];

		A->rows = (int *) rrd_malloc(A->nrows * sizeof(int));
		A->row_offsets = (int *) rrd_malloc((A->nrows + 1) * sizeof(int));

		int *row_coords = (int *) rrl_malloc(A->nrows * 3 * sizeof(int));

		init_offsets_task(row_coords,
		                  A->rows,
		                  A->row_offsets,
		                  global_nodes,
		                  mesh,
		                  mesh->ompss2_ids_to_rows,
		                  mesh->ids_to_rows_size,
		                  A->global_nrows,
		                  A->nrows,
		                  &(A->nnz),
		                  &(A->first_row));


		#pragma oss taskwait
		A->packed_cols = (int *) rrd_malloc(A->nnz * sizeof(int));
		A->packed_coefs = (double *) rrd_malloc(A->nnz * sizeof(double));

		init_matrix_task(row_coords,
		                 global_nodes,
		                 A->global_nrows,
		                 mesh,
		                 mesh->ompss2_ids_to_rows,
		                 mesh->ids_to_rows_size,
		                 A->row_offsets,
		                 A->packed_cols,
		                 A->packed_coefs,
		                 A->nrows,
		                 A->nnz);

		// TODO: Taskwait here
		rrl_free(row_coords, A->nrows * 3 * sizeof(int));

	}


	void matvec_task(const CSRMatrix *A, const Vector *x, Vector *y)
	{
		int *Arow_offsets = A->row_offsets;
		size_t Anrows = A->nrows;
		int *Apacked_cols = A->packed_cols;
		double *Apacked_coefs = A->packed_coefs;
		double *xcoefs = x->coefs;
		size_t xlocal_size = x->local_size;
		double *ycoefs = y->coefs;
		size_t ylocal_size = y->local_size;
		size_t Annz = A->nnz;

		assert(xlocal_size >= Anrows);
		assert(ylocal_size >= Anrows);

		#pragma oss task				\
			in(*A)					\
			in(Arow_offsets[0; Anrows + 1])	\
			in(Apacked_cols[0; Annz])		\
			in(Apacked_coefs[0; Annz])		\
			in(*x)					\
			in(xcoefs[0; xlocal_size])		\
			in(*y)					\
			inout(ycoefs[0; ylocal_size])
		{
			const double beta = 0.0;  // I don't really understand what is this for

			for (size_t row = 0; row < Anrows; ++row) {
				double sum = beta * ycoefs[row];

				for(int i = Arow_offsets[row]; i < Arow_offsets[row + 1]; ++i) {
					const int col = Apacked_cols[i];
					assert((size_t)i < Annz);
					sum += Apacked_coefs[i] * xcoefs[col];
				}

				//std::cout << "row[" << row << "] = " << sum << std::endl;
				y->coefs[row] = sum;
			}
		}
	}


	inline void write_task(std::string filename, const CSRMatrix &Min, size_t id)
	{
		CSRMatrix Mcopy(Min);  // This is a work around for the dependency issue

		#pragma oss task					\
			in(Mcopy)						\
			in(Mcopy.rows[0; Mcopy.nrows])			\
			in(Mcopy.row_offsets[0; Mcopy.nrows + 1])	\
			in(Mcopy.packed_cols[0; Mcopy.nnz])		\
			in(Mcopy.packed_coefs[0; Mcopy.nnz])		\
			in(Mcopy.recv_neighbors[0; Mcopy.nrecv_neighbors]) \
			in(Mcopy.recv_length[0; Mcopy.nrecv_neighbors])	\
			in(Mcopy.external_index[0; Mcopy.nexternals])	\
			in(Mcopy.send_neighbors[0; Mcopy.nsend_neighbors]) \
			in(Mcopy.send_length[0; Mcopy.nsend_neighbors])	\
			in(Mcopy.elements_to_send[0; Mcopy.nelements_to_send])
		{
			std::ofstream stream;

			if (id == 0)
				stream.open(filename, std::ofstream::out);
			else
				stream.open(filename, std::ofstream::app);

			Mcopy.write(stream);

			stream.close();
		}

		#pragma oss taskwait
	}


	#pragma oss task						\
		in(*mesh)		 				\
		in(mesh_ompss2_ids_to_rows[0; mesh_i_ids_to_rows_size])	\
		inout(*A)						\
		inout(A_rows[0; A_nrows])				\
		inout(A_row_offsets[0; A_nrows + 1])			\
		inout(A_packed_cols[0; A_nnz])				\
		inout(A_packed_coefs[0; A_nnz])				\
		inout(*b)						\
		inout(b_coefs[0; b_local_size])
	inline void assemble_FE_data_task(
		const simple_mesh_description *mesh,
		const std::pair<int,int> *mesh_ompss2_ids_to_rows,
		size_t mesh_i_ids_to_rows_size,
		CSRMatrix *A,
		int *A_rows,
		int *A_row_offsets,
		size_t A_nrows,
		int *A_packed_cols,
		double *A_packed_coefs,
		size_t A_nnz,
		Vector *b,
		double *b_coefs,
		size_t b_local_size)
	{
		Box local_elem_box(mesh->local_box);

		if (local_elem_box.get_num_ids() < 1)
			return;

		//
		//We want the element-loop to loop over our (processor-local) domain plus a
		//ghost layer, so we can assemble the complete linear-system without doing
		//any communication.
		//
		const int ghost = 1;
		if (local_elem_box[0][0] > 0)
			local_elem_box[0][0] -= ghost;
		if (local_elem_box[1][0] > 0)
			local_elem_box[1][0] -= ghost;
		if (local_elem_box[2][0] > 0)
			local_elem_box[2][0] -= ghost;
		if (local_elem_box[0][1] < mesh->global_box[0][1])
			local_elem_box[0][1] += ghost;
		if (local_elem_box[1][1] < mesh->global_box[1][1])
			local_elem_box[1][1] += ghost;
		if (local_elem_box[2][1] < mesh->global_box[2][1])
			local_elem_box[2][1] += ghost;

		perform_element_loop(*mesh, local_elem_box, *A, *b);
	}

	#pragma oss task						\
		inout(*A)						\
		inout(A_rows[0; A_nrows])				\
		inout(A_row_offsets[0; A_nrows + 1])			\
		inout(A_packed_cols[0; A_nnz])				\
		inout(A_packed_coefs[0; A_nnz])				\
		inout(*b)						\
		inout(b_coefs[0; b_local_size])				\
		inout(bc_rows_array[0: bc_rows_size])
	void impose_dirichlet(
		double prescribed_value,
		CSRMatrix *A,
		int *A_rows,
		int *A_row_offsets,
		size_t A_nrows,
		int *A_packed_cols,
		double *A_packed_coefs,
		size_t A_nnz,
		Vector *b,
		double *b_coefs,
		size_t b_local_size,
		int global_nx, int global_ny, int global_nz,
		const int *bc_rows_array,
		size_t bc_rows_size)
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

