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

namespace miniFE
{

	class CSRMatrix {
	public:
		int id;
		bool has_local_indices, is_copy;
		int *rows;
		int *row_offsets;
		int *row_coords;
		//int *row_offsets_external;
		size_t global_nrows, nrows, nnz;

		int *packed_cols;
		double *packed_coefs;

		int first_row, stop_row;
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
			  rows(nullptr), row_offsets(nullptr), row_coords(nullptr),
			  global_nrows(-1), nrows(-1),
			  nnz(0), packed_cols(nullptr), packed_coefs(nullptr),
			  first_row(-1), stop_row(-1), num_cols(0),

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
			  stop_row(in.stop_row),
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
			cols = nullptr;
			coefs = nullptr;
			row_length = 0;

			// This solves a race condition in the original code
			// and optimizes a very frequent case.
			if (row < first_row || row > stop_row || nrows == 0)
				return;

			ptrdiff_t local_row = -1;

			//first see if we can get the local-row index using fast direct lookup:
			assert(row >= first_row);

			const size_t idx = row - first_row;
			if (idx < nrows && rows[idx] == row) {
				local_row = idx;
			} else {
				//if we didn't get the local-row index using direct lookup, try a
				//more expensive binary-search:

				auto row_iter = std::lower_bound(rows, &rows[nrows], row);

				// I keep  this because it is in the original code.
				// But with the first check this should never happen.
				//if we still haven't found row, it's not local so jump out:
				if (row_iter == &rows[nrows] || *row_iter != row) {
					row_length = 0;
					return;
				}

				local_row = row_iter - rows;
			}

			assert(local_row < nrows);

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

		void write(std::ostream &stream = std::cout) const
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

			dbvprint_vector("\nrecv_neighbors", nrecv_neighbors, recv_neighbors, stream);
			dbvprint_vector("\nrecv_length", nrecv_neighbors, recv_length, stream);
			dbvprint_vector("\nexternals", nexternals, external_index, stream);

			dbvprint_vector("\nsend_neighbors", nsend_neighbors, send_neighbors, stream);
			dbvprint_vector("\nsend_length", nsend_neighbors, send_length, stream);
			dbvprint_vector("\nnelements_to_send", nelements_to_send, elements_to_send, stream);

			#endif

		}

	}; // class CSRMatrix

	void generate_matrix_structure_all(CSRMatrix *A_array,
	                                   const simple_mesh_description *mesh_array,
	                                   singleton *sing,
	                                   size_t numboxes)
	{
		for (int i = 0; i < numboxes; ++i) {

			CSRMatrix *A = &A_array[i];
			const simple_mesh_description *mesh = &mesh_array[i];

			dbvwrite(mesh);


			A->id = i;
			A->nrows = mesh->extended_box.get_num_ids();

			//num-owned-nodes in each dimension is num-elems+1
			//only if num-elems > 0 in that dimension *and*
			//we are at the high end of the global range in that dimension:
			//A->rows = (int *) rrd_malloc(A->nrows * sizeof(int));
			//A->row_offsets = (int *) rrd_malloc((A->nrows + 1) * sizeof(int));
			//A->row_coords = (int *) rrd_malloc(A->nrows * 3 * sizeof(int));

		}

		sing->allocate_rows(A_array);

		for (int i = 0; i < numboxes; ++i) {

			CSRMatrix *A = &A_array[i];
			const simple_mesh_description *mesh = &mesh_array[i];
			int global_nodes[3] = {
				mesh->global_box[0][1] + 1,
				mesh->global_box[1][1] + 1,
				mesh->global_box[2][1] + 1 };

			A->global_nrows = global_nodes[0] * global_nodes[1] * global_nodes[2];

			init_offsets_task(A->row_coords,
			                  A->rows,
			                  A->row_offsets,
			                  global_nodes,
			                  mesh,
			                  mesh->ompss2_ids_to_rows,
			                  mesh->ids_to_rows_size,
			                  A->global_nrows,
			                  A->nrows,
			                  &(A->nnz),
			                  &(A->first_row),
			                  &(A->stop_row));
		}

		#pragma oss taskwait

		for (int i = 0; i < numboxes; ++i) {
			CSRMatrix *A = &A_array[i];
			const simple_mesh_description *mesh = &mesh_array[i];

			A->packed_cols = (int *) rrd_malloc(A->nnz * sizeof(int));
			A->packed_coefs = (double *) rrd_malloc(A->nnz * sizeof(double));

			int global_nodes[3] = {
				mesh->global_box[0][1] + 1,
				mesh->global_box[1][1] + 1,
				mesh->global_box[2][1] + 1 };

			init_matrix_task(A->row_coords,
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
		}
		#pragma oss taskwait
	}

	void matvec_task(CSRMatrix *A, Vector *x, Vector *y,
	                 size_t id = 0, bool print  = false)
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
			in(A[0])				\
			in(Arow_offsets[0; Anrows + 1])		\
			in(Apacked_cols[0; Annz])		\
			in(Apacked_coefs[0; Annz])		\
			in(x[0])				\
			in(xcoefs[0; xlocal_size])		\
			in(y[0])				\
			out(ycoefs[0; ylocal_size])
		{
			if (print) {
				std::string filename = "VERB_matvec_" + std::to_string(id) + ".verb";
				std::ofstream stream(filename);

				A->write(stream);
				x->write(stream);
				stream.close();
			}

			const double beta = 0.0;  // I don't really understand what is this for

			for (size_t row = 0; row < Anrows; ++row) {
				double sum = beta * ycoefs[row];

				for(int i = Arow_offsets[row]; i < Arow_offsets[row + 1]; ++i) {
					const int col = Apacked_cols[i];
					assert((size_t)i < Annz);
					sum += Apacked_coefs[i] * xcoefs[col];
				}

				//std::cout << "row[" << row << "] = " << sum << std::endl;
				ycoefs[row] = sum;
			}
		}
	}

	inline void assemble_FE_data_task(size_t id,
		const simple_mesh_description *mesh,
		CSRMatrix *A,
		Vector *b)
	{

		const std::pair<int,int> *mesh_ompss2_ids_to_rows = mesh->ompss2_ids_to_rows;
		size_t mesh_i_ids_to_rows_size = mesh->ids_to_rows_size;
		int *Arows = A->rows;
		int *Arow_offsets = A->row_offsets;
		size_t Anrows = A->nrows;
		int *Apacked_cols = A->packed_cols;
		double *Apacked_coefs = A->packed_coefs;
		size_t Annz = A->nnz;
		double *bcoefs = b->coefs;
		size_t blocal_size = b->local_size;

		#pragma oss task					\
			in(mesh[0])					\
			in(mesh_ompss2_ids_to_rows[0; mesh_i_ids_to_rows_size])	\
			inout(A[0])					\
			inout(Arows[0; Anrows])				\
			inout(Arow_offsets[0; Anrows + 1])		\
			inout(Apacked_cols[0; Annz])			\
			inout(Apacked_coefs[0; Annz])			\
			inout(b[0])					\
			inout(bcoefs[0; blocal_size])
		{
			dbvwrite(mesh);
			dbvwrite(A);
			dbvwrite(b);

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

			dbvwrite(mesh);
			dbvwrite(A);
			dbvwrite(b);
		}
	}


	void impose_dirichlet_task(
		double prescribed_value,
		const CSRMatrix *A,
		const Vector *b,
		int global_nx, int global_ny, int global_nz,
		const int *bc_rows_array,
		size_t bc_rows_size)
	{
		int *Arows = A->rows;
		int *Arow_offsets = A->row_offsets;
		size_t Anrows = A->nrows;
		int *Apacked_cols = A->packed_cols;
		double *Apacked_coefs = A->packed_coefs;
		size_t Annz = A->nnz;
		double *bcoefs = b->coefs;
		size_t blocal_size = b->local_size;

		#pragma oss task					\
			in(A[0])					\
			in(Arows[0; Anrows])				\
			in(Arow_offsets[0; Anrows + 1])		\
			in(Apacked_cols[0; Annz])			\
			inout(Apacked_coefs[0; Annz])			\
			in(b[0])					\
			inout(bcoefs[0; blocal_size])			\
			in(bc_rows_array[0: bc_rows_size])
		{
			const int first_local_row = Anrows > 0 ? Arows[0] : 0;
			const int last_local_row  = Anrows > 0 ? Arows[Anrows - 1] : -1;

			for (size_t i = 0; i < bc_rows_size; ++i) {
				const int row = bc_rows_array[i];

				if (row >= first_local_row && row <= last_local_row) {
					const size_t local_row = row - first_local_row;
					bcoefs[local_row] = prescribed_value;
					zero_row_and_put_1_on_diagonal(A, row);
				}
			}

			for (size_t i = 0; i < Anrows; ++i) {
				const int row = Arows[i];

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

				bcoefs[i] -= sum * prescribed_value;
			}
		}
	}






}//namespace miniFE

#endif

