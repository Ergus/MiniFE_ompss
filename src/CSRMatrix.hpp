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
			  send_neighbors(send_neighbors),
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
			                  nnz,
			                  first_row);

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

	};


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

	template <typename T>
	void print_vector(std::string vname, size_t size,const  T *vect ,std::ostream &stream)
	{
		stream << vname << "["  << size << "]={";
		for (size_t i = 0; i < size; ++i) {
			if (i > 0)
				stream << "; ";

			stream << vect[i];
		}
		stream << " }\n";
	}

	inline std::ostream& operator <<(std::ostream &stream, const CSRMatrix &Min)
	{
		CSRMatrix Mcopy(Min);

		#pragma oss task					\
			in(Min)						\
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
			stream << "has_local_indices" << "="<< Min.has_local_indices << "\n";
			stream << "global_nrows" << "="<< Min.global_nrows << "\n";
			stream << "nrows" << "="<< Min.nrows << "\n";
			stream << "nnz" << "="<< Min.nnz << "\n";
			stream << "first_row" << "="<< Min.first_row << "\n";
			stream << "num_cols" << "="<< Min.num_cols << "\n";
			stream << "nrecv_neighbors" << "="<< Min.nrecv_neighbors << "\n";
			stream << "nexternals" << "="<< Min.nexternals << "\n";
			stream << "nsend_neighbors" << "="<< Min.nsend_neighbors << "\n";
			stream << "nelements_to_send" << "="<< Min.nelements_to_send << "\n";

			for (size_t i = 0; i < Min.nrows; ++i) {
				int *Acols = NULL;
				double *Acoefs = NULL;
				size_t row_len = 0;
				Min.get_row_pointers(Min.rows[i], row_len, Acols, Acoefs);

				stream << i << ":" << Min.rows[i] << ":" << Min.row_offsets[i] << "={";

				for (size_t j = 0; j < row_len; ++j) {
					if (j > 0)
						stream << "; ";
					stream << "<" << Acols[j] << ";" << Acoefs[j] << ">";
				}
				stream << "}\n";
			}

			print_vector("recv_neighbors", Min.nrecv_neighbors, Min.recv_neighbors, stream);
			print_vector("recv_length", Min.nrecv_neighbors, Min.recv_length, stream);
			print_vector("externals", Min.nexternals, Min.external_index, stream);

			print_vector("send_neighbors", Min.nsend_neighbors, Min.send_neighbors, stream);
			print_vector("send_length", Min.nsend_neighbors, Min.send_length, stream);
			print_vector("nelements_to_send", Min.nelements_to_send, Min.elements_to_send, stream);

			return stream;
		}
	}

	// TODO: pretty sure this is an out in y->coefs
	inline void write_all(std::string &filename, const CSRMatrix *A_array, size_t numboxes)
	{
		for (size_t id = 0; id < numboxes; ++id) {

			std::ofstream myfile;
			if (id == 0)
				myfile.open(filename, std::ofstream::out);
			else
				myfile.open(filename, std::ofstream::app);

			myfile << A_array[id] << std::endl;
			myfile.close();
		}
	}


}//namespace miniFE

#endif

