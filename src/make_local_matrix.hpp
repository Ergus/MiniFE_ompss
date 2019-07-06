#ifndef _make_local_matrix_hpp_
#define _make_local_matrix_hpp_
#include <assert.h>

#include <fstream>
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

#include "utils.hpp"
#include "ompss_utils.hpp"
#include "CSRMatrix.hpp"
#include "singleton.hpp"

#include <vector>
#include <map>

namespace miniFE {

	#pragma oss task						\
		inout(A[0])						\
		in(nrows_array[0; numboxes])				\
		in(start_row_array[0; numboxes])			\
		in(stop_row_array[0; numboxes])				\
		in(rows[0; nrows])					\
		in(row_offsets[0; nrows + 1])				\
		in(packed_cols[0; nnz])					\
		out(recv_neighbors[0; numboxes])			\
		out(recv_length[0; numboxes])
	void get_recv_info_task(CSRMatrix *A, size_t id, size_t numboxes,
	                        const size_t *nrows_array,
	                        const int *start_row_array,
	                        const int *stop_row_array,
	                        int *rows,
	                        int *row_offsets,
	                        int *packed_cols,
	                        int *recv_neighbors,
	                        int *recv_length,
	                        size_t nnz, size_t nrows)
	{

		// First count and find the external elements
		// And invert the externals
		size_t num_external = 0;
		std::map<int,int> externals;
		std::vector<int> external_index_vector; // temporal.
		const int start_row = start_row_array[id];
		const int stop_row = stop_row_array[id];
		const int local_nrow = nrows_array[id];

		for (size_t i = 0; i < nrows; ++i) {
			int *Acols = NULL;
			double *Acoefs = NULL;
			size_t row_len = 0;
			A->get_row_pointers(A->rows[i], row_len, Acols, Acoefs);

			for (size_t j = 0; j < row_len; ++j) {
				const int cur_ind = Acols[j];

				if (start_row <= cur_ind &&
				    cur_ind <= stop_row) {
					Acols[j] -= start_row;
				} else { // Must find out if we have already set up this point
					if (externals.find(cur_ind) == externals.end()) {
						externals[cur_ind] = num_external++;
						external_index_vector.push_back(cur_ind);
					}
					// Mark index as external by adding 1 and negating it
					Acols[j] = -(Acols[j] + 1);
				}
			}
		}

		// Go through list of externals and find the processor that owns each
		std::vector<int> external_processor_vector(num_external, -1);

		for (size_t i = 0; i < num_external; ++i) {
			const int cur_ind = external_index_vector[i];
			for (int j = numboxes - 1; j >= 0; --j) {
				if (0 <= start_row_array[j] &&
				    start_row_array[j] <= cur_ind && cur_ind <= stop_row_array[j]) {
					external_processor_vector[i] = j;
					break;
				}
			}
			// test if a processor was found for this ind
			assert(external_processor_vector[i] >= 0);
		}

		// Filling the externals
		int *recv_list_local = (int *) rrl_malloc(numboxes * sizeof(int));
		int *recv_length_local = (int *) rrl_malloc(numboxes * sizeof(int));
		int *external_local_index_local = (int *) rrl_malloc(num_external * sizeof(int));
		int *external_grouped_index_local = (int *) rrl_malloc(num_external * sizeof(int));

		for (size_t i = 0; i < num_external; ++ i)
			external_local_index_local[i] = -1; // initialize to -1


		size_t count = 0, count_proc = 0;
		for (size_t i = 0; i < num_external; ++i) {
			if (external_local_index_local[i] < 0) {
				external_local_index_local[i] = count + nrows_array[id];
				external_grouped_index_local[count] = external_index_vector[i];
				++count;
				int count_i = 1;

				for(size_t j = i + 1; j < num_external; ++j) {
					if (external_processor_vector[j] ==
					    external_processor_vector[i]) {
						external_local_index_local[j]
							= local_nrow + count;
						external_grouped_index_local[count]
							= external_index_vector[j];

						++count;
						++count_i;
					}
				}
				recv_list_local[count_proc] = external_processor_vector[i];

				recv_length_local[count_proc] = count_i;
				++count_proc;
			}
		}

		assert(count_proc <= numboxes);

		A->nrecv_neighbors = count_proc;
		A->nexternals = num_external;
		ompss_memcpy_task(recv_neighbors, recv_list_local, count_proc * sizeof(int));
		ompss_memcpy_task(recv_length, recv_length_local, count_proc * sizeof(int));

		A->external_index = (int *) rrd_malloc(num_external * sizeof(int));
		ompss_memcpy_task(A->external_index, external_grouped_index_local,
		                  num_external * sizeof(int));

		//Change index of externals
		for (size_t i = 0; i < A->nrows; ++i) {
			int *Acols = NULL;
			double *Acoefs = NULL;
			size_t row_len = 0;
			A->get_row_pointers(A->rows[i], row_len, Acols, Acoefs);

			for (size_t j = 0; j < row_len; ++j) {
				if (Acols[j] < 0) { // Change index values of externals
					const int cur_ind = -Acols[j] - 1;
					Acols[j] = external_local_index_local[externals[cur_ind]];
				}
			}
		}

		#pragma oss taskwait
		// Release local memory
		// Filling the externals
		rrl_free(recv_list_local, numboxes * sizeof(int));
		rrl_free(recv_length_local, numboxes * sizeof(int));
		rrl_free(external_local_index_local, num_external * sizeof(int));
		rrl_free(external_grouped_index_local, num_external * sizeof(int));
	}

	// TODO try to substitute A with first private
	#pragma oss task inout(*A)					\
		in(nrecv_neighbors_global[0; numboxes])			\
		in(recv_neighbors_global[0; global_nrecv_neighbors])	\
		in(recv_length_global[0; global_nrecv_neighbors])	\
		out(send_neighbors_local[0; numboxes])			\
		out(send_length_local[0; numboxes])
	void get_send_info_task(CSRMatrix *A, size_t id,
		size_t numboxes,
		const int *nrecv_neighbors_global,

		int global_nrecv_neighbors,
		const int *recv_neighbors_global,
		const int *recv_length_global,

		int *send_neighbors_local,
		int *send_length_local)
	{
		size_t nsend_neighbors = 0;
		int nelements_to_send = 0;

		const int * it_recv_neighbors = recv_neighbors_global;
		const int * it_recv_length = recv_length_global;

		for (size_t i = 0; i < numboxes; ++i) {
			const int nrecv_neighbors_i = nrecv_neighbors_global[i];

			if (i != id) { // Not send to myself

				for (int j = 0; j < nrecv_neighbors_i; ++j) {
					if ((size_t)it_recv_neighbors[j] == id) { // I am the sender.
						const int send_i = it_recv_length[j];

						send_neighbors_local[nsend_neighbors] = i;
						send_length_local[nsend_neighbors] = send_i;

						nelements_to_send += send_i;
						nsend_neighbors++;

						assert(nsend_neighbors < numboxes);
						break;
					}
				}
			}
			it_recv_neighbors += nrecv_neighbors_i;
			it_recv_length += nrecv_neighbors_i;
		}

		A->nsend_neighbors = nsend_neighbors;
		A->nelements_to_send = nelements_to_send;
	}

	// TODO: possible error here
	#pragma oss task in(A_array[0: id - 1])				\
		inout(A_array[id])					\
		in(A_array[id + 1: numboxes - 1])			\
									\
		in(send_neighbors_local[0; nsend_neighbors_local])	\
		in(send_length_local[0; nsend_neighbors_local])		\
									\
		in(recv_neighbors_global[0; global_nrecv_neighbors])	\
		in(recv_length_global[0; global_nrecv_neighbors])	\
		weakout(recv_ptr_global[0; global_nrecv_neighbors])	\
									\
		in(external_index_global[0; global_nexternals_global])	\
									\
		out(elements_to_send_local[0; nelements_to_send_local])
	void set_send_info_task(CSRMatrix *A_array, size_t id,
	                        size_t numboxes,

	                        int nsend_neighbors_local,
	                        const int *send_neighbors_local,
	                        const int *send_length_local,

	                        int global_nrecv_neighbors,
	                        const int *recv_neighbors_global,
	                        const int *recv_length_global,
	                        double **recv_ptr_global,

	                        int global_nexternals_global,
	                        const int *external_index_global,

	                        int nelements_to_send_local,
	                        int *elements_to_send_local)
	{

		CSRMatrix *A = &A_array[id];
		// Remember this creates tasks internally

		// Fill the elements_to_send array
		const int nrows = A->nrows;
		const int start_row = nrows > 0 ? A->first_row : -1;
		const int stop_row = nrows > 0 ? A->rows[nrows - 1] : -1;

		int *ptr_local = A->elements_to_send; 	// Where I will put the elements for this neighbors
		double *ptr_for_remote = A->send_buffer;

		for (int i = 0; i < A->nsend_neighbors; ++i) { // Iterate over neighbors list
			const size_t send_neighbor_i = A->send_neighbors[i]; // neighbor

			assert(send_neighbor_i < numboxes);

			const CSRMatrix *A_i = &A_array[send_neighbor_i];

			// Where the neighbor has the ids I will send to it
			const int *ptr_remote = A_i->external_index;

			for (int j = 0; j < A_i->nrecv_neighbors; ++j) {
				if ((size_t)A_i->recv_neighbors[j] == id) {  // If I am the sender

					// How many elements I will send to it
					const int nsend_to_i = A_i->recv_length[j];

					#pragma oss task		\
						out(A_i->recv_ptr[j])	\
						in(ptr_remote[0; nsend_to_i]) \
						out(ptr_local[0; nsend_to_i])
					{
						// inform the remote about my pointer
						A_i->recv_ptr[j] = ptr_for_remote;

						for (int k = 0; k < nsend_to_i; ++k) {
							const int id_to_send_global = ptr_remote[k];

							// Assert I have this element
							assert(start_row <= id_to_send_global);
							assert(id_to_send_global <= stop_row);

							ptr_local[k] = id_to_send_global - start_row;
						}

					}
					// move local pointers
					ptr_local += nsend_to_i;
					ptr_for_remote += nsend_to_i;
				}
				// move remote pointer because it is not me
				ptr_remote += A_i->recv_length[j];
			}
		}

		// The elements_id is filled and moved to relative indices.
		A->num_cols = nrows + A->nexternals;
		A->has_local_indices = true;
	}


	// This will be called from driver ================================

	void make_local_matrix(CSRMatrix *A_array, singleton *sing, size_t numboxes)
	{
		if (numboxes < 2) {
			A_array[0].num_cols = A_array[0].nrows;
			A_array[0].has_local_indices = true;
			return;
		}

		// OmpSs arrays to access from tasks
		size_t *nrows_array = (size_t *) rrd_malloc(numboxes * sizeof(size_t));
		int *start_row_array = (int *) rrd_malloc(numboxes * sizeof(int));
		int *stop_row_array = (int *) rrd_malloc(numboxes * sizeof(int));

		// Array for reduction task
		int *recv_neighbors_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int)); // Process
		int *recv_length_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int)); // Process

		int *send_neighbors_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int)); // Process
		int *send_length_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int));

		// Boundary information
		for (size_t id = 0; id < numboxes; ++id) {
			#pragma oss task				\
				in(A_array[id])				\
				in(A_array[id].rows[A_array[id].nrows - 1]) \
				out(start_row_array[id])		\
				out(stop_row_array[id])			\
				out(nrows_array[id])
			{
				const size_t local_nrow = A_array[id].nrows;
				nrows_array[id] = local_nrow;

				if (local_nrow > 0) {
					start_row_array[id] = A_array[id].first_row;
					stop_row_array[id] = A_array[id].rows[local_nrow - 1];
				} else {
					start_row_array[id] = -1;
					stop_row_array[id] = -1;
				}
			}
		}

		// Find the external elements (recv information).
		// Scan the indices and transform to local
		for (size_t id = 0; id < numboxes; ++id) {
			//int *tmp_neighbors_local = &tmp_neighbors_global[numboxes * id];
			CSRMatrix *A = &A_array[id];
			A->recv_neighbors = &recv_neighbors_global[id * numboxes];
			A->recv_length = &recv_length_global[id * numboxes];

			// This is already a task (look before)
			get_recv_info_task(A, id, numboxes,
			                   nrows_array,
			                   start_row_array,
			                   stop_row_array,
			                   A->rows,
			                   A->row_offsets,
			                   A->packed_cols,
			                   A->recv_neighbors,
			                   A->recv_length,
			                   A->nnz, A->nrows);
		}

		#pragma oss taskwait

		{// This allocates the recv information in a single huge array
			int global_nrecv_neighbors = 0;
			int *nrecv_neighbors_offset = (int *) alloca(numboxes * sizeof(int));

			int nexternals = 0;
			int *nexternals_offset = (int *) alloca(numboxes * sizeof(int));;

			for (size_t id = 0; id < numboxes; ++id) {
				sing->nrecv_neighbors[id] = A_array[id].nrecv_neighbors;
				sing->nexternals[id] = A_array[id].nexternals;

				nrecv_neighbors_offset[id] = global_nrecv_neighbors;
				global_nrecv_neighbors += A_array[id].nrecv_neighbors;

				nexternals_offset[id] = nexternals;
				nexternals += A_array[id].nexternals;
			}

			// this allocates the arrays internally with dmalloc
			sing->allocate_recv(global_nrecv_neighbors, nexternals);

			// The next extra copy enables a latter optimization, it
			// is in parallel
			for (size_t id = 0; id < numboxes; ++id) {

				ompss_memcpy_task(&(sing->recv_neighbors[nrecv_neighbors_offset[id]]),
				                  A_array[id].recv_neighbors,
				                  sing->nrecv_neighbors[id] * sizeof(int));

				ompss_memcpy_task(&(sing->recv_length[nrecv_neighbors_offset[id]]),
				                  A_array[id].recv_length,
				                  sing->nrecv_neighbors[id] * sizeof(int));

				ompss_memcpy_task(&(sing->external_index[nexternals_offset[id]]),
				                  A_array[id].external_index,
				                  sing->nexternals[id] * sizeof(int));

			}

			#pragma oss taskwait

			for (size_t id = 0; id < numboxes; ++id) {
				rrd_free(A_array[id].external_index, sing->nexternals[id] * sizeof(int));

				A_array[id].recv_neighbors =
					&(sing->recv_neighbors[nrecv_neighbors_offset[id]]);
				A_array[id].recv_length =
					&(sing->recv_length[nrecv_neighbors_offset[id]]);
				A_array[id].recv_ptr =
					&(sing->recv_ptr[nrecv_neighbors_offset[id]]);
				A_array[id].external_index =
					&(sing->external_index[nexternals_offset[id]]);
			}

		}

		// Fill send Information
		for (size_t id = 0; id < numboxes; ++id) {
			CSRMatrix *A = &A_array[id];

			A->send_neighbors = &send_neighbors_global[id * numboxes];
			A->send_length = &send_length_global[id * numboxes];

			get_send_info_task(A, id,
			                   numboxes,
			                   sing->nrecv_neighbors,
			                   sing->global_nrecv_neighbors,
			                   sing->recv_neighbors,
			                   sing->recv_length,
			                   A->send_neighbors,
			                   A->send_length);
		}

		#pragma oss taskwait

		{// This allocates the recv information in a single huge array
			int global_nsend_neighbors = 0;
			int *nsend_neighbors_offset = (int *) alloca(numboxes * sizeof(int));

			int global_nelements_to_send = 0;
			int *nelements_to_send_offset = (int *) alloca(numboxes * sizeof(int));

			for (size_t id = 0; id < numboxes; ++id) {
				sing->nsend_neighbors[id] = A_array[id].nsend_neighbors;
				sing->nelements_to_send[id] = A_array[id].nelements_to_send;

				nsend_neighbors_offset[id] = global_nsend_neighbors;
				global_nsend_neighbors += A_array[id].nsend_neighbors;

				nelements_to_send_offset[id] = global_nelements_to_send;
				global_nelements_to_send += A_array[id].nelements_to_send;

			}

			// this allocated the arrays internally
			sing->allocate_send(global_nsend_neighbors, global_nelements_to_send);

			for (size_t id = 0; id < numboxes; ++id) {
				ompss_memcpy_task(&(sing->send_neighbors[nsend_neighbors_offset[id]]),
				                  A_array[id].send_neighbors,
				                  A_array[id].nsend_neighbors * sizeof(int));

				ompss_memcpy_task(&(sing->send_length[nsend_neighbors_offset[id]]),
				                  A_array[id].send_length,
				                  A_array[id].nsend_neighbors * sizeof(int));
			}

			#pragma oss taskwait

			for (size_t id = 0; id < numboxes; ++id) {
				A_array[id].send_neighbors =
					&(sing->send_neighbors[nsend_neighbors_offset[id]]);
				A_array[id].send_length =
					&(sing->send_length[nsend_neighbors_offset[id]]);
				A_array[id].elements_to_send =
					&(sing->elements_to_send[nelements_to_send_offset[id]]);
				A_array[id].send_buffer =
					&(sing->send_buffer[nelements_to_send_offset[id]]);
			}

		}

		for (size_t id = 0; id < numboxes; ++id) {

			CSRMatrix *A = &A_array[id];

			set_send_info_task(A_array, id,
			                   numboxes,

			                   A->nsend_neighbors,
			                   A->send_neighbors,
			                   A->send_length,

			                   sing->global_nrecv_neighbors,
			                   sing->recv_neighbors,
			                   sing->recv_length,
			                   sing->recv_ptr,

			                   sing->global_nexternals,
			                   sing->external_index,

			                   A->nelements_to_send,
			                   A->elements_to_send);

		}

		#pragma oss taskwait

		rrd_free(nrows_array, numboxes * sizeof(size_t));
		rrd_free(start_row_array, numboxes * sizeof(int));
		rrd_free(stop_row_array, numboxes * sizeof(int));

		rrd_free(recv_neighbors_global, numboxes * numboxes * sizeof(int));
		rrd_free(recv_length_global, numboxes * numboxes * sizeof(int));

		rrd_free(send_neighbors_global, numboxes * numboxes * sizeof(int));
		rrd_free(send_length_global, numboxes * numboxes * sizeof(int));


	}

}//namespace miniFE

#endif
