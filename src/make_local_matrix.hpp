#ifndef _make_local_matrix_hpp_
#define _make_local_matrix_hpp_
#include <assert.h>

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

#include <vector>
#include <map>

namespace miniFE {

	#pragma oss task						\
		inout(*A)						\
		in(rows[0; nrows])				\
		in(row_offsets[0; nrows + 1])		\
		in(packed_cols[0; nnz])			\
		inout(packed_coefs[0; nnz])			\
		in(nrows_array[0; numboxes])				\
		in(start_row_array[0; numboxes])			\
		in(stop_row_array[0; numboxes])				\
		out(recv_list_local[0; numboxes])			\
		out(nrecv_list_local)					\
		out(recv_length_local[0; numboxes])			\
		out(nrecv_length_local)
	void fill_recv_task(CSRMatrix *A,
	                    size_t id,
	                    size_t numboxes,
	                    const size_t *nrows_array,
	                    const int *start_row_array,
	                    const int *stop_row_array,
	                    int *recv_list_local,
	                    int &nrecv_list_local,
	                    int *recv_length_local,
	                    int &nrecv_length_local,
	                    int *rows,
	                    int *row_offsets,
	                    int *packed_cols,
	                    double *packed_coefs,
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
		int *external_local_index_local = (int *) rrl_malloc(num_external * sizeof(int));
		for (size_t i = 0; i < num_external; ++ i)
			external_local_index_local[i] = -1; // initialize to -1

		int *external_grouped_index_local = (int *) rrl_malloc(num_external * sizeof(int));

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

		// Copies to simplify tasks
		nrecv_list_local = count_proc;
		nrecv_length_local = num_external;

		// To write in A
		A->nrecv_neighbors = count_proc;
		A->nexternals = num_external;
		A->recv_neighbors = (int *) rrd_malloc(count_proc * sizeof(int));
		A->recv_ptr = (double **) rrd_malloc(count_proc * sizeof(double *));
		A->recv_length = (int *) rrd_malloc(count_proc * sizeof(int));
		A->external_index = (int *) rrd_malloc(num_external * sizeof(int));

		// TODO: reassign pointers instead of copy
		// This copy creates tasks internally
		copy_local_to_global_task(A->recv_neighbors,
		                          recv_list_local, count_proc);
		copy_local_to_global_task(A->recv_length,
		                          recv_length_local, count_proc);
		copy_local_to_global_task(A->external_index,
		                          external_grouped_index_local, num_external);

		//Change index of externals
		int count_assert = 0;
		for (size_t i = 0; i < A->nrows; ++i) {
			int *Acols = NULL;
			double *Acoefs = NULL;
			size_t row_len = 0;
			A->get_row_pointers(A->rows[i], row_len, Acols, Acoefs);

			for (size_t j = 0; j < row_len; ++j) {
				if (Acols[j] < 0) { // Change index values of externals
					const int cur_ind = -Acols[j] - 1;
					Acols[j] = external_local_index_local[externals[cur_ind]];
					++count_assert;
				}
			}
		}

		#pragma oss taskwait
		// Release local memory
		rrl_free(external_local_index_local, num_external * sizeof(int));
		rrl_free(external_grouped_index_local, num_external * sizeof(int));
	}

	#pragma oss task						\
		inout(*A)						\
		in(recv_list_global[0; numboxes * numboxes])		\
		in(nrecv_list_global[0; numboxes])			\
		in(recv_length_global[0; numboxes * numboxes])		\
		out(send_list_local[0; numboxes])			\
		out(nsend_list_local)					\
		out(send_length_local[0; numboxes])			\
		out(nsend_length_local)
	void get_send_info_task(CSRMatrix *A,
	                        size_t id,
	                        size_t numboxes,
	                        const int *recv_list_global,
	                        const int *nrecv_list_global,
	                        const int *recv_length_global,
	                        int *send_list_local,
	                        int &nsend_list_local,
	                        int *send_length_local,
	                        int &nsend_length_local)
	{
		///////////////////////////////////////////////////////////////////////
		///
		// Make a list of the neighbors that will send information to update our
		// external elements (in the order that we will receive this information).
		///
		///////////////////////////////////////////////////////////////////////

		// Construct send_list (substitutes send-recv code)
		// Send a 0 length message to each of our recv neighbors

		nsend_list_local = 0;
		nsend_length_local = 0;
		for (size_t i = 0; i < numboxes; ++i) {
			if (i != id) { // Not send to myself
				const int nrecv_list_remote_i = nrecv_list_global[i];
				const int *recv_list_i = &recv_list_global[i * numboxes];
				const int *recv_length_i = &recv_length_global[i * numboxes];

				for (int j = 0; j < nrecv_list_remote_i; ++j) {
					if ((size_t)recv_list_i[j] == id) {
						const int send_i = recv_length_i[j];
						send_list_local[nsend_list_local] = i;
						send_length_local[nsend_list_local] = send_i;
						nsend_length_local += send_i;
						++nsend_list_local;
						assert(nsend_list_local < (int)numboxes);
						break;
					}
				}
			}
		}

		A->nsend_neighbors = nsend_list_local;
		A->nelements_to_send = nsend_length_local;

		// Just allocate write in another task
		A->send_neighbors = (int *) rrd_malloc(nsend_list_local * sizeof(int));
		A->send_length = (int *) rrd_malloc(nsend_list_local * sizeof(int));
		A->elements_to_send = (int *) rrd_malloc(nsend_length_local * sizeof(int));
		A->send_buffer = (double *) rrd_malloc(nsend_length_local * sizeof(double));
	}

	// TODO: possible error here
	#pragma oss task						\
		in(A_array[0; numboxes])				\
		in(nrows_array[id])					\
		in(start_row_array[id])					\
		in(stop_row_array[id])					\
		in(recv_list_global[0; numboxes * numboxes])		\
		in(nrecv_list_global[0; numboxes])			\
		in(recv_length_global[0; numboxes * numboxes])		\
		in(nrecv_length_global[0; numboxes])			\
		in(send_list_local[0; numboxes])			\
		in(send_length_local[0; numboxes])			\
		out(Asend_neighbors[0; nsend_list_local])		\
		out(Asend_length[0; nsend_list_local])			\
		out(Aelements_to_send[0; nsend_length_local])
	void set_send_info_task(CSRMatrix *A_array,
	                        size_t id,
	                        size_t numboxes,
	                        const size_t *nrows_array,
	                        const int *start_row_array,
	                        const int *stop_row_array,
	                        const int *recv_list_global,
	                        const int *nrecv_list_global,
	                        const int *recv_length_global,
	                        const int *nrecv_length_global,
	                        const int *send_list_local,
	                        const int nsend_list_local,
	                        const int *send_length_local,
	                        const int nsend_length_local,
	                        int *Asend_neighbors,
	                        int *Asend_length,
	                        int *Aelements_to_send)
	{

		CSRMatrix *A = &A_array[id];
		// Remember this creates tasks internally
		memcpy(A->send_neighbors, send_list_local, A->nsend_neighbors);
		memcpy(A->send_length, send_length_local, A->nsend_neighbors);

		// Fill the elements_to_send array
		const int start_row = start_row_array[id];
		const int stop_row = stop_row_array[id];
		int start_local = 0;
		for (int i = 0; i < A->nsend_neighbors; ++i) { // Iterate over neighbors list
			const size_t send_neighbor_id = send_list_local[i]; // neighbor

			assert(send_neighbor_id < numboxes);

			// number of neighbors that will send to neighbors
			const int nrecv_list_id =
				nrecv_list_global[send_neighbor_id];

			// list of neighbors of the send_neighbors size (nrecv_list_id)
			const int *recv_list_id =
				&recv_list_global[send_neighbor_id * numboxes];

			// remote lengths
			const int *recv_length_id =
				&recv_length_global[send_neighbor_id * numboxes];

			// Where the neighbor has the ids I will send to it
			int *ptr_remote = A_array[send_neighbor_id].external_index;

			for (int j = 0; j < nrecv_list_id; ++j) {
				if ((size_t)recv_list_id[j] == id) {  // If I am the sender

					// Where I will put the elements for this neighbors
					int *ptr_local = &(A->elements_to_send[start_local]);

					// How many elements I will send to it
					const int nsend_to_id = recv_length_id[j];

					// TODO: errors come from  HERE.
					#pragma oss task		\
						in(A)			\
						in(ptr_remote[0; nsend_to_id]) \
						out(A_array[send_neighbor_id].recv_ptr[j]) \
						out(ptr_local[0; nsend_to_id])
					{
						// inform the remote about my pointer
						A_array[send_neighbor_id].recv_ptr[j] =
							&A->send_buffer[start_local];

						for (int k = 0; k < nsend_to_id; ++k) {
							const int id_to_send_global = ptr_remote[k];

							// Assert I have this element
							assert(start_row <= id_to_send_global);
							assert(id_to_send_global <= stop_row);

							ptr_local[k] = id_to_send_global - start_row;
						}

					}

					start_local += nsend_to_id;
				}
				ptr_remote += recv_length_id[j];
			}
		}

		// The elements_id is filled and moved to relative indices.
		A->num_cols = nrows_array[id] + nrecv_length_global[id];
		A->has_local_indices = true;
	}


	// This will be called from driver ================================

	void make_local_matrix(CSRMatrix *A_array, size_t numboxes)
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
		// int *tmp_neighbors_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int));

		int *nrecv_list_global = (int *) rrd_malloc(numboxes * sizeof(int));
		int *recv_list_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int));
		int *recv_length_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int));
		int *nrecv_length_global = (int *) rrd_malloc(numboxes * sizeof(int));

		int *nsend_list_global = (int *) rrl_malloc(numboxes * sizeof(int));
		int *send_list_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int)); // Process
		int *send_length_global = (int *) rrd_malloc(numboxes * numboxes * sizeof(int));
		int *nsend_length_global = (int *) rrl_malloc(numboxes * sizeof(int));

		// Boundary information
		for (size_t id = 0; id < numboxes; ++id) {
			#pragma oss task		\
				in(A_array[id])		\
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
			int *recv_list_local = &recv_list_global[id * numboxes];
			int *recv_length_local = &recv_length_global[id * numboxes];
			CSRMatrix *A = &A_array[id];

			// This is already a task (look before)
			fill_recv_task(A, id, numboxes,
			               nrows_array,
			               start_row_array,
			               stop_row_array,
			               recv_list_local,
			               nrecv_list_global[id],
			               recv_length_local,
			               nrecv_length_global[id],
			               A->rows, A->row_offsets,
			               A->packed_cols, A->packed_coefs,
			               A->nnz,
			               A->nrows);


		}

		// Fill send Information
		for (size_t id = 0; id < numboxes; ++id) {

			int *send_list_local = &send_list_global[id * numboxes];
			int *send_length_local = &send_length_global[id * numboxes];


			get_send_info_task(&A_array[id],
			                   id,
			                   numboxes,
			                   recv_list_global,
			                   nrecv_list_global,
			                   recv_length_global,
			                   send_list_local,
			                   nsend_list_global[id],
			                   send_length_local,
			                   nsend_length_global[id]);

		}

		#pragma oss taskwait

		for (size_t id = 0; id < numboxes; ++id) {

			CSRMatrix *A = &A_array[id];
			int *send_list_local = &send_list_global[id * numboxes];
			int *send_length_local = &send_length_global[id * numboxes];

			set_send_info_task(A_array,
			                   id,
			                   numboxes,
			                   nrows_array,
			                   start_row_array,
			                   stop_row_array,
			                   recv_list_global,
			                   nrecv_list_global,
			                   recv_length_global,
			                   nrecv_length_global,
			                   send_list_local,
			                   nsend_list_global[id],
			                   send_length_local,
			                   nsend_list_global[id],
			                   A->send_neighbors,
			                   A->send_length,
			                   A->elements_to_send);
		}

		for (size_t id = 0; id < numboxes; ++id) {
			std::cout  << "id: " << id << std::endl;
			std::cout << A_array[id] << std::endl;
		}

		#pragma oss taskwait

		rrd_free(nrows_array, numboxes * sizeof(size_t));
		rrd_free(start_row_array, numboxes * sizeof(int));
		rrd_free(stop_row_array, numboxes * sizeof(int));

		rrd_free(recv_list_global, numboxes * numboxes * sizeof(int));
		rrd_free(nrecv_list_global, numboxes * sizeof(int));

		rrd_free(recv_length_global, numboxes * numboxes * sizeof(int));
		rrd_free(nrecv_length_global, numboxes * sizeof(int));

		rrd_free(send_list_global, numboxes * numboxes * sizeof(int));
		rrl_free(nsend_list_global, numboxes * sizeof(int));

		rrd_free(send_length_global, numboxes * numboxes * sizeof(int));
		rrl_free(nsend_length_global, numboxes * sizeof(int));

	}

}//namespace miniFE

#endif
