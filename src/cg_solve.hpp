#ifndef _cg_solve_hpp_
#define _cg_solve_hpp_

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

#include <cmath>
#include <limits>

#include "Vector_functions.hpp"
#include "mytimer.hpp"

#include "CSRMatrix.hpp"
#include "outstream.hpp"

#include "singleton.hpp"

namespace miniFE {

	#pragma oss task					   \
		in(A_elements_to_send[0; A_nelements_to_send])	   \
		out(A_send_buffer[0; A_nelements_to_send])	   \
		in(x_coefs[0; x_local_size])
	void copy_to_send_task(
		int *A_elements_to_send,
		double *A_send_buffer,
		size_t A_nelements_to_send,
		const double *x_coefs,
		size_t x_local_size)
	{
		for (size_t i = 0; i < A_nelements_to_send; ++i) {
			assert((size_t)A_elements_to_send[i] < x_local_size);
			A_send_buffer[i] = x_coefs[A_elements_to_send[i]];
		}
	}

	#pragma oss task						\
		in(Arecv_ptr[0; Anrecv_neighbors])			\
		in(Arecv_length[0; Anrecv_neighbors])			\
		weakout(x_start[0; Anexternals])				\
		weakin(global_send_buffer[0; global_nelements_to_send])
	void copy_from_senders(
		CSRMatrix *A,
		double **Arecv_ptr,
		int *Arecv_length,
		int Anrecv_neighbors,
		double *x_start,
		int Anexternals,
		double *global_send_buffer,
		int global_nelements_to_send)
	{
		double *local_iter = x_start;

		for (int i = 0; i < Anrecv_neighbors; ++i) {

			// This creates task internally
			#ifndef NDEBUG
			if (!Arecv_ptr[i]) {
				dbprintf("Error in vector A.recv_ptr[%d]\n", i);
				dbvprint_vector("Arecv_ptr", Anrecv_neighbors, Arecv_ptr, std::cerr);
				std::cerr<<std::endl;
			}
			assert(Arecv_ptr[i] != NULL);
			#endif

			ompss_memcpy_task(local_iter, Arecv_ptr[i], Arecv_length[i] * sizeof(double));

			local_iter += Arecv_length[i];

			// boundary check
			assert(local_iter - x_start <= Anexternals);
		}
	}


	void exchange_externals_all(CSRMatrix *A_array, Vector *x_array, size_t numboxes,
	                            const singleton *sing)
	{
		dbvprintf("Executing exchange externals\n");

		if (numboxes < 2)
			return;

		for (size_t id = 0; id < numboxes; ++id) {
			CSRMatrix *A = &A_array[id];
			Vector *x = &x_array[id];

			copy_to_send_task(
				A->elements_to_send,
				A->send_buffer,
				A->nelements_to_send,
				x->coefs,
				x->local_size);
		}

		for (size_t id = 0; id < numboxes; ++id) {

			CSRMatrix *A = &A_array[id];
			Vector *x = &x_array[id];

			copy_from_senders(A,
			                  A->recv_ptr,
			                  A->recv_length,
			                  A->nrecv_neighbors,
			                  &(x->coefs[A->nrows]),
			                  A->nexternals,
			                  sing->send_buffer,
			                  sing->global_nelements_to_send);
		}
	}

	inline int breakdown(double inner, double vnorm2, double wnorm2)
	{
		const double vnorm = std::sqrt(vnorm);
		const double wnorm = std::sqrt(wnorm);

		return std::abs(inner) <= 100 * vnorm * wnorm * std::numeric_limits<double>::epsilon();
	}


	void cg_solve_all(CSRMatrix *A_array, size_t numboxes,
	                  const Vector *b_array, Vector *x_array,
	                  int max_iter,
	                  const double tolerance,
	                  int &num_iters,
	                  double &normr,
	                  singleton *sing)
	{

		Vector *r_array = new Vector[numboxes];
		Vector *p_array = new Vector[numboxes];
		Vector *Ap_array = new Vector[numboxes];

		{
			init_vector_all(r_array, A_array, sing, numboxes);
			init_vector_all(Ap_array, A_array, sing, numboxes);

			for (size_t i = 0; i < numboxes; ++i)
				p_array[i].init(i, 0, A_array[i].num_cols);

			sing->allocate_vectors(p_array);
		}

		normr = 0;
		double rtrans[numboxes];
		double oldrtrans;
		double rtrans_global = 0.0;
		double beta = 0;

		int print_freq = max_iter / 10;
		if (print_freq > 50)
			print_freq = 50;
		if (print_freq < 1)
			print_freq = 1;

		for (size_t i = 0; i < numboxes; ++i) {
			rtrans[i] = 0;
			assert(A_array[i].has_local_indices);

			waxpby_task(1.0, &x_array[i], 0.0, &x_array[i], &p_array[i]);
		}

		exchange_externals_all(A_array, p_array, numboxes, sing);

		for (size_t i = 0; i < numboxes; ++i) {
			matvec_task(&A_array[i], &p_array[i], &Ap_array[i]);

			waxpby_dot_task(1.0, &b_array[i], -1.0, &Ap_array[i], &r_array[i], &rtrans[i]);
		}

		reduce_sum_task(&rtrans_global, rtrans, numboxes);
		#pragma oss taskwait

		normr = std::sqrt(rtrans_global);
		printf("Initial Residual = %g\n", normr);

		double brkdown_tol = std::numeric_limits<double>::epsilon();

		dbprintf("brkdown_tol = %g\n", brkdown_tol);

		for (int k = 1; k <= max_iter && normr > tolerance; ++k) {

			for (size_t i = 0; i < numboxes; ++i)
				waxpby_task(1.0, &r_array[i], beta, &p_array[i], &p_array[i]);


			// rtrans_global is here because of tw above
			normr = std::sqrt(rtrans_global);
			if (k % print_freq == 0 || k == max_iter)
				printf("Iteration = %d Residual = %g\n", k, normr);

			double p_ap_dot[numboxes], p_ap_dot_global = 0.0;
			double Ap2[numboxes];
			double p2[numboxes];

			// This creates tasks internally
			exchange_externals_all(A_array, p_array, numboxes, sing);
			for (size_t i = 0; i < numboxes; ++i)
				matvec_dot_task(&A_array[i], &p_array[i], &Ap_array[i],
				                &p_ap_dot[i], &Ap2[i], &p2[i]);

			reduce_sum_task(&p_ap_dot_global, p_ap_dot, numboxes);
			#pragma oss taskwait

			const double alpha = rtrans_global / p_ap_dot_global;

			dbprintf("iter: %d p_ap_dot: %g rtrans: %g alpha: %g\n",
			         k, p_ap_dot_global, rtrans_global, alpha);

			if (p_ap_dot_global < brkdown_tol) {

				int breakdown_global = 0;
				for (size_t i = 0; i < numboxes; ++i)
					breakdown_global += breakdown(p_ap_dot_global, Ap2[i], p2[i]);

				if (p_ap_dot_global < 0 || breakdown_global) {

					perror("miniFE::cg_solve ERROR, numerical breakdown!\n");

					dbprintf("ERROR: p_ap_dot_global = %lf && breakdown_global = %d\n",
					         p_ap_dot_global, breakdown_global);

					//update the timers before jumping out.
					delete [] r_array;
					delete [] p_array;
					delete [] Ap_array;

					return;
				} else {
					brkdown_tol = 0.1 * p_ap_dot_global;
				}
			}

			oldrtrans = rtrans_global;
			for (size_t i = 0; i < numboxes; ++i) {
				waxpby_task(1.0, &x_array[i], alpha, &p_array[i], &x_array[i]);
				waxpby_dot_task(1.0, &r_array[i], -alpha, &Ap_array[i], &r_array[i],
				            &rtrans[i]);
			}

			reduce_sum_task(&rtrans_global, rtrans, numboxes);
			#pragma oss taskwait

			beta = rtrans_global / oldrtrans;

			num_iters = k;
		} // for k

		delete [] r_array;
		delete [] p_array;
		delete [] Ap_array;
	}

}//namespace miniFE

#endif

