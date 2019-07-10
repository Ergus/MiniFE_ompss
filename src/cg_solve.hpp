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
		in(*A)							\
		in(Arecv_ptr[0; Anrecv_neighbors])			\
		in(Arecv_length[0; Anrecv_neighbors])			\
		out(x_start[0; Anexternals])				\
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
			ompss_memcpy_task(local_iter, Arecv_ptr[i], Arecv_length[i] * sizeof(double));

			local_iter += Arecv_length[i];

			// boundary check
			assert(local_iter - x_start <= Anexternals);
		}

	}


	void exchange_externals_all(CSRMatrix *A_array, Vector *x_array, size_t numboxes,
	                            const singleton *sing)
	{

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

	int breakdown(double inner, const Vector *v, const Vector *w)
	{
		//This is code that was copied from Aztec, and originally written
		//by my hero, Ray Tuminaro.
		//
		//Assuming that inner = <v,w> (inner product of v and w),
		//v and w are considered orthogonal if
		//  |inner| < 100 * ||v||_2 * ||w||_2 * epsilon

		double vnorm = 0;
		double wnorm = 0;

		dot2_task(v, &vnorm);
		dot2_task(w, &wnorm);
		#pragma oss taskwait

		vnorm = std::sqrt(vnorm);
		wnorm = std::sqrt(wnorm);

		return std::abs(inner) <= 100 * vnorm * wnorm * std::numeric_limits<double>::epsilon();
	}


	void cg_solve_all(CSRMatrix *A_array, size_t numboxes,
	                  const Vector *b_array, Vector *x_array,
	                  int max_iter,
	                  const double tolerance,
	                  int &num_iters,
	                  double &normr,
	                  timer_type* my_cg_times,
	                  singleton *sing)
	{

		timer_type t0 = 0, tWAXPY[numboxes], tDOT[numboxes], tMATVEC[numboxes];
		timer_type total_time = mytimer();

		Vector *r_array = new Vector[numboxes];
		Vector *p_array = new Vector[numboxes];
		Vector *Ap_array = new Vector[numboxes];

		{
			int start[numboxes];
			int length[numboxes];

			int p_start[numboxes];
			int p_length[numboxes];

			for (size_t i = 0; i < numboxes; ++i) {
				start[i] = A_array[i].first_row;
				length[i] = A_array[i].nrows;

				p_start[i] = 0;
				p_length[i] = A_array[i].num_cols;

				tWAXPY[i] = 0.0;
				tDOT[i] = 0.0;
				tMATVEC[i] = 0.0;
			}

			init_vector_all(r_array, sing, numboxes, start, length);
			init_vector_all(Ap_array, sing, numboxes, start, length);
			init_vector_all(p_array, sing, numboxes, p_start, p_length);
		}


		for (size_t i = 0; i < numboxes; ++i) {
		}


		normr = 0;
		double rtrans[numboxes];
		double oldrtrans;
		double rtrans_global = 0.0;

		int print_freq = max_iter / 10;
		if (print_freq > 50)
			print_freq = 50;
		if (print_freq < 1)
			print_freq = 1;

		//TICK();
		for (size_t i = 0; i < numboxes; ++i) {
			rtrans[i] = 0;
			assert(A_array[i].has_local_indices);

			waxpby_task(1.0, &x_array[i], 0.0, &x_array[i], &p_array[i]);

			#ifdef VERBOSE
			std::string filename = "VERB_p_1_"  + std::to_string(i) + ".verb";
			write_task(filename ,p_array[i], i);
			#endif

		}
		//TOCK(tWAXPY[0]);

		// This creates tasks internally
		exchange_externals_all(A_array, p_array, numboxes, sing);

		for (size_t i = 0; i < numboxes; ++i) {
			// Tasks here
			// in A_array[i] (full)
			//TICK();
			#ifdef VERBOSE
			std::string filename1 = "VERB_p_2_"  + std::to_string(i) + ".verb";
			write_task(filename1 ,p_array[i], i);
			std::string filename2 = "VERB_A_2_"  + std::to_string(i) + ".verb";
			write_task(filename2 ,A_array[i], i);
			#endif

			matvec_task(&A_array[i], &p_array[i], &Ap_array[i]);
			//TOCK(tMATVEC[i]);

			//TICK();
			waxpby_task(1.0, &b_array[i], -1.0, &Ap_array[i], &r_array[i]);
			//TOCK(tWAXPY[i]);


			//TICK();
			dot2_task(&r_array[i], &rtrans[i]);
			//TOCK(tDOT[i]);
		}

		// TODO: taskwait here
		reduce_sum_task(&rtrans_global, rtrans, numboxes);
		#pragma oss taskwait

		normr = std::sqrt(rtrans_global);
		std::cout << "Initial Residual = "<< normr << std::endl;

		double brkdown_tol = std::numeric_limits<double>::epsilon();

		#ifdef MINIFE_DEBUG
		std::ostream& os = outstream();
		os << "brkdown_tol = " << brkdown_tol << std::endl;
		#endif

		for (int k = 1; k <= max_iter && normr > tolerance; ++k) {

			if (k == 1) {
				for (size_t i = 0; i < numboxes; ++i)
					waxpby_task(1.0, &r_array[i], 0.0, &r_array[i], &p_array[i]);
			} else {

				oldrtrans = rtrans_global;
				for (size_t i = 0; i < numboxes; ++i)
					dot2_task(&r_array[i], &rtrans[i]);


				reduce_sum_task(&rtrans_global, rtrans, numboxes);
				#pragma oss taskwait

				const double beta = rtrans_global / oldrtrans;

				for (size_t i = 0; i < numboxes; ++i) {
					{
						//TICK();
						waxpby_task(1.0, &r_array[i], beta, &p_array[i], &p_array[i]);
						//TOCK(tWAXPY[i]);
					}
				}
			}

			// rtrans_global is here because of tw above
			normr = std::sqrt(rtrans_global);
			if (k % print_freq == 0 || k == max_iter)
				std::cout << "Iteration = " << k << "   Residual = " << normr << std::endl;

			double p_ap_dot[numboxes];
			double p_ap_dot_global = 0.0;

			// This creates tasks internally
			exchange_externals_all(A_array, p_array, numboxes, sing);
			for (size_t i = 0; i < numboxes; ++i) {
				//TICK();
				matvec_task(&A_array[i], &p_array[i], &Ap_array[i]);
				//TOCK(tMATVEC[i]);

				//TICK();
				dot_task(&Ap_array[i], &p_array[i], &p_ap_dot[i]);
				//TOCK(tDOT[i]);
			}

			reduce_sum_task(&p_ap_dot_global, p_ap_dot, numboxes);
			#pragma oss taskwait

			#ifdef MINIFE_DEBUG
			os << "iter " << k << ", p_ap_dot = " << p_ap_dot_global;
			os.flush();
			#endif

			if (p_ap_dot_global < brkdown_tol) {

				int breakdown_array[numboxes];
				int breakdown_global;

				for (size_t i = 0; i < numboxes; ++i) {
					breakdown_array[i] =
						breakdown(p_ap_dot_global, &Ap_array[i], &p_array[i]);
				}

				// TODO taskwait here, because this must run locally.
				reduce_sum_task(&breakdown_global, breakdown_array, numboxes);
				#pragma oss taskwait

				if (p_ap_dot_global < 0 || breakdown_global) {

					fprintf(stderr, "miniFE::cg_solve ERROR, numerical breakdown!\n");

					dbprintf("ERROR: p_ap_dot_global = %lf && breakdown_global = %d\n",
					         p_ap_dot_global, breakdown_global);

					//update the timers before jumping out.
					reduce_sum_task(&my_cg_times[WAXPY], tWAXPY, numboxes);
					reduce_sum_task(&my_cg_times[DOT], tDOT, numboxes);
					reduce_sum_task(&my_cg_times[MATVEC], tMATVEC, numboxes);
					my_cg_times[TOTAL] = mytimer() - total_time;

					#pragma oss taskwait
					delete [] r_array;
					delete [] p_array;
					delete [] Ap_array;

					return;
				} else {
					brkdown_tol = 0.1 * p_ap_dot_global;
				}
			}

			const double alpha = rtrans_global / p_ap_dot_global;
			#ifdef MINIFE_DEBUG
			os << ", rtrans = " << rtrans_global << ", alpha = " << alpha << std::endl;
			#endif

			for (size_t i = 0; i < numboxes; ++i) {
				// Task here
				TICK();
				waxpby_task(1.0, &x_array[i], alpha, &p_array[i], &x_array[i]);

				// Task here
				waxpby_task(1.0, &r_array[i], -alpha, &Ap_array[i], &r_array[i]);
				TOCK(tWAXPY[i]);
			}

			#pragma oss taskwait
			num_iters = k;
		} // for k

		reduce_sum_task(&my_cg_times[WAXPY], tWAXPY, numboxes);
		reduce_sum_task(&my_cg_times[DOT], tDOT, numboxes);
		reduce_sum_task(&my_cg_times[MATVEC], tMATVEC, numboxes);
		#pragma oss taskwait

		my_cg_times[TOTAL] = mytimer() - total_time;

		delete [] r_array;
		delete [] p_array;
		delete [] Ap_array;
	}

}//namespace miniFE

#endif

