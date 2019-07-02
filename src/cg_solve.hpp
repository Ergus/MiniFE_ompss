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

#include <Vector_functions.hpp>
#include <mytimer.hpp>

#include <outstream.hpp>

namespace miniFE {

	inline int breakdown(double inner, const Vector &v, const Vector &w)
	{
		//This is code that was copied from Aztec, and originally written
		//by my hero, Ray Tuminaro.
		//
		//Assuming that inner = <v,w> (inner product of v and w),
		//v and w are considered orthogonal if
		//  |inner| < 100 * ||v||_2 * ||w||_2 * epsilon

		const double vnorm = std::sqrt(v.dot2());
		const double wnorm = std::sqrt(w.dot2());
		return std::abs(inner) <= 100 * vnorm * wnorm * std::numeric_limits<double>::epsilon();
	}


	template<typename MatrixType>
	void cg_solve_all(MatrixType *A_array, size_t numboxes,
	                  const Vector *b_array, Vector *x_array,
	                  int max_iter,
	                  const double tolerance,
	                  int &num_iters,
	                  double &normr,
	                  timer_type* my_cg_times)
	{
		timer_type t0 = 0, tWAXPY[numboxes] = {}, tDOT[numboxes] = {}, tMATVEC[numboxes] = {};
		timer_type total_time = mytimer();

		Vector *r_array = new Vector[numboxes];
		Vector *p_array = new Vector[numboxes];
		Vector *Ap_array = new Vector[numboxes];

		for (size_t i = 0; i < numboxes; ++i) {
			const int first_row_i = b_array[i].startIndex;
			const int local_nrows_i = A_array[i].nrows;
			const int Anumcols = A_array[i].num_cols;

			//TODO: Task here
			r_array[i].init(first_row_i, local_nrows_i);
			Ap_array[i].init(first_row_i, local_nrows_i);
			p_array[i].init(0, Anumcols);
		}

		normr = 0;
		double rtrans[numboxes] = {};
		double oldrtrans;
		double rtrans_global = 0.0;

		int print_freq = max_iter / 10;
		if (print_freq > 50)
			print_freq = 50;
		if (print_freq < 1)
			print_freq = 1;

		for (size_t i = 0; i < numboxes; ++i) {
			assert(A_array[i].has_local_indices);

			// TODO: task here
			{
				TICK();
				waxpby(1.0, x_array[i], 0.0, x_array[i], p_array[i]);
				TOCK(tWAXPY[i]);
			}
		}

		// This creates tasks internally
		exchange_externals_all(A_array, p_array, numboxes);

		for (size_t i = 0; i < numboxes; ++i) {
			// Tasks here
			// in A_array[i] (full)
			{
				TICK();
				A_array[i].matvec(p_array[i], Ap_array[i]);
				TOCK(tMATVEC[i]);

				TICK();
				waxpby(1.0, b_array[i], -1.0, Ap_array[i], r_array[i]);
				TOCK(tWAXPY[i]);


				TICK();
				rtrans[i] = r_array[i].dot2();
				TOCK(tDOT[i]);
			}
		}

		// TODO: taskwait here
		reduce_sum(rtrans_global, rtrans, numboxes);

		normr = std::sqrt(rtrans_global);
		std::cout << "Initial Residual = "<< normr << std::endl;

		double brkdown_tol = std::numeric_limits<double>::epsilon();

		#ifdef MINIFE_DEBUG
		std::ostream& os = outstream();
		os << "brkdown_tol = " << brkdown_tol << std::endl;
		#endif

		// TODO: this is the main part of the code, the tasks must go here,
		// not in the nested functions (if possible)
		for (int k = 1; k <= max_iter && normr > tolerance; ++k) {

			if (k == 1) {
				for (size_t i = 0; i < numboxes; ++i) {
					// TODO: Task here
					// in: r_array[i] + array
					// out p_array[i] + array
					TICK();
					waxpby(1.0, r_array[i], 0.0, r_array[i], p_array[i]);
					TOCK(tWAXPY[i]);
				}
			} else {

				oldrtrans = rtrans_global;
				for (size_t i = 0; i < numboxes; ++i)
				{
					// TODO: Task here
					// out: oldrtrans[i]
					// inout rtrans[i]
					TICK();
					rtrans[i] = r_array[i].dot2();
					TOCK(tDOT[i]);
				}


				// TODO: task_here
				// in rtrans[0;numboxes]
				// out rtrans_global
				{
					reduce_sum(rtrans_global, rtrans, numboxes);
				}

				for (size_t i = 0; i < numboxes; ++i) {
					// in rtrans_global
					// in oldrtrans
					{
						TICK();
						double beta = rtrans_global / oldrtrans;
						waxpby(1.0, r_array[i], beta, p_array[i], p_array[i]);
						TOCK(tWAXPY[i]);
					}
				}
			}

			// TODO: taskwait (or task with in rtrans_global, in k)
			{
				normr = std::sqrt(rtrans_global);
				if (k % print_freq == 0 || k == max_iter)
					std::cout << "Iteration = " << k << "   Residual = " << normr << std::endl;
			}

			double p_ap_dot[numboxes];
			double p_ap_dot_global = 0.0;

			// This creates tasks internally
			exchange_externals_all(A_array, p_array, numboxes);
			for (size_t i = 0; i < numboxes; ++i) {
				// TODO:  task here
				// in A_array (full)
				// in p_array[i]
				// out Ap_array[i]
				// out p_ap_dot[i]
				{
					TICK();
					A_array[i].matvec(p_array[i], Ap_array[i]);
					TOCK(tMATVEC[i]);

					TICK();
					p_ap_dot[i] = Ap_array[i].dot(p_array[i]);
					TOCK(tDOT[i]);
				}
			}

			// TODO: task here (or taskwait)
			// in k
			// in p_ap_dot[0; numboxes]
			// out p_ap_dot_global
			{
				reduce_sum(p_ap_dot_global, p_ap_dot, numboxes);
				#ifdef MINIFE_DEBUG
				os << "iter " << k << ", p_ap_dot = " << p_ap_dot_global;
				os.flush();
				#endif
			}

			if (p_ap_dot_global < brkdown_tol) {

				int breakdown_array[numboxes];
				int breakdown_global;

				for (size_t i = 0; i < numboxes; ++i) {
					// TODO: Task here
					// in Ap_array[i]
					// in p_array[i]
					// out breakdown[i]
					{
						breakdown_array[i] = breakdown(p_ap_dot_global, Ap_array[i], p_array[i]);
					}
				}

				// TODO taskwait here, because this must run locally.
				reduce_sum(breakdown_global, breakdown_array, numboxes);

				if (p_ap_dot_global < 0 || breakdown_global) {
					std::cerr << "miniFE::cg_solve ERROR, numerical breakdown!"<<std::endl;
					#ifdef MINIFE_DEBUG
					os << "ERROR, numerical breakdown!"<<std::endl;
					#endif
					//update the timers before jumping out.
					reduce_sum(my_cg_times[WAXPY], tWAXPY, numboxes);
					reduce_sum(my_cg_times[DOT], tDOT, numboxes);
					reduce_sum(my_cg_times[MATVEC], tMATVEC, numboxes);
					my_cg_times[TOTAL] = mytimer() - total_time;
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
				waxpby(1.0, x_array[i], alpha, p_array[i], x_array[i]);

				// Task here
				waxpby(1.0, r_array[i], -alpha, Ap_array[i], r_array[i]);
				TOCK(tWAXPY[i]);
			}
			num_iters = k;
		} // for k

		// TODO: taskwait here

		reduce_sum(my_cg_times[WAXPY], tWAXPY, numboxes);
		reduce_sum(my_cg_times[DOT], tDOT, numboxes);
		reduce_sum(my_cg_times[MATVEC], tMATVEC, numboxes);
		my_cg_times[TOTAL] = mytimer() - total_time;

		delete [] r_array;
		delete [] p_array;
		delete [] Ap_array;
	}

}//namespace miniFE

#endif

