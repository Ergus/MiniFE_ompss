#ifndef _driver_hpp_
#define _driver_hpp_

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
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "box_utils.hpp"
#include "Vector.hpp"

#include "CSRMatrix.hpp"

#include "simple_mesh_description.hpp"

#include "SparseMatrix_functions.hpp"


#include "verify_solution.hpp"

#include "singleton.hpp"
#include "compute_matrix_stats.hpp"
#include "make_local_matrix.hpp"
#include "imbalance.hpp"
#include "cg_solve.hpp"
#if MINIFE_KERNELS != 0
#include "time_kernels.hpp"
#endif
#include "outstream.hpp"
#include "utils.hpp"
#include "mytimer.hpp"
#include "YAML_Doc.hpp"

//This program assembles finite-element matrices into a global matrix and
//vector, then solves the linear-system using Conjugate Gradients.
//Each finite-element is a hexahedron with 8 vertex-nodes.
//
//Notes:
//- In finite-element terms, the box dimensions are in elements, not nodes.
//  In other words, a 2x2x2 box describes 8 elements, each of which has 8 nodes,
//  so it is a 3x3x3 node domain (27 nodes).
//  The assembled linear system will have 1 equation for each finite element node.
//
//- The coordinate origin is at the corner of the global box where x=0,
//  y=0, z=0, and the box extends along the positive x-axis, positive y-axis,
//  and the positive z-axis.
//
//- Some aspects of matrix-structure generation and finite-element assembly
//  are convenient to do using global node identifiers.
//  A global identifier for each node is obtained from coordinates plus
//  global box dimensions. See the function 'get_id' in box_utils.hpp.
//
//- Each node corresponds to a row in the matrix. The RCB partitioning method
//  we use to split the global box among processors results in some
//  processors owning non-contiguous blocks of global node identifiers.
//  Since it is convenient for matrices and vectors to store contiguously-
//  numbered blocks of rows, we map global node identifiers to a separate
//  space of row numbers such that each processor's nodes correspond to a
//  contiguous block of row numbers.
//

namespace miniFE
{

	inline int driver(const Box &global_box,
	                  Box *local_boxes_array, const size_t numboxes,
	                  Parameters &params, YAML_Doc &ydoc)
	{
		const int global_nx = global_box[0][1];
		const int global_ny = global_box[1][1];
		const int global_nz = global_box[2][1];

		// TODO: imbalance is implemented but not tested yet
		// if (params.load_imbalance > 0)
		// 	add_imbalance<GlobalOrdinal>(global_box, local_boxes, numboxes,
		// 	                             params.load_imbalance, ydoc);

		// float largest_imbalance = 0, std_dev = 0;
		// compute_imbalance<GlobalOrdinal>(global_box, my_box, largest_imbalance,
		//                                  std_dev, ydoc, true);

		//Create a representation of the mesh:
		//Note that 'simple_mesh_description' is a virtual or conceptual
		//mesh that doesn't actually store mesh data.

		timer_type t_start = mytimer();
		timer_type t0 = mytimer();

		// Tasks here (is not really needed, just to paralelize)
		Box *local_node_box_array = new Box[numboxes];
		for (size_t i = 0; i < numboxes; ++i) {
			local_node_box_array[i] = local_boxes_array[i];

			for (int j = 0; j < 3; ++j) {
				if (local_boxes_array[i][j][1] > local_boxes_array[i][j][0] &&
				    local_boxes_array[i][j][1] == global_box[j][1])
					++local_node_box_array[i][j][1];
			}
		}

		global_box.write();

		simple_mesh_description *mesh_array = new simple_mesh_description[numboxes];

		for (size_t i = 0; i < numboxes; ++i) {
			init_mesh_task(&mesh_array[i],
			               &global_box,
			               local_boxes_array, \
			               local_node_box_array,
			               i,
			               numboxes);

		}
		//#pragma oss taskwait


		timer_type mesh_fill = mytimer() - t0;
		timer_type t_total = mytimer() - t_start;

		std::cout << mesh_fill << "s, total time: " << t_total << std::endl;

		//Declare matrix object array
		CSRMatrix *A_array = new CSRMatrix[numboxes];
		{
			std::cout << "generating matrix structure..." << std::endl;
			timer_type gen_structure = mytimer();

			for (size_t i = 0; i < numboxes; ++i) {
				generate_matrix_structure_task(&A_array[i],
				                               &mesh_array[i],
				                               mesh_array[i].ompss2_ids_to_rows,
				                               mesh_array[i].ids_to_rows_size,
				                               i);

			}

			//#pragma oss taskwait
			REGISTER_ELAPSED_TIME(gen_structure, t_total);

			ydoc.add("Matrix structure generation","");
			ydoc.get("Matrix structure generation")->add("Mat-struc-gen Time", gen_structure);
		}

		assert(numboxes > 0);
		singleton sing(numboxes);

		// Declare vector objects array
		Vector *b_array = new Vector[numboxes];
		Vector *x_array = new Vector[numboxes];


		// TODO: Task Here gene

		{
			int Astart[numboxes];
			int Alength[numboxes];

			for (size_t i = 0; i < numboxes; ++i) {
				Astart[i] = A_array[i].first_row;
				Alength[i] = A_array[i].nrows;
			}

			init_vector_all(b_array, &sing, numboxes, Astart, Alength);
			init_vector_all(x_array, &sing, numboxes, Astart, Alength);
		}

		//Assemble finite-element sub-matrices and sub-vectors into the global linear system:
		{
			std::cout << "assembling FE data..." << std::endl;
			timer_type fe_assembly = mytimer();

			for (size_t i = 0; i < numboxes; ++i) {
				simple_mesh_description *mesh_i = &mesh_array[i];
				CSRMatrix *A_i = &A_array[i];
				Vector *b_i = &b_array[i];

				assemble_FE_data_task(i,
				                      mesh_i,
				                      mesh_i->ompss2_ids_to_rows,
				                      mesh_i->ids_to_rows_size,
				                      A_i,
				                      A_i->rows,
				                      A_i->row_offsets,
				                      A_i->nrows,
				                      A_i->packed_cols,
				                      A_i->packed_coefs,
				                      A_i->nnz,
				                      b_i,
				                      b_i->coefs,
				                      b_i->local_size);
			}
			#pragma oss taskwait

			REGISTER_ELAPSED_TIME(fe_assembly, t_total);

			ydoc.add("FE assembly", "");
			ydoc.get("FE assembly")->add("FE assembly Time",fe_assembly);
		}


		//Now apply dirichlet boundary-conditions
		//(Apply the 0-valued surfaces first, then the 1-valued surface last.)
		{
			std::cout << "imposing Dirichlet BC..." << std::endl;
			timer_type dirbc_time  = mytimer();;

			// This dependencies are inout and individual per impose_dirichlet.
			for (size_t i = 0; i < numboxes; ++i) {
				simple_mesh_description *mesh_i = &mesh_array[i];
				CSRMatrix *A_i = &A_array[i];
				Vector *b_i = &b_array[i];

				impose_dirichlet_task(0.0,
				                 A_i,
				                 A_i->rows,
				                 A_i->row_offsets,
				                 A_i->nrows,
				                 A_i->packed_cols,
				                 A_i->packed_coefs,
				                 A_i->nnz,
				                 b_i,
				                 b_i->coefs,
				                 b_i->local_size,
				                 global_nx + 1, global_ny + 1, global_nz + 1,
				                 mesh_i->ompss2_bc_rows_0,
				                 mesh_i->bc_rows_0_size);

				impose_dirichlet_task(1.0,
				                 A_i,
				                 A_i->rows,
				                 A_i->row_offsets,
				                 A_i->nrows,
				                 A_i->packed_cols,
				                 A_i->packed_coefs,
				                 A_i->nnz,
				                 b_i,
				                 b_i->coefs,
				                 b_i->local_size,
				                 global_nx + 1, global_ny + 1, global_nz + 1,
				                 mesh_i->ompss2_bc_rows_1,
				                 mesh_i->bc_rows_1_size);

			}

			#pragma oss taskwait

			REGISTER_ELAPSED_TIME(dirbc_time, t_total);
		}



		//Transform global indices to local, set up communication information:
		{
			std::cout << "making matrix indices local..." << std::endl;
			timer_type make_local_time = mytimer();;

			// TODO: weak task here
			make_local_matrix(A_array, &sing, numboxes);
			REGISTER_ELAPSED_TIME(make_local_time, t_total);
		}

		#ifdef MINIFE_DEBUG
		A.write_matrix("A_local.mtx");
		b.write_vector("b_local.vec");
		#endif

		size_t global_nnz;
		compute_matrix_stats(A_array, numboxes, ydoc, global_nnz);

		//Prepare to perform conjugate gradient solve:

		const int max_iters = 200;
		int num_iters = 0;
		double rnorm = 0.0;
		const double tol = std::numeric_limits<double>::epsilon();

		timer_type cg_times[NUM_TIMERS];

		t_total = mytimer() - t_start;

		int verify_result = 0;

		#if MINIFE_KERNELS != 0
		std::cout.width(30);
		std::cout << "Starting kernel timing loops ..." << std::endl;

		max_iters = 500;
		x.coefs[0] = 0.9;

		time_kernels(A, b, x,
		             matvec<MatrixType>(),
		             max_iters, rnorm, cg_times);

		num_iters = max_iters;
		std::string title("Kernel timings");
		#else

		std::cout << "Starting CG solver ... " << std::endl;

		cg_solve_all(A_array, numboxes, b_array, x_array,
		             max_iters, tol, num_iters, rnorm, cg_times, &sing);

		std::cout << "Final Resid Norm: " << rnorm << std::endl;

		if (params.verify_solution > 0) {
			#ifndef NDEBUG
			bool verify_whole_domain = true;
			std::cout << "verifying solution..." << std::endl;
			#else
			bool verify_whole_domain = false;
			std::cout << "verifying solution at ~ (0.5, 0.5, 0.5) ..." << std::endl;
			#endif

			verify_result = verify_solution(&global_box, local_node_box_array,
			                                x_array, numboxes, 0.06, verify_whole_domain);
		}

		#ifdef MINIFE_DEBUG
		write_vector("x.vec", x);
		#endif // MINIFE_DEBUG
		std::string title("CG solve");
		#endif

		ydoc.get("Global Run Parameters")->add("ScalarType","double");
		ydoc.get("Global Run Parameters")->add("GlobalOrdinalType","int");
		ydoc.get("Global Run Parameters")->add("LocalOrdinalType","int");
		ydoc.add(title,"");
		ydoc.get(title)->add("Iterations",num_iters);
		ydoc.get(title)->add("Final Resid Norm",rnorm);

		const int global_nrows = global_nx * global_ny * global_nz;

		//flops-per-mv, flops-per-dot, flops-per-waxpy:
		double mv_flops = global_nnz * 2.0;
		double dot_flops = global_nrows * 2.0;
		double waxpy_flops = global_nrows * 3.0;

		#if MINIFE_KERNELS == 0
		//if MINIFE_KERNELS == 0 then we did a CG solve, and in that case
		//there were num_iters+1 matvecs, num_iters*2 dots, and num_iters*3+2 waxpys.
		mv_flops *= (num_iters + 1);
		dot_flops *= (2 * num_iters);
		waxpy_flops *= (3 * num_iters+2);
		#else
		//if MINIFE_KERNELS then we did one of each operation per iteration.
		mv_flops *= num_iters;
		dot_flops *= num_iters;
		waxpy_flops *= num_iters;
		#endif

		const double total_flops = mv_flops + dot_flops + waxpy_flops;

		double mv_mflops = -1;
		if (cg_times[MATVEC] > 1.e-4)
			mv_mflops = 1.e-6 * (mv_flops/cg_times[MATVEC]);

		double dot_mflops = -1;
		if (cg_times[DOT] > 1.e-4)
			dot_mflops = 1.e-6 * (dot_flops/cg_times[DOT]);

		double waxpy_mflops = -1;
		if (cg_times[WAXPY] > 1.e-4)
			waxpy_mflops = 1.e-6 *  (waxpy_flops/cg_times[WAXPY]);

		double total_mflops = -1;
		if (cg_times[TOTAL] > 1.e-4)
			total_mflops = 1.e-6 * (total_flops/cg_times[TOTAL]);

		ydoc.get(title)->add("WAXPY Time",cg_times[WAXPY]);
		ydoc.get(title)->add("WAXPY Flops",waxpy_flops);
		if (waxpy_mflops >= 0)
			ydoc.get(title)->add("WAXPY Mflops",waxpy_mflops);
		else
			ydoc.get(title)->add("WAXPY Mflops","inf");

		ydoc.get(title)->add("DOT Time",cg_times[DOT]);
		ydoc.get(title)->add("DOT Flops",dot_flops);
		if (dot_mflops >= 0)
			ydoc.get(title)->add("DOT Mflops",dot_mflops);
		else
			ydoc.get(title)->add("DOT Mflops","inf");

		ydoc.get(title)->add("MATVEC Time",cg_times[MATVEC]);
		ydoc.get(title)->add("MATVEC Flops",mv_flops);
		if (mv_mflops >= 0)
			ydoc.get(title)->add("MATVEC Mflops",mv_mflops);
		else
			ydoc.get(title)->add("MATVEC Mflops","inf");

		#ifdef MINIFE_FUSED
		ydoc.get(title)->add("MATVECDOT Time",cg_times[MATVECDOT]);
		ydoc.get(title)->add("MATVECDOT Flops",mv_flops);
		if (mv_mflops >= 0)
			ydoc.get(title)->add("MATVECDOT Mflops",mv_mflops);
		else
			ydoc.get(title)->add("MATVECDOT Mflops","inf");
		#endif

		#if MINIFE_KERNELS == 0
		ydoc.get(title)->add("Total","");
		ydoc.get(title)->get("Total")->add("Total CG Time",cg_times[TOTAL]);
		ydoc.get(title)->get("Total")->add("Total CG Flops",total_flops);
		if (total_mflops >= 0)
			ydoc.get(title)->get("Total")->add("Total CG Mflops",total_mflops);
		else
			ydoc.get(title)->get("Total")->add("Total CG Mflops","inf");
		ydoc.get(title)->add("Time per iteration",cg_times[TOTAL]/num_iters);
		#endif

		#pragma oss taskwait
		delete [] mesh_array;
		delete [] local_node_box_array;

		delete [] A_array;
		delete [] b_array;
		delete [] x_array;

		return verify_result;
	}

}//namespace miniFE

#endif

