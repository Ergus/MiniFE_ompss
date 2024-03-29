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
#include "CSRMatrix.hpp"
#include "Vector.hpp"
#include "simple_mesh_description.hpp"
#include "SparseMatrix_functions.hpp"
#include "verify_solution.hpp"
#include "ompss_utils.hpp"

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
		assert(numboxes > 0);

		singleton sing(numboxes);

		const int global_nx = global_box[0][1];
		const int global_ny = global_box[1][1];
		const int global_nz = global_box[2][1];

		timer_type t_total = 0.0, t0 = mytimer();

		std::cout << "Global Box: " << global_box << std::endl;

		double *ptr = (double *) rrd_malloc(1000 * numboxes * sizeof(double));
		double *ptr2 = (double *) rrd_malloc(100000 * numboxes * sizeof(double));
		double *ptr3 = (double *) rrd_malloc(100000 * numboxes * sizeof(double));
		double *ptr4 = (double *) rrd_malloc(100000 * numboxes * sizeof(double));
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

		Box *local_node_box_array = new Box[numboxes];
		for (size_t i = 0; i < numboxes; ++i) {
			local_node_box_array[i] = local_boxes_array[i];

			for (int j = 0; j < 3; ++j) {
				if (local_boxes_array[i][j][1] > local_boxes_array[i][j][0] &&
				    local_boxes_array[i][j][1] == global_box[j][1])
					++local_node_box_array[i][j][1];
			}
		}


		simple_mesh_description *mesh_array = new simple_mesh_description[numboxes];
		{
			timer_type t_mesh_fill = mytimer();
			init_mesh_all(mesh_array,
			              ptr,
			              &global_box,
			              local_boxes_array, \
			              local_node_box_array,
			              numboxes);

			REGISTER_ELAPSED_TIME(t_mesh_fill, t_total);
		}

		//Declare matrix object array
		CSRMatrix *A_array = new CSRMatrix[numboxes];
		{
			std::cout << "generating matrix structure..." << std::endl;
			timer_type t_gen_structure = mytimer();

			generate_matrix_structure_all(ptr3,
			                              A_array,
			                              mesh_array,
			                              &sing,
			                              numboxes);

			REGISTER_ELAPSED_TIME(t_gen_structure, t_total);

			ydoc.add("Matrix structure generation","");
			ydoc.get("Matrix structure generation")->add("Mat-struc-gen Time", t_gen_structure);
		}

		// Declare vector objects array
		Vector *b_array = new Vector[numboxes];
		Vector *x_array = new Vector[numboxes];

		{
			init_vector_all(b_array, A_array, &sing, numboxes);
			init_vector_all(x_array, A_array, &sing, numboxes);
		}

		//Assemble finite-element sub-matrices and sub-vectors into the global linear system:
		{
			std::cout << "assembling FE data..." << std::endl;
			timer_type t_fe_assembly = mytimer();

			assemble_FE_data_task(ptr4, mesh_array, A_array, b_array, numboxes);

			REGISTER_ELAPSED_TIME(t_fe_assembly, t_total);

			ydoc.add("FE assembly", "");
			ydoc.get("FE assembly")->add("FE assembly Time", t_fe_assembly);
		}


		//Now apply dirichlet boundary-conditions
		//(Apply the 0-valued surfaces first, then the 1-valued surface last.)
		{
			std::cout << "imposing Dirichlet BC..." << std::endl;
			timer_type t_dirbc_time  = mytimer();;

			// This dependencies are inout and individual per impose_dirichlet.
			for (size_t i = 0; i < numboxes; ++i) {
				simple_mesh_description *mesh_i = &mesh_array[i];
				CSRMatrix *A_i = &A_array[i];
				Vector *b_i = &b_array[i];

				impose_dirichlet_task(0.0,
				                      A_i,
				                      b_i,
				                      global_nx + 1, global_ny + 1, global_nz + 1,
				                      mesh_i->ompss2_bc_rows_0,
				                      mesh_i->bc_rows_0_size);

				impose_dirichlet_task(1.0,
				                      A_i,
				                      b_i,
				                      global_nx + 1, global_ny + 1, global_nz + 1,
				                      mesh_i->ompss2_bc_rows_1,
				                      mesh_i->bc_rows_1_size);

			}

			REGISTER_ELAPSED_TIME(t_dirbc_time, t_total);
		}



		//Transform global indices to local, set up communication information:
		{
			std::cout << "making matrix indices local..." << std::endl;
			timer_type t_make_local_time = mytimer();

			// TODO: weak task here
			make_local_matrix(ptr, A_array, &sing, numboxes);

			REGISTER_ELAPSED_TIME(t_make_local_time, t_total);
		}


		size_t global_nnz;
		compute_matrix_stats_task(A_array, numboxes, ydoc, &global_nnz);

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

		{
			cg_times[TOTAL] = mytimer();

			cg_solve_all(ptr2, A_array, numboxes, b_array, x_array,
			             max_iters, tol, num_iters, rnorm, &sing);

			REGISTER_ELAPSED_TIME(cg_times[TOTAL], t_total);

			printf("Final Resid Norm: %g solver_time %g iterations %d\n",
			       rnorm, cg_times[TOTAL], num_iters);
		}


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

		std::string title("CG solve");
		#endif

		ydoc.get("Global Run Parameters")->add("ScalarType", "double");
		ydoc.get("Global Run Parameters")->add("GlobalOrdinalType", "int");
		ydoc.get("Global Run Parameters")->add("LocalOrdinalType", "int");
		ydoc.add(title,"");
		ydoc.get(title)->add("Iterations", num_iters);
		ydoc.get(title)->add("Final Resid Norm", rnorm);

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

		const double total_op_flops = mv_flops + dot_flops + waxpy_flops;

		const double total_mflops = get_mflops(cg_times[TOTAL], total_op_flops);


		#if MINIFE_KERNELS == 0
		ydoc.get(title)->add("Total","");
		ydoc.get(title)->get("Total")->add("Total CG Time",cg_times[TOTAL]);
		ydoc.get(title)->get("Total")->add("Total CG Flops",total_op_flops);

		if (total_mflops >= 0)
			ydoc.get(title)->get("Total")->add("Total CG Mflops",total_mflops);
		else
			ydoc.get(title)->get("Total")->add("Total CG Mflops","inf");

		ydoc.get(title)->add("Time per iteration", cg_times[TOTAL] / num_iters);
		printf("Time per iteration: %g mflops: %g t_total: %g\n",
		       cg_times[TOTAL] / num_iters, total_mflops, t_total);
		#endif

		delete [] mesh_array;
		delete [] local_node_box_array;

		delete [] A_array;
		delete [] b_array;
		delete [] x_array;

		rrd_free(ptr, (1000 * numboxes * sizeof(double)));
		rrd_free(ptr2, (100000 * numboxes * sizeof(double)));
		rrd_free(ptr3, (100000 * numboxes * sizeof(double)));
		rrd_free(ptr4, (100000 * numboxes * sizeof(double)));

		return verify_result;
	}

}//namespace miniFE

#endif

