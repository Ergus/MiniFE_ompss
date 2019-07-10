#ifndef _perform_element_loop_hpp_
#define _perform_element_loop_hpp_

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

#include <BoxIterator.hpp>
#include <simple_mesh_description.hpp>
#include "CSRMatrix.hpp"
#include <SparseMatrix_functions.hpp>
#include <box_utils.hpp>
#include <Hex8_box_utils.hpp>
#include <Hex8_ElemData.hpp>

namespace miniFE
{

	template <typename MatrixType>
	void perform_element_loop(const simple_mesh_description &mesh,
	                          const Box &local_elem_box,
	                          MatrixType &A,
	                          Vector &b)
	{
		int global_elems_x = mesh.global_box[0][1];
		int global_elems_y = mesh.global_box[1][1];
		int global_elems_z = mesh.global_box[2][1];

		//We will iterate the local-element-box (local portion of the mesh), and
		//get element-IDs in preparation for later assembling the FE operators
		//into the global sparse linear-system.

		int num_elems = local_elem_box.get_num_ids();
		std::vector<int> elemIDs(num_elems);

		BoxIterator iter = BoxIterator::begin(local_elem_box);
		BoxIterator end  = BoxIterator::end(local_elem_box);

		for(size_t i = 0; iter != end; ++iter, ++i) {
			elemIDs[i] = get_id(global_elems_x, global_elems_y, global_elems_z,
			                    iter.x, iter.y, iter.z);

			dbv2printf("elem ID %d (%d, %d, %d)\n", elemIDs[i], iter.x, iter.y, iter.z);
		}

		//Now do the actual finite-element assembly loop:

		ElemData elem_data;

		compute_gradient_values(elem_data.grad_vals);

		timer_type t_gn = 0, t_ce = 0, t_si = 0;
		timer_type t0 = 0;
		for (size_t i = 0; i < elemIDs.size(); ++i) {
			//Given an element-id, populate elem_data with the
			//element's node_ids and nodal-coords:

			TICK();
			elem_data.get_elem_nodes_and_coords(mesh, elemIDs[i]);
			TOCK(t_gn);

			//Next compute element-diffusion-matrix and element-source-vector:

			TICK();
			elem_data.compute_element_matrix_and_vector();
			TOCK(t_ce);

			//Now assemble the (dense) element-matrix and element-vector into the
			//global sparse linear system:

			TICK();
			sum_into_global_linear_system(elem_data, A, b);
			TOCK(t_si);
		}

		dbvprintf("get-nodes: %lf\n", t_gn);
		dbvprintf("compute-elems: %lf\n", t_ce);
		dbvprintf("sum-in: %lf\n", t_si);

	}



}//namespace miniFE

#endif

