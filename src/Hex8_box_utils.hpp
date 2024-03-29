#ifndef _Hex8_box_utils_hpp_
#define _Hex8_box_utils_hpp_

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

#include <stdexcept>

#include "box_utils.hpp"
#include "ElemData.hpp"
#include "simple_mesh_description.hpp"
#include "Hex8.hpp"
#include "ompss_utils.hpp"

namespace miniFE
{


	template<typename GlobalOrdinal>
	void get_hex8_node_ids(int nx, int ny,
	                       GlobalOrdinal node0,
	                       GlobalOrdinal* elem_node_ids)
	{
		//Given box dimensions nx and ny, and a starting node
		//(local-node-0 for a hex8), compute the other nodes
		//of the hex8 using the exodus ordering convention.
		elem_node_ids[0] = node0;
		elem_node_ids[1] = node0 + 1;
		elem_node_ids[2] = node0 + nx + 1;
		elem_node_ids[3] = node0 + nx;
		elem_node_ids[4] = node0 +     nx*ny;
		elem_node_ids[5] = node0 + 1 + nx*ny;
		elem_node_ids[6] = node0 + nx + nx*ny + 1;
		elem_node_ids[7] = node0 + nx + nx*ny;
	}

	template<typename Scalar>
	void get_hex8_node_coords_3d(Scalar x, Scalar y, Scalar z,
	                             Scalar hx, Scalar hy, Scalar hz,
	                             Scalar* elem_node_coords)
	{
		//Input: x,y,z are the coordinates of local-node 0 for a Hex8.
		//'hx', 'hy', 'hz' are the lengths of the sides of the element
		//in each direction.

		elem_node_coords[0] = x;
		elem_node_coords[1] = y;
		elem_node_coords[2] = z;

		elem_node_coords[3] = x + hx;
		elem_node_coords[4] = y;
		elem_node_coords[5] = z;

		elem_node_coords[6] = x + hx;
		elem_node_coords[7] = y + hy;
		elem_node_coords[8] = z;

		elem_node_coords[9]  = x;
		elem_node_coords[10] = y + hy;
		elem_node_coords[11] = z;

		elem_node_coords[12] = x;
		elem_node_coords[13] = y;
		elem_node_coords[14] = z + hz;

		elem_node_coords[15] = x + hx;
		elem_node_coords[16] = y;
		elem_node_coords[17] = z + hz;

		elem_node_coords[18] = x + hx;
		elem_node_coords[19] = y + hy;
		elem_node_coords[20] = z + hz;

		elem_node_coords[21] = x;
		elem_node_coords[22] = y + hy;
		elem_node_coords[23] = z + hz;
	}


	inline void get_elem_nodes_and_coords(const simple_mesh_description &mesh, int elemID,
	                                      int *node_ords, double *node_coords)
	{
		const int global_nodes_x = mesh.global_box[0][1] + 1;
		const int global_nodes_y = mesh.global_box[1][1] + 1;
		const int global_nodes_z = mesh.global_box[2][1] + 1;

		if (elemID < 0) //I don't think this can happen, but check for the sake of paranoia...
			throw std::runtime_error("get_elem_nodes_and_coords ERROR, negative elemID");

		int elem_int_x, elem_int_y, elem_int_z;
		get_int_coords(elemID,
		               global_nodes_x - 1, global_nodes_y - 1, global_nodes_z - 1,
		               elem_int_x, elem_int_y, elem_int_z);

		int nodeID = get_id(global_nodes_x, global_nodes_y, global_nodes_z,
		                    elem_int_x, elem_int_y, elem_int_z);

		dbv2printf("elemID: %d, nodeID: %d\n", elemID, nodeID);

		get_hex8_node_ids(global_nodes_x, global_nodes_y, nodeID, node_ords);

		//Map node-IDs to rows because each processor may have a non-contiguous block of
		//node-ids, but needs a contiguous block of row-numbers:

		dbv2printf("elem %d nodes: ", elemID);
		for(int i = 0; i < Hex8::numNodesPerElem; ++i) {
			dbv2printf("%d ", node_ords[i]);

			node_ords[i] = mesh.map_id_to_row(node_ords[i]);
		}
		dbv2printf("\n");

		double ix, iy, iz;

		get_coords(nodeID, global_nodes_x,global_nodes_y,global_nodes_z,
		           ix,iy,iz);

		const double hx = 1.0 / mesh.global_box[0][1];
		const double hy = 1.0 / mesh.global_box[1][1];
		const double hz = 1.0 / mesh.global_box[2][1];

		get_hex8_node_coords_3d(ix, iy, iz, hx, hy, hz, node_coords);

		#if VERBOSE == 2
		int offset = 0;
		for(int i = 0; i < Hex8::numNodesPerElem; ++i) {
			dbv2printf("(%lf, %lf,%lf)",
			           node_coords[offset++],
			           node_coords[offset++],
			           node_coords[offset++]);
		}
		std::cout << std::endl;
		#endif
	}


}//namespace miniFE

#endif
