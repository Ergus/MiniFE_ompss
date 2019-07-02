
#ifndef _simple_mesh_description_hpp_
#define _simple_mesh_description_hpp_

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

#include <utils.hpp>

#include "ompss_utils.hpp"

#include "Box.hpp"
#include "box_utils.hpp"

namespace miniFE
{
	class simple_mesh_description {
	public:
		simple_mesh_description():
			ompss2_bc_rows_0(nullptr), ompss2_bc_rows_1(nullptr),
			bc_rows_0_size(0), bc_rows_1_size(0),
			ompss2_ids_to_rows(nullptr), ids_to_rows_size(0)
		{}

		~simple_mesh_description()
		{
			rrd_free(ompss2_bc_rows_0, bc_rows_0_size * sizeof(int));
			rrd_free(ompss2_bc_rows_1, bc_rows_1_size * sizeof(int));
			rrd_free(ompss2_ids_to_rows, ids_to_rows_size * sizeof(int));

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
			dbvprintf("Calling: %s, address %p size: %lu\n",
			          __PRETTY_FUNCTION__, ptr, sz);
			return rrl_free(ptr, sz);
		}

		void init(const Box &global_box_in,
		          const Box *local_boxes_array,    // Global boxes
		          const Box *local_node_box_array, // Boxes resized
		          size_t id, size_t total)
		{
			assert(id < total);

			// Internal class copies of the boxes
			global_box = global_box_in;
			local_box = local_boxes_array[id];

			// Needed for this function
			Box local_node_box(local_node_box_array[id]);

			const int max_node_x = global_box[0][1] + 1;
			const int max_node_y = global_box[1][1] + 1;
			const int max_node_z = global_box[2][1] + 1;

			// Local set and map. Will be copied later
			std::set<int> bc_rows_0, bc_rows_1;
			std::map<int, int> map_ids_to_rows =
				create_map_id_to_row(max_node_x, max_node_y, max_node_z,
				                     local_node_box_array, id, total);

			//As described in analytic_soln.hpp,
			//we will impose a 0 boundary-condition on faces x=0, y=0, z=0, y=1, z=1
			//we will impose a 1 boundary-condition on face x=1

			#ifdef MINIFE_DEBUG
			std::cout << std::endl;
			#endif
			const int X = 0;
			const int Y = 1;
			const int Z = 2;

			const int x1 = max_node_x - 1;
			const int y1 = max_node_y - 1;
			const int z1 = max_node_z - 1;

			//if we're on the x=0 face:
			if (global_box[X][0] == local_box[X][0]) {
				int miny = local_node_box[Y][0];
				int minz = local_node_box[Z][0];
				int maxy = local_node_box[Y][1];
				int maxz = local_node_box[Z][1];
				//expand y and z dimensions to include ghost layer
				if (local_node_box[Y][0] > 0) --miny;
				if (local_node_box[Z][0] > 0) --minz;
				if (local_node_box[Y][1] < max_node_y) ++maxy;
				if (local_node_box[Z][1] < max_node_z) ++maxz;

				for (int iz = minz; iz < maxz; ++iz) {
					for(int iy = miny; iy < maxy; ++iy) {
						const int nodeID = get_id(max_node_x, max_node_y, max_node_z,
						                          0, iy, iz);
						#ifdef MINIFE_DEBUG
						std::cout << "x = 0 BC, node " << nodeID
						          << ", (0," << iy << "," << iz << ")"
						          << std::endl;
						#endif
						bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
					}
				}
			}

			//if we're on the y=0 face:
			if (global_box[Y][0] == local_box[Y][0]) {
				int minx = local_node_box[X][0];
				int minz = local_node_box[Z][0];
				int maxx = local_node_box[X][1];
				int maxz = local_node_box[Z][1];
				//expand x and z dimensions to include ghost layer
				if (local_node_box[X][0] > 0) --minx;
				if (local_node_box[Z][0] > 0) --minz;
				if (local_node_box[X][1] < max_node_x) ++maxx;
				if (local_node_box[Z][1] < max_node_z) ++maxz;

				for(int iz=minz; iz<maxz; ++iz) {
					for(int ix=minx; ix<maxx; ++ix) {
						int nodeID = get_id(max_node_x, max_node_y, max_node_z,
						                    ix, 0, iz);
						#ifdef MINIFE_DEBUG
						std::cout << "y=0 BC, node " << nodeID
						          << ", (" << ix << ",0," << iz << ")"
						          << std::endl;
						#endif
						int row = find_row_for_id(nodeID, map_ids_to_rows);

						if (row < 0)
							std::cout<< "on the y==0 face (ix="<<ix
							         <<", iz="<<iz<<") "
							         << "ERROR: found negative row (" <<row
							         <<") for nodeID=" << nodeID <<std::endl;

						bc_rows_0.insert(row);
					}
				}
			}

			//if we're on the z=0 face:
			if (global_box[Z][0] == local_box[Z][0]) {
				int minx = local_node_box[X][0];
				int miny = local_node_box[Y][0];
				int maxx = local_node_box[X][1];
				int maxy = local_node_box[Y][1];
				//expand x and y dimensions to include ghost layer
				if (local_node_box[X][0] > 0) --minx;
				if (local_node_box[Y][0] > 0) --miny;
				if (local_node_box[X][1] < max_node_x) ++maxx;
				if (local_node_box[Y][1] < max_node_y) ++maxy;

				for(int iy = miny; iy < maxy; ++iy) {
					for(int ix = minx; ix < maxx; ++ix) {
						int nodeID = get_id(max_node_x, max_node_y, max_node_z,
							            ix, iy, 0);
						#ifdef MINIFE_DEBUG
						std::cout << "z = 0 BC, node " << nodeID
						          << ", (" << ix << "," << iy << ",0)"
						          << std::endl;
						#endif
						bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
					}
				}
			}

			//if we're on the x=1 face:
			if (global_box[X][1] == local_box[X][1]) {
				int minz = local_node_box[Z][0];
				int miny = local_node_box[Y][0];
				int maxz = local_node_box[Z][1];
				int maxy = local_node_box[Y][1];
				//expand z and y dimensions to include ghost layer
				if (local_node_box[Z][0] > 0) --minz;
				if (local_node_box[Y][0] > 0) --miny;
				if (local_node_box[Z][1] < max_node_z) ++maxz;
				if (local_node_box[Y][1] < max_node_y) ++maxy;

				for(int iy=miny; iy<maxy; ++iy) {
					for(int iz=minz; iz<maxz; ++iz) {
						int nodeID = get_id(max_node_x, max_node_y, max_node_z,
							            x1, iy, iz);
						int row = find_row_for_id(nodeID, map_ids_to_rows);
						#ifdef MINIFE_DEBUG
						std::cout << "x=1 BC, node " << nodeID
						          << ", row " << row
						          <<", ("<<x1<<","<<iy<<","<<iz<<")"
						          << std::endl;
						#endif
						bc_rows_1.insert(row);
					}
				}
			}

			//if we're on the y=1 face:
			if (global_box[Y][1] == local_box[Y][1]) {
				int minz = local_node_box[Z][0];
				int minx = local_node_box[X][0];
				int maxz = local_node_box[Z][1];
				int maxx = local_node_box[X][1];
				//expand z and x dimensions to include ghost layer
				if (local_node_box[Z][0] > 0) --minz;
				if (local_node_box[X][0] > 0) --minx;
				if (local_node_box[Z][1] < max_node_z) ++maxz;
				if (local_node_box[X][1] < max_node_x) ++maxx;

				for(int ix=minx; ix<maxx; ++ix) {
					for(int iz=minz; iz<maxz; ++iz) {
						int nodeID = get_id(max_node_x, max_node_y, max_node_z,
						                         ix, y1, iz);
						#ifdef MINIFE_DEBUG
						std::cout << "y=1 BC, node " << nodeID
						          << ", (" << ix << "," << y1 << "," << iz << ")"
						          << std::endl;
						#endif
						bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
					}
				}
			}

			//if we're on the z=1 face:
			if (global_box[Z][1] == local_box[Z][1]) {
				int miny = local_node_box[Y][0];
				int minx = local_node_box[X][0];
				int maxy = local_node_box[Y][1];
				int maxx = local_node_box[X][1];
				//expand x and y dimensions to include ghost layer
				if (local_node_box[Y][0] > 0) --miny;
				if (local_node_box[X][0] > 0) --minx;
				if (local_node_box[Y][1] < max_node_y) ++maxy;
				if (local_node_box[X][1] < max_node_x) ++maxx;

				for(int ix = minx; ix < maxx; ++ix) {
					for(int iy = miny; iy < maxy; ++iy) {
						int nodeID =
							get_id(max_node_x, max_node_y, max_node_z,
							       ix, iy, z1);
						#ifdef MINIFE_DEBUG
						std::cout << "z=1 BC, node " << nodeID
						          << ", (" << ix << "," << iy << "," << z1 << ")"
						          << std::endl;
						#endif
						bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
					}
				}
			}

			// Copy the sets to ompss supported ones.
			bc_rows_0_size = stl_to_global_task(&ompss2_bc_rows_0, bc_rows_0);
			bc_rows_1_size = stl_to_global_task(&ompss2_bc_rows_1, bc_rows_1);

			ids_to_rows_size =
				stl_to_global_task(&ompss2_ids_to_rows, map_ids_to_rows);

		}

		int map_id_to_row(int id) const
		{
			return find_row_for_id(id,
			                       ompss2_ids_to_rows,
			                       &ompss2_ids_to_rows[ids_to_rows_size]);
		}

		int *ompss2_bc_rows_0;
		int *ompss2_bc_rows_1;
		size_t bc_rows_0_size, bc_rows_1_size;

		std::pair<int,int> *ompss2_ids_to_rows;
		size_t ids_to_rows_size;

		Box global_box;
		Box local_box;
	}; //class simple_mesh_description

}//namespace miniFE

#endif
