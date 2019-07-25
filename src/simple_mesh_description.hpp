
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
#include <cstdio>

#include "ompss_utils.hpp"

#include "Box.hpp"
#include "box_utils.hpp"

namespace miniFE
{
	class simple_mesh_description {
	public:
		simple_mesh_description():
			id(-1), is_copy(false),
			ompss2_bc_rows_0(nullptr), ompss2_bc_rows_1(nullptr),
			bc_rows_0_size(0), bc_rows_1_size(0),
			ompss2_ids_to_rows(nullptr), ids_to_rows_size(0)
		{}

		simple_mesh_description(const simple_mesh_description &in):
			id(in.id), is_copy(true),
			ompss2_bc_rows_0(in.ompss2_bc_rows_0),
			ompss2_bc_rows_1(in.ompss2_bc_rows_1),
			bc_rows_0_size(in.bc_rows_0_size),
			bc_rows_1_size(in.bc_rows_1_size),
			ompss2_ids_to_rows(in.ompss2_ids_to_rows),
			ids_to_rows_size(in.ids_to_rows_size)
		{}

		~simple_mesh_description()
		{
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

		int map_id_to_row(int id) const
		{
			return find_row_for_id(id,
			                       ompss2_ids_to_rows,
			                       &ompss2_ids_to_rows[ids_to_rows_size]);
		}

		void write(std::ostream &stream = std::cout) const
		{
			stream << "mesh: " << id << "\n";
			stream << "global_box: " <<  global_box << "\n";
			stream << "local_box: " << local_box << "\n";
			stream << "extended_box: " << extended_box << "\n";

			dbvprint_vector("bc_rows_0", bc_rows_0_size, ompss2_bc_rows_0, stream);
			stream << std::endl;
			dbvprint_vector("bc_rows_1", bc_rows_1_size, ompss2_bc_rows_1, stream);
			stream << std::endl;
			dbvprint_vector("ids_to_rows", ids_to_rows_size, ompss2_ids_to_rows, stream);
			stream << std::endl;
		}

		int id;
		bool is_copy;

		int *ompss2_bc_rows_0;
		int *ompss2_bc_rows_1;
		size_t bc_rows_0_size, bc_rows_1_size;

		std::pair<int,int> *ompss2_ids_to_rows;
		size_t ids_to_rows_size;

		Box global_box;
		Box local_box;
		Box extended_box;
	}; //class simple_mesh_description


	void init_mesh_all(simple_mesh_description *mesh_array,
	                   double *ptr,
	                   const Box *global_box_in,
	                   const Box *local_boxes_array,    // Global boxes
	                   const Box *local_node_box_array, // Boxes resized
	                   size_t numboxes)
	{

		for (size_t id = 0; id < numboxes; ++id) {

			simple_mesh_description *mesh = &mesh_array[id];

			#pragma oss task				\
				out(mesh[0])				\
				in(ptr[id * 1000; 1000])			\
				in(global_box_in[0])			\
				in(local_boxes_array[0; numboxes])	\
				in(local_node_box_array[0; numboxes])
			{
				mesh->id = id;
				mesh->global_box = *global_box_in;
				mesh->local_box = local_boxes_array[id];
				// Internal class copies of the boxes
				mesh->extended_box = local_node_box_array[id];
				const int max_node_x = mesh->global_box[0][1] + 1;


				const int max_node_y = mesh->global_box[1][1] + 1;
				const int max_node_z = mesh->global_box[2][1] + 1;

				// Local set and map. Will be copied later
				std::set<int> bc_rows_0, bc_rows_1;
				std::map<int, int> map_ids_to_rows =
					create_map_id_to_row(max_node_x, max_node_y, max_node_z,
					                     local_node_box_array, id, numboxes);

				//As described in analytic_soln.hpp,
				//we will impose a 0 boundary-condition on faces x=0, y=0, z=0, y=1, z=1
				//we will impose a 1 boundary-condition on face x=1

				const int X = 0;
				const int Y = 1;
				const int Z = 2;

				const int x1 = max_node_x - 1;
				const int y1 = max_node_y - 1;
				const int z1 = max_node_z - 1;

				//if we're on the x=0 face:
				if (mesh->global_box[X][0] == mesh->local_box[X][0]) {
					int miny = mesh->extended_box[Y][0];
					int minz = mesh->extended_box[Z][0];
					int maxy = mesh->extended_box[Y][1];
					int maxz = mesh->extended_box[Z][1];
					//expand y and z dimensions to include ghost layer
					if (mesh->extended_box[Y][0] > 0) --miny;
					if (mesh->extended_box[Z][0] > 0) --minz;
					if (mesh->extended_box[Y][1] < max_node_y) ++maxy;
					if (mesh->extended_box[Z][1] < max_node_z) ++maxz;

					for (int iz = minz; iz < maxz; ++iz) {
						for(int iy = miny; iy < maxy; ++iy) {
							const int nodeID = get_id(max_node_x, max_node_y, max_node_z,
							                          0, iy, iz);

							dbv2printf("x = 0 BC, node %d, (%d, %d, %d)\n",
							           nodeID, 0, iy, iz);

							bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
						}
					}
				}

				//if we're on the y=0 face:
				if (mesh->global_box[Y][0] == mesh->local_box[Y][0]) {
					int minx = mesh->extended_box[X][0];
					int minz = mesh->extended_box[Z][0];
					int maxx = mesh->extended_box[X][1];
					int maxz = mesh->extended_box[Z][1];
					//expand x and z dimensions to include ghost layer
					if (mesh->extended_box[X][0] > 0) --minx;
					if (mesh->extended_box[Z][0] > 0) --minz;
					if (mesh->extended_box[X][1] < max_node_x) ++maxx;
					if (mesh->extended_box[Z][1] < max_node_z) ++maxz;

					for(int iz=minz; iz<maxz; ++iz) {
						for(int ix=minx; ix<maxx; ++ix) {
							int nodeID = get_id(max_node_x, max_node_y, max_node_z,
							                    ix, 0, iz);
							dbv2printf("y = 0 BC, node %d, (%d, %d, %d)\n",
							           nodeID, ix, 0, iz);

							int row = find_row_for_id(nodeID, map_ids_to_rows);

							if (row < 0)
								fprintf(stderr,
								        "on the y==0 face (ix=%d, iz=%d) \n"
								        "ERROR: found negative row (%d) for nodeID=%d\n",
								        ix, iz, row, nodeID);

							bc_rows_0.insert(row);
						}
					}
				}

				//if we're on the z=0 face:
				if (mesh->global_box[Z][0] == mesh->local_box[Z][0]) {
					int minx = mesh->extended_box[X][0];
					int miny = mesh->extended_box[Y][0];
					int maxx = mesh->extended_box[X][1];
					int maxy = mesh->extended_box[Y][1];
					//expand x and y dimensions to include ghost layer
					if (mesh->extended_box[X][0] > 0) --minx;
					if (mesh->extended_box[Y][0] > 0) --miny;
					if (mesh->extended_box[X][1] < max_node_x) ++maxx;
					if (mesh->extended_box[Y][1] < max_node_y) ++maxy;

					for(int iy = miny; iy < maxy; ++iy) {
						for(int ix = minx; ix < maxx; ++ix) {
							int nodeID = get_id(max_node_x, max_node_y, max_node_z,
							                    ix, iy, 0);
							dbv2printf("z = 0 BC, node %d, (%d, %d, %d)\n",
							           nodeID, ix, iy, 0);

							bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
						}
					}
				}

				//if we're on the x=1 face:
				if (mesh->global_box[X][1] == mesh->local_box[X][1]) {
					int minz = mesh->extended_box[Z][0];
					int miny = mesh->extended_box[Y][0];
					int maxz = mesh->extended_box[Z][1];
					int maxy = mesh->extended_box[Y][1];
					//expand z and y dimensions to include ghost layer
					if (mesh->extended_box[Z][0] > 0) --minz;
					if (mesh->extended_box[Y][0] > 0) --miny;
					if (mesh->extended_box[Z][1] < max_node_z) ++maxz;
					if (mesh->extended_box[Y][1] < max_node_y) ++maxy;

					for(int iy=miny; iy<maxy; ++iy) {
						for(int iz=minz; iz<maxz; ++iz) {
							int nodeID = get_id(max_node_x, max_node_y, max_node_z,
							                    x1, iy, iz);
							int row = find_row_for_id(nodeID, map_ids_to_rows);

							dbv2printf("x = 1 BC, node %d, row %d, (%d, %d, %d)\n",
							           nodeID, row, x1, iy, iz);

							bc_rows_1.insert(row);
						}
					}
				}

				//if we're on the y=1 face:
				if (mesh->global_box[Y][1] == mesh->local_box[Y][1]) {
					int minz = mesh->extended_box[Z][0];
					int minx = mesh->extended_box[X][0];
					int maxz = mesh->extended_box[Z][1];
					int maxx = mesh->extended_box[X][1];
					//expand z and x dimensions to include ghost layer
					if (mesh->extended_box[Z][0] > 0) --minz;
					if (mesh->extended_box[X][0] > 0) --minx;
					if (mesh->extended_box[Z][1] < max_node_z) ++maxz;
					if (mesh->extended_box[X][1] < max_node_x) ++maxx;

					for(int ix=minx; ix<maxx; ++ix) {
						for(int iz=minz; iz<maxz; ++iz) {
							int nodeID = get_id(max_node_x, max_node_y, max_node_z,
							                    ix, y1, iz);

							dbv2printf("y = 1 BC, node %d, (%d, %d, %d)\n",
							           nodeID, ix, y1, iz);

							bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
						}
					}
				}

				//if we're on the z=1 face:
				if (mesh->global_box[Z][1] == mesh->local_box[Z][1]) {
					int miny = mesh->extended_box[Y][0];
					int minx = mesh->extended_box[X][0];
					int maxy = mesh->extended_box[Y][1];
					int maxx = mesh->extended_box[X][1];
					//expand x and y dimensions to include ghost layer
					if (mesh->extended_box[Y][0] > 0) --miny;
					if (mesh->extended_box[X][0] > 0) --minx;
					if (mesh->extended_box[Y][1] < max_node_y) ++maxy;
					if (mesh->extended_box[X][1] < max_node_x) ++maxx;

					for(int ix = minx; ix < maxx; ++ix) {
						for(int iy = miny; iy < maxy; ++iy) {
							int nodeID =
								get_id(max_node_x, max_node_y, max_node_z,
								       ix, iy, z1);

							dbv2printf("z = 1 BC, node %d, (%d, %d, %d)\n",
							           nodeID, ix, iy, z1);

							bc_rows_0.insert(find_row_for_id(nodeID, map_ids_to_rows));
						}
					}
				}

				// Copy the sets to ompss supported ones.
				// TODO: this is a workaround due to the dmalloc issue.
				// I am using an lmalloc here that is never released to avoid double copy and taskwaits
				mesh->ompss2_bc_rows_0 = stl_to_local<int,
				                                      std::set<int>>(mesh->bc_rows_0_size, bc_rows_0);
				mesh->ompss2_bc_rows_1 = stl_to_local<int,
				                                      std::set<int>>(mesh->bc_rows_1_size, bc_rows_1);

				mesh->ompss2_ids_to_rows = stl_to_local<std::pair<int, int>,
				                                        std::map<int,int>>(mesh->ids_to_rows_size, map_ids_to_rows);

			}
		}


	}

}//namespace miniFE

#endif
