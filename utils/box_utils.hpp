#ifndef _box_utils_hpp_
#define _box_utils_hpp_

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

#include <vector>
#include <set>
#include <map>

#include <cassert>

#include <TypeTraits.hpp>
#include <Box.hpp>

namespace miniFE
{

	inline void get_int_coords(int ID,
	                           int nx, int ny, int nz,
	                           int &x, int &y, int &z)
	{
		z = ID / (nx * ny);
		y = (ID % (nx * ny)) / nx;
		x = ID % nx;
	}

	inline void get_coords(int ID, int nx, int ny, int nz,
	                       double &x, double &y, double &z)
	{
		const int xdiv = nx > 1 ? nx - 1 : 1;
		const int ydiv = ny > 1 ? ny - 1 : 1;
		const int zdiv = nz > 1 ? nz - 1 : 1;

		//This code assumes that ID is 0-based.
		//
		//compute coordinates that lie on (or in) the unit cube.
		//that's why we're dividing by nz,ny,nx:
		z = (1.0 * (ID / (nx * ny))) / zdiv;
		y = 1.0 * ((ID % (nx * ny)) / nx) / ydiv;
		x = 1.0 * (ID % nx) / xdiv;
	}


	inline void print_box(const char* name, const Box& box,
	                      const char* name2, const Box& box2)
	{
		std::cout << name
		          << " ("<<box[0][0]<<","<<box[0][1]<<") "
		          << " ("<<box[1][0]<<","<<box[1][1]<<") "
		          << " ("<<box[2][0]<<","<<box[2][1]<<") "
		          << name2
		          << " ("<<box2[0][0]<<","<<box2[0][1]<<") "
		          << " ("<<box2[1][0]<<","<<box2[1][1]<<") "
		          << " ("<<box2[2][0]<<","<<box2[2][1]<<") "
		          << std::endl;
	}


	inline std::map<int,int> create_map_id_to_row(
		int global_nx, int global_ny, const int global_nz,
		const Box *local_node_box_array, size_t id, size_t numboxes)
	{

		assert(id < numboxes);

		const Box &box = local_node_box_array[id];

		std::vector<int> all_ids = box.get_ids(global_nx, global_ny, global_nz, false);

		// This substituted the all_gather
		int global_offsets[numboxes];
		for (size_t i = 0; i < numboxes; ++i)
			global_offsets[i] = local_node_box_array[i].get_num_ids();

		int offset = 0;
		for(size_t i = 0; i < numboxes; ++i) {
			const int tmp = global_offsets[i];
			global_offsets[i] = offset;
			offset += tmp;
		}

		const int my_first_row = global_offsets[id];

		std::map<int,int> id_to_row;

		if (all_ids.size() > 0)
			id_to_row.insert(std::make_pair(all_ids[0], my_first_row));

		for (size_t i = 1; i < all_ids.size(); ++i)
			if (all_ids[i] != all_ids[i - 1] + 1)
				id_to_row[all_ids[i]] = my_first_row + i;


		for (size_t i = 0; i < numboxes; ++i) {
			if (i == id)
				continue;

			const Box &box_i = local_node_box_array[i];

			if (!box.is_neighbor(box_i))
				continue;

			all_ids = box_i.get_ids(global_nx, global_ny, global_nz, false);

			const int first_row_i = global_offsets[i];
			if (all_ids.size() > 0)
				id_to_row[all_ids[0]] = first_row_i;

			for (size_t j = 1; j < all_ids.size(); ++j)
				if (all_ids[j] != all_ids[j - 1] + 1)
					id_to_row[all_ids[j]] = first_row_i + j;
		}

		return id_to_row;
	}

}//namespace miniFE

#endif

