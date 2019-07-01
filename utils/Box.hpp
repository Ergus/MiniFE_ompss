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

#ifndef _Box_hpp_
#define _Box_hpp_

#include <iostream>
#include <cstring>
#include <vector>

#include "ompss_utils.hpp"
/**
  * a 'Box' is 3 pairs of ints, where each pair specifies a lower
  * and upper bound for one of the 3 spatial dimensions.
  *
  * This struct stores the 3 pairs as a simple array of 6 ints,
  * but defines the bracket operator so that it can be referenced
  * using 2-dimensional array notation like this:
  * int xmin = box[0][0]; int xmax = box[0][1];
  * int ymin = box[1][0]; int ymax = box[1][1];
  * int zmin = box[2][0]; int zmax = box[2][1];
  */

inline int get_id(int nx, int ny, int nz, int x, int y, int z)
{
	if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz)
		return -1;

	return x + nx * y + nx * ny * z;
}

class Box {
public:
	int ranges[6];

	Box()
	{
		memset(ranges, 0, 6 * sizeof(int));
	}

	Box(const Box &in)
	{
		memcpy(ranges, in.ranges, 6 * sizeof(int));
	}

	Box(int nx, int ny, int nz)
	{
		ranges[0] = 0;
		ranges[1] = nx;
		ranges[2] = 0;
		ranges[3] = ny;
		ranges[4] = 0;
		ranges[5] = nz;
	}

	// Operator []
	int *operator[](const int xyz)
	{
		return &ranges[xyz * 2];
	}

	const int *operator[](const int xyz) const
	{
		return &ranges[xyz * 2];
	}

	// Operator =
	void operator= (const Box &in)
	{
		memcpy(ranges, in.ranges, 6 * sizeof(int));
	}

	// Arrays allocations (boxes can go in local memory)
	static void* operator new[](std::size_t sz)
	{
		void * const tmp = rrl_malloc(sz);
		dbprintf("Calling: %s, size: %lu\n", __PRETTY_FUNCTION__, sz);
		return tmp;
	}

	static void operator delete[](void* ptr, std::size_t sz)
	{
		printf("Calling: %s, address %p size: %lu\n", __PRETTY_FUNCTION__, ptr, sz);
		return rrl_free(ptr, sz);
	}


	bool is_neighbor(const Box &otherbox) const
	{
		//neighbors in the x dimension if:
		bool x_neighbor = ((*this)[0][1] == otherbox[0][0]) || ((*this)[0][0] == otherbox[0][1]) || // min matches max
			((*this)[0][0] == otherbox[0][0]) ||
			((*this)[0][1] == otherbox[0][1]) || // mins or maxs match
			((*this)[0][0] >  otherbox[0][0]  &&  (*this)[0][1] <  otherbox[0][1]) || // range contains other
			(otherbox[0][0] >  (*this)[0][0]  &&  otherbox[0][1] <  (*this)[0][1]) || // range contains other
			((*this)[0][0] >  otherbox[0][0]  &&  (*this)[0][0] <  otherbox[0][1]) || // min contained in rng
			(otherbox[0][0] >  (*this)[0][0]  &&  otherbox[0][0] <  (*this)[0][1]);   // min contained in rng
		if (!x_neighbor)
			x_neighbor = ((*this)[0][1] == otherbox[0][0]-1) || ((*this)[0][0] == otherbox[0][1]+1);

		bool y_neighbor = ((*this)[1][1] == otherbox[1][0]) || ((*this)[1][0] == otherbox[1][1]) || // min matches max
			((*this)[1][0] == otherbox[1][0]) || ((*this)[1][1] == otherbox[1][1]) || // mins or maxs match
			((*this)[1][0] >  otherbox[1][0]  &&  (*this)[1][1] <  otherbox[1][1]) || // range contains other
			(otherbox[1][0] >  (*this)[1][0]  &&  otherbox[1][1] <  (*this)[1][1]) || // range contains other
			((*this)[1][0] >  otherbox[1][0]  &&  (*this)[1][0] <  otherbox[1][1]) || // min contained in rng
			(otherbox[1][0] >  (*this)[1][0]  &&  otherbox[1][0] <  (*this)[1][1]);   // min contained in rng
		if (!y_neighbor)
			y_neighbor = ((*this)[1][1] == otherbox[1][0]-1) || ((*this)[1][0] == otherbox[1][1]+1);

		bool z_neighbor = ((*this)[2][1] == otherbox[2][0]) || ((*this)[2][0] == otherbox[2][1]) || // min matches max
			((*this)[2][0] == otherbox[2][0]) || ((*this)[2][1] == otherbox[2][1]) || // mins or maxs match
			((*this)[2][0] >  otherbox[2][0]  &&  (*this)[2][1] <  otherbox[2][1]) || // range contains other
			(otherbox[2][0] >  (*this)[2][0]  &&  otherbox[2][1] <  (*this)[2][1]) || // range contains other
			((*this)[2][0] >  otherbox[2][0]  &&  (*this)[2][0] <  otherbox[2][1]) || // min contained in rng
			(otherbox[2][0] >  (*this)[2][0]  &&  otherbox[2][0] <  (*this)[2][1]);   // min contained in rng
		if (!z_neighbor)
			z_neighbor = ((*this)[2][1] == otherbox[2][0]-1) || ((*this)[2][0] == otherbox[2][1]+1);

		return x_neighbor && y_neighbor && z_neighbor;
	}

	int get_num_ids() const
	{
		const int nx = (*this)[0][1] - (*this)[0][0];
		const int ny = (*this)[1][1] - (*this)[1][0];
		const int nz = (*this)[2][1] - (*this)[2][0];
		return nx * ny * nz;
	}

	std::vector<int> get_ids(int nx, int ny, int nz,
	                         bool include_ghost_layer = false) const
	{
		int minz = (*this)[2][0];
		int maxz = (*this)[2][1];
		int miny = (*this)[1][0];
		int maxy = (*this)[1][1];
		int minx = (*this)[0][0];
		int maxx = (*this)[0][1];

		if (include_ghost_layer) {
			if (minz > 0) minz--;
			if (miny > 0) miny--;
			if (minx > 0) minx--;
			if (maxz < nz) maxz++;
			if (maxy < ny) maxy++;
			if (maxx < nx) maxx++;
		}

		const size_t ids_size = ((maxz - minz) * (maxy - miny)) * (maxx - minx);
		std::vector<int> ids(ids_size);

		size_t idx = 0;
		for(int z = minz; z < maxz; ++z)
			for(int y = miny; y < maxy; ++y)
				for(int x = minx; x < maxx; ++x)
					ids[idx++] = get_id(nx, ny, nz, x, y, z);

		return ids;
	}


	void get_ghost_ids(int nx, int ny, int nz, std::vector<int> &ids) const
	{
		ids.clear();
		int minz,maxz,miny,maxy,minx,maxx;
		int orig_minz = minz = (*this)[2][0];
		int orig_maxz = maxz = (*this)[2][1];
		int orig_miny = miny = (*this)[1][0];
		int orig_maxy = maxy = (*this)[1][1];
		int orig_minx = minx = (*this)[0][0];
		int orig_maxx = maxx = (*this)[0][1];

		if (minz > 0) minz--;
		if (miny > 0) miny--;
		if (minx > 0) minx--;
		if (maxz < nz) maxz++;
		if (maxy < ny) maxy++;
		if (maxx < nx) maxx++;

		for(int z = minz; z < maxz; ++z) {
			for(int y = miny; y<maxy; ++y) {
				for(int x = minx; x < maxx; ++x) {
					const bool x_in_ghost_layer = (x < orig_minx) || (x >= orig_maxx);
					const bool y_in_ghost_layer = (y < orig_miny) || (y >= orig_maxy);
					const bool z_in_ghost_layer = (z < orig_minz) || (z >= orig_maxz);
					//we are in the ghost layer if any one of x,y,z are in the ghost layer
					if (!x_in_ghost_layer && !y_in_ghost_layer && !z_in_ghost_layer)
						continue;
					ids.push_back(get_id(nx, ny, nz, x, y, z));
				}
			}
		}
	}
};

inline std::ostream& operator <<(std::ostream &stream, const Box &in) {
	for (int i = 0; i < 6; ++i)
		stream << in.ranges[i] << "; ";
	stream << std::endl;
	return stream;
}

#endif

