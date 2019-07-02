#ifndef _verify_solution_hpp_
#define _verify_solution_hpp_

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

#include <sstream>
#include <stdexcept>
#include <map>
#include <algorithm>

#include <simple_mesh_description.hpp>
#include <analytic_soln.hpp>
#include <box_utils.hpp>
#include <utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace miniFE
{

	struct err_info {
		double err;
		double computed;
		double analytic;
		double coords[3];

		err_info() : err(0.0)
		{}
	};

	inline std::ostream& operator <<(std::ostream &stream, const err_info &in) {
		stream << "max absolute error is " << in.err << "\n"
		       << "   at position ("
		       << in.coords[0] << ","
		       << in.coords[1] << ","
		       << in.coords[2] << "), \n"
		       << " computed solution: " << in.computed
		       << " analytic solution: " << in.analytic
		       << std::endl;
		return stream;
	}

	const inline err_info &get_max_error(const err_info *in, size_t numboxes) {
		double max_err = in[0].err;
		size_t idx = 0;
		for (size_t i = 1; i < numboxes; ++i)
			if (in[i].err > max_err)
				idx = i;
		return in[idx];
	}

	inline int verify_solution(const Box &global_box,
	                           const Box *local_node_box_array,
	                           const Vector *x_array, int numboxes,
	                           double tolerance,
	                           bool verify_whole_domain = false)
	{

		const int global_nodes_x = global_box[0][1]+1;
		const int global_nodes_y = global_box[1][1]+1;
		const int global_nodes_z = global_box[2][1]+1;
		err_info max_error[numboxes];

		for (int id = 0; id < numboxes; ++id) {
			const Box &box = local_node_box_array[id];

			std::vector<int> rows;
			std::vector<double> row_coords;

			int roffset = 0;
			for(int iz=box[2][0]; iz<box[2][1]; ++iz) {
				for(int iy = box[1][0]; iy < box[1][1]; ++iy) {
					for(int ix = box[0][0]; ix < box[0][1]; ++ix) {
						int row_id =
							get_id(global_nodes_x, global_nodes_y, global_nodes_z,
							       ix, iy, iz);
						double x, y, z;
						get_coords(row_id, global_nodes_x, global_nodes_y, global_nodes_z, x, y, z);

						if (verify_whole_domain || (std::abs(x - 0.5) < 0.05 &&
						                            std::abs(y - 0.5) < 0.05 &&
						                            std::abs(z - 0.5) < 0.05)) {

							rows.push_back(roffset);
							row_coords.push_back(x);
							row_coords.push_back(y);
							row_coords.push_back(z);

						}
						++roffset;
					}
				}
			}

			for(size_t i = 0; i < rows.size(); ++i) {
				double computed_soln = x_array[id].coefs[rows[i]];
				double x = row_coords[i * 3];
				double y = row_coords[i * 3 + 1];
				double z = row_coords[i * 3 + 2];
				double analytic_soln = 0.0;
				//set exact boundary-conditions:
				//x==1 is first, we want soln to be 1 even around the edges
				//of the x==1 plane where y and/or z may be 0 or 1...
				if (x == 1.0)
					analytic_soln = 1;
				else if (x == 0.0 || y == 0.0 || z == 0.0)
					analytic_soln = 0;
				else if (y == 1.0 || z == 1.0)
					analytic_soln = 0;
				else
					analytic_soln = soln(x, y, z, 300, 300);

				#ifdef MINIFE_DEBUG_VERBOSE
				std::cout << "(" << x << "," << y << "," << z << ")"
				          << " row " << rows[i]
				          << ": computed: " << computed_soln
				          << ",  analytic: " << analytic_soln
				          << std::endl;
				#endif
				double err = std::abs(analytic_soln - computed_soln);
				if (err > max_error[id].err) {
					max_error[id].err = err;
					max_error[id].computed = computed_soln;
					max_error[id].analytic = analytic_soln;
					max_error[id].coords[0] = x;
					max_error[id].coords[1] = y;
					max_error[id].coords[2] = z;
				}
			}
		}

		const err_info &global_max_error = get_max_error(max_error, numboxes);

		std::cout << global_max_error << std::endl;

		return (global_max_error.err > tolerance);
	}

}//namespace miniFE

#endif

