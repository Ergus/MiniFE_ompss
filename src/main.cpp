
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

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>

#include <miniFE_version.h>

#include <outstream.hpp>

#include <Box.hpp>
#include <BoxPartition.hpp>
#include <box_utils.hpp>

#include <CSRMatrix.hpp>
#include <Parameters.hpp>
#include <utils.hpp>
#include <driver.hpp>
#include <YAML_Doc.hpp>

#include <miniFE_info.hpp>

//The following macros should be specified as compile-macros in the
//makefile. They are defaulted here just in case...
#ifndef MINIFE_SCALAR
#define MINIFE_SCALAR double
#endif
#ifndef MINIFE_LOCAL_ORDINAL
#define MINIFE_LOCAL_ORDINAL int
#endif
#ifndef MINIFE_GLOBAL_ORDINAL
#define MINIFE_GLOBAL_ORDINAL int
#endif

// ************************************************************************

//
//We will create a 'box' of size nx X ny X nz, partition it among processors,
//then call miniFE::driver which will use the partitioned box as the domain
//from which to assemble finite-element matrices into a global matrix and
//vector, then solve the linear-system using Conjugate Gradients.
//


// numboxes will be like the numtasks
// myproc will be like an ordinal taskid
int main(int argc, char** argv)
{
	miniFE::Parameters params(argc, argv);

	miniFE::timer_type start_time = miniFE::mytimer();

	dbprintf("Num boxes: %d\n", params.numboxes);

	Box global_box(params.nx, params.ny, params.nz);

	//This arrays is in local memory
	Box *local_boxes = new Box[params.numboxes];

	// This one is recursive, play with weak here (but it is a cheap)
	box_partition(0, params.numboxes, 2, global_box, local_boxes);

	dbarray_to_stream(local_boxes, params.numboxes);

	std::ostringstream osstr;
	osstr << "miniFE." << params.nx << "x" << params.ny << "x" << params.nz;
	osstr << ".";
	if (params.name != "")
	 	osstr << params.name << ".";

	YAML_Doc doc("miniFE", MINIFE_VERSION, ".", osstr.str());
	doc.add_params_to_yaml(params);
	doc.add_configuration_to_yaml(params.numboxes);
	doc.add_timestring_to_yaml();

	// //Most of the program is performed in the 'driver' function, which is
	// //templated on < Scalar, LocalOrdinal, GlobalOrdinal >.
	// //To run miniFE with float instead of double, or 'long long' instead of int,
	// //etc., change these template-parameters by changing the macro definitions in
	// //the makefile or on the make command-line.

	// This can be a weak task
	miniFE::driver(global_box, local_boxes, params.numboxes, params, doc);

	// TODO: taskwait here
	miniFE::timer_type total_time = miniFE::mytimer() - start_time;

	doc.add("Total Program Time", total_time);
	doc.generateYAML();

	#pragma oss taskwait
	delete [] local_boxes;

	return 0;
}

