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

#ifndef _parameters_hpp_
#define _parameters_hpp_

#include <string>
#include <sstream>
#include <fstream>

#include <param_utils.hpp>

namespace miniFE {

	class Parameters {
	private:
		/**
		 * Concatenate command-line arguments into a single string.
		 *
		 * Note: this function is purely serial. If argc and argv have different
		 * values on different MPI processes, then you need to resolve that by
		 * broadcasting arg_string's contents.
		 */

		/**
		 * Read the contents of a text-file into a single string.
		 *
		 * Note: this function is purely serial. If you want file_contents on multiple
		 * MPI processes, you need to broadcast it (or call this function on each
		 * MPI process...).
		 */
		static std::string read_args_into_string(int argc, char** argv)
		{
			std::string arg_string = argv[0];
			for (int i = 1; i < argc; ++i)
				arg_string += " " + std::string(argv[i]);

			return arg_string;
		}


		/**
		 * Parse a named parameter value from input 'arg_string'.
		 *
		 * Search 'arg_string' for an occurrence of param_name and attempt to parse
		 * a value into the return-type. If param_name is not found, then default_value
		 * is returned.
		 *
		 * Example:
		 * arg_string = "foo = 3.14159";
		 * float foo = parse_parameter<float>(arg_string, "foo", -999.9);
		 * //foo should now contain the value 3.14159; if 'foo' was not found in
		 * //arg_string, then -999.9 would have been returned.
		 *
		 * Other legal name-value separators are ':' and ' '. Extra spaces are also ok,
		 * e.g. "foo : 3.114159".
		 *
		 * Note that if a YAML file is read into a string, that would be a valid input
		 * string for this function.
		 */
		static std::string read_file_into_string(const std::string& filename)
		{
			std::string file_contents;
			std::ifstream ifs(filename.c_str());
			char line[256];
			while(!ifs.eof()) {
				ifs.getline(line, 256);
				file_contents += " " + std::string(line);
			}
			return file_contents;
		}

	public:
		Parameters(int argc, char** argv)
			: nx(5), ny(nx), nz(nx), numboxes(1),
			  use_locking(0),
			  load_imbalance(0), name(), elem_group_size(1),
			  use_elem_mat_fields(1), verify_solution(0),
			  device(0),num_devices(2),skip_device(9999),numa(1)
		{
			std::string argstring = read_args_into_string(argc, argv);

			std::string garbage("garbage");
			std::string filename =
				Mantevo::parse_parameter<std::string>(argstring, "input_file", garbage);

			if (filename != garbage)
				read_file_into_string(filename);

			nx = Mantevo::parse_parameter<int>(argstring, "nx", 10);
			ny = Mantevo::parse_parameter<int>(argstring, "ny", nx);
			nz = Mantevo::parse_parameter<int>(argstring, "nz", ny);
			load_imbalance = Mantevo::parse_parameter<float>(argstring, "load_imbalance", 0);
			numboxes = Mantevo::parse_parameter<int>(argstring, "numboxes", 1);

			use_locking = Mantevo::parse_parameter<int>(argstring, "use_locking", 0);
			name = Mantevo::parse_parameter<std::string>(argstring, "name","");
			elem_group_size = Mantevo::parse_parameter<int>(argstring, "elem_group_size", 1);
			use_elem_mat_fields = Mantevo::parse_parameter<int>(argstring, "use_elem_mat_fields", 1);
			verify_solution = Mantevo::parse_parameter<int>(argstring, "verify_solution", 0);
			device = Mantevo::parse_parameter<int>(argstring, "device", 0);
			num_devices = Mantevo::parse_parameter<int>(argstring, "num_devices", 2);
			skip_device = Mantevo::parse_parameter<int>(argstring, "skip_device", 9999);
			numa = Mantevo::parse_parameter<int>(argstring, "numa", 1);

		}

		int nx, ny, nz;
		int numboxes;
		int use_locking;
		float load_imbalance;
		std::string name;
		int elem_group_size;
		int use_elem_mat_fields;
		int verify_solution;
		int device;
		int num_devices;
		int skip_device;
		int numa;
	}; //struct Parameters

}//namespace miniFE

#endif
