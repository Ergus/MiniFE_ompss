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

#ifndef _imbalance_hpp_
#define _imbalance_hpp_

#include <cmath>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <box_utils.hpp>
#include <utils.hpp>
#include <YAML_Doc.hpp>

namespace miniFE
{

	const int X = 0;
	const int Y = 1;
	const int Z = 2;
	const int NONE = 3;

	const int LOWER = 0;
	const int UPPER = 1;

	inline void compute_imbalance(const Box &global_box, const Box *local_boxes,
	                              const size_t numboxes,
	                              float &largest_imbalance, float &std_dev,
	                              YAML_Doc &doc, const bool record_in_doc)
	{
		int min_nrows = 0, max_nrows = 0, global_nrows = 0;
		int min_proc = 0, max_proc = 0;

		int local_nrows[numboxes];
		for (size_t i = 0; i < numboxes; ++i) {
			// TODO: task here if not within a task all the function
			local_nrows[i] = local_boxes[i].get_num_ids();
		}

		get_global_min_max(global_nrows,
		                   min_nrows, min_proc,
		                   max_nrows, max_proc,
		                   local_nrows, numboxes);

		const float avg_nrows = (float) global_nrows / (float) numboxes;

		// largest_imbalance will be the difference between the min (or max)
		// rows-per-processor and avg_nrows, represented as a percentage:
		largest_imbalance = percentage_difference<float>(min_nrows, avg_nrows);

		const float tmp = percentage_difference<float>(max_nrows, avg_nrows);
		if (tmp > largest_imbalance)
			largest_imbalance = tmp;

		std_dev = compute_std_dev_as_percentage<int, float>(local_nrows, numboxes, avg_nrows);

		if (record_in_doc) {
			doc.add("Rows-per-proc Load Imbalance","");
			doc.get("Rows-per-proc Load Imbalance")->add("Largest (from avg, %)",
			                                             largest_imbalance);
			doc.get("Rows-per-proc Load Imbalance")->add("Std Dev (%)", std_dev);
		}
	}

	inline std::pair<int,int> decide_how_to_grow(const Box& global_box,
	                                             const Box& local_box)
	{
		std::pair<int,int> result(NONE,UPPER);

		if (local_box[Z][UPPER] < global_box[Z][UPPER]) {
			result.first = Z;
			result.second = UPPER;
			return result;
		}
		if (local_box[Z][LOWER] > global_box[Z][LOWER]) {
			result.first = Z;
			result.second = LOWER;
			return result;
		}
		if (local_box[Y][UPPER] < global_box[Y][UPPER]) {
			result.first = Y;
			result.second = UPPER;
			return result;
		}
		if (local_box[Y][LOWER] > global_box[Y][LOWER]) {
			result.first = Y;
			result.second = LOWER;
			return result;
		}
		if (local_box[X][UPPER] < global_box[X][UPPER]) {
			result.first = X;
			result.second = UPPER;
			return result;
		}
		if (local_box[X][LOWER] > global_box[X][LOWER]) {
			result.first = X;
			result.second = LOWER;
			return result;
		}
		return result;
	}

	inline std::pair<int,int> decide_how_to_shrink(const Box &global_box,
	                                               const Box &local_box)
	{
		std::pair<int,int> result(NONE,UPPER);

		if (local_box[Z][UPPER] < global_box[Z][UPPER] &&
		    local_box[Z][UPPER]-local_box[Z][LOWER] > 2) {
			result.first = Z;
			result.second = UPPER;
			return result;
		}
		if (local_box[Z][LOWER] > global_box[Z][LOWER] &&
		    local_box[Z][UPPER]-local_box[Z][LOWER] > 2) {
			result.first = Z;
			result.second = LOWER;
			return result;
		}
		if (local_box[Y][UPPER] < global_box[Y][UPPER] &&
		    local_box[Y][UPPER]-local_box[Y][LOWER] > 2) {
			result.first = Y;
			result.second = UPPER;
			return result;
		}
		if (local_box[Y][LOWER] > global_box[Y][LOWER] &&
		    local_box[Y][UPPER]-local_box[Y][LOWER] > 2) {
			result.first = Y;
			result.second = LOWER;
			return result;
		}
		if (local_box[X][UPPER] < global_box[X][UPPER] &&
		    local_box[X][UPPER]-local_box[X][LOWER] > 2) {
			result.first = X;
			result.second = UPPER;
			return result;
		}
		if (local_box[X][LOWER] > global_box[X][LOWER] &&
		    local_box[X][UPPER]-local_box[X][LOWER] > 2) {
			result.first = X;
			result.second = LOWER;
			return result;
		}
		return result;
	}

	inline void add_imbalance(const Box &global_box, Box *local_boxes,
	                          const size_t numboxes,
	                          const float imbalance, YAML_Doc &doc)
	{
		if (numboxes == 1)
			return;

		float cur_imbalance = 0, cur_std_dev = 0;
		compute_imbalance(global_box, local_boxes, numboxes,
		                  cur_imbalance, cur_std_dev,
		                  doc, false);

		while (cur_imbalance < imbalance) {

			int min_nrows = 0, max_nrows = 0, global_nrows = 0;
			int min_proc = 0, max_proc = 0;

			int local_nrows[numboxes];
			for (size_t i = 0; i < numboxes; ++i) {
				// TODO: Task here
				// in local_boxes[i]
				// out local_nrows[i]
				local_nrows[i] = local_boxes[i].get_num_ids();
			}

			get_global_min_max(global_nrows,
			                   min_nrows, min_proc,
			                   max_nrows, max_proc,
			                   local_nrows, numboxes);


			std::pair<int, int> grow(NONE, UPPER);
			Box &max_proc_box = local_boxes[max_proc];
			grow = decide_how_to_grow(global_box, max_proc_box);

			std::pair<int, int> shrink(NONE, UPPER);
			Box &min_proc_box = local_boxes[min_proc];
			shrink = decide_how_to_shrink(global_box, min_proc_box);

			if (grow.first == NONE && shrink.first == NONE)
				break;

			const int grow_axis_val = (grow.first != NONE) ?
				max_proc_box[grow.first][grow.second] : -1;

			const int shrink_axis_val = (shrink.first != NONE) ?
				min_proc_box[shrink.first][shrink.second] : -1;

			for (size_t i = 0; i < numboxes; ++i) {
				Box &local_box = local_boxes[i];
				if (grow.first != NONE) {
					if (local_box[grow.first][1] - local_box[grow.first][0] < 2) {
						if (grow.second != LOWER &&
						    local_box[grow.first][0] == grow_axis_val)
							++local_box[grow.first][0];

						else if (local_box[grow.first][1] == grow_axis_val)
							--local_box[grow.first][1];

					}
				}
				if (shrink.first != NONE) {
					if (local_box[shrink.first][1] - local_box[shrink.first][0] < 2) {
						if (shrink.second == LOWER &&
						    local_box[shrink.first][0] == shrink_axis_val)
							++local_box[shrink.first][0];

						else if (local_box[shrink.first][1] == shrink_axis_val)
							--local_box[shrink.first][1];
					}
				}
			}

			compute_imbalance(global_box, local_boxes,
			                  numboxes,
			                  cur_imbalance, cur_std_dev,
			                  doc, false);
		}
	}

}//namespace miniFE

#endif

