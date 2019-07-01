#ifndef _utils_hpp_
#define _utils_hpp_

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

#include <cstdlib>
#include <cmath>
#include <climits>
#include <vector>
#include <map>
#include <iostream>


#include <Parameters.hpp>

#include <algorithm>


namespace miniFE
{

	inline void sort_if_needed(int *list, int list_len)
	{
		bool need_to_sort = false;
		for(int i = list_len - 1; i >= 1; --i) {
			if (list[i] < list[i-1]) {
				need_to_sort = true;
				break;
			}
		}

		if (need_to_sort)
			std::sort(list,list+list_len);
	}


	template<typename Scalar>
	Scalar percentage_difference(Scalar value, Scalar average)
	{
		//result will be the difference between value and average, represented as
		//a percentage of average.
		//Examples:
		//  if value=100 and average=50, result is 100%
		//  if value=500 and average=400, result is 25%

		//Note: if average is 0, result is undefined. We'll return -1.0;

		const Scalar result = std::abs(value - average);

		if (std::abs(average) > 1.e-5)
			return result * 100 / average ;

		return -1;
	}

	template<typename T>
	void get_global_min_max(T &global_n,
	                        T &min_n, int &min_box,
	                        T &max_n, int &max_box,
	                        const T *input,
	                        const size_t numboxes)
	{
		global_n = 0;
		min_n = INT_MAX;
		min_box = 0;
		max_n = 0;
		max_box = 0;

		for(size_t i = 0; i < numboxes; ++i) {
			global_n += input[i];
			if (input[i] < min_n) {
				min_n = input[i];
				min_box = i;
			}

			if (input[i] >= max_n) {
				max_n = input[i];
				max_box = i;
			}
		}
	}

	template<typename Tin, typename Tout>
	Tout compute_std_dev_as_percentage(const Tin *input, const size_t numboxes,
	                                const Tin avg_nrows)
	{

		if (numboxes <= 1)
			return 0.0;

		//turn all_nrows contents into deviations, add to sum-of-squares-of-deviations:
		Tout sum_sqr_dev = 0;
		for(size_t i = 0; i < numboxes; ++i) {
			const Tout tmp = input[i] - avg_nrows;
			sum_sqr_dev += (tmp * tmp);
		}
		Tout std_dev = std::sqrt(sum_sqr_dev / (numboxes - 1));

		return std_dev * 100.0 / avg_nrows;
	}

	template <typename maptype>
	int find_row_for_id(int id, const maptype &ids_to_rows)
	{
		auto iter = ids_to_rows.lower_bound(id);

		if (iter == ids_to_rows.end() || iter->first != id) {
			if (ids_to_rows.size() > 0)
				--iter;
			else {
				std::cerr << "ERROR, failed to map id to row."
				          << std::endl;
				return -99;
			}
		}

		if (iter->first == id)
			return iter->second;

		if (iter == ids_to_rows.begin() && iter->first > id) {
			std::cerr << "ERROR, id:" << id
			          << ", ids_to_rows.begin(): " << iter->first
			          << std::endl;
			return -99;
		}

		const int offset = id - iter->first;

		if (offset < 0) {
			std::cerr << "ERROR, negative offset in find_row_for_id for id="<< id
			          << std::endl;
			return -99;
		}

		return iter->second + offset;
	}

}//namespace miniFE

#endif

