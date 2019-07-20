/*
 * Copyright (C) 2019  Jimmy Aguilar Mena
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OMPSS_UTILS_HPP
#define OMPSS_UTILS_HPP

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <set>
#include <map>
#include <cstring>
#include <fstream>      // std::ofstream
#include <functional>

// nanos conditional macros here

#ifdef NANOS6
#include "nanos6.h"
#define dmalloc(size) nanos6_dmalloc(size, nanos6_equpart_distribution, 0, NULL)
#define lmalloc(size) nanos6_lmalloc(size)
#define dfree(ptr, size) nanos6_dfree(ptr, size)
#define lfree(ptr, size) nanos6_lfree(ptr, size)
#else
#define nanos6_get_cluster_node_id() 0
#define nanos6_get_cluster_nodes() 1

#define dmalloc(size) malloc(size)
#define lmalloc(size) malloc(size)
#define dfree(ptr, size) free(ptr)
#define lfree(ptr, size) free(ptr)
#endif


// General use macros Here

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

inline std::string prefix(std::string caller, std::string file, int line, std::string varname) {
	std::string ret = "VERB_" + caller + "_"  + std::to_string(line) + "_" + varname + "_";
	return ret;
}

#define REGISTER_ELAPSED_TIME(time_inc, time_total)			\
	{								\
		time_inc = mytimer() - time_inc;			\
		time_total += time_inc;					\
	}

// Debug conditional macros here.

#ifndef NDEBUG
#define dbprintf(...) fprintf(stderr, __VA_ARGS__)
#define dbprint_vector(...) print_vector(__VA_ARGS__)
#else
#define dbprintf(...)
#define dbprint_vector(...)
#endif

#ifdef VERBOSE
#define dbvprintf(ARG, ...) fprintf(stderr, "Node: %d: " ARG, nanos6_get_cluster_node_id(), ##__VA_ARGS__)
#define dbvwrite(var) write(var, prefix(__func__, __FILE__, __LINE__, #var) )
#define dbvprint_vector(...) print_vector(__VA_ARGS__)
#define rrd_malloc(size) _rrd_malloc(size, __FILE__ ":"  STR(__LINE__))
#define rrl_malloc(size) _rrl_malloc(size, __FILE__ ":"  STR(__LINE__))
#define rrd_free(var, size) _rrd_free(var, size, __FILE__ ":"  STR(__LINE__) "(" #var ")")
#define rrl_free(var, size) _rrl_free(var, size, __FILE__ ":"  STR(__LINE__) "(" #var ")")
#else
#define dbvprintf(...)
#define dbvwrite(var)
#define dbvprint_vector(...)
#define rrd_malloc(size) _rrd_malloc(size)
#define rrl_malloc(size) _rrl_malloc(size)
#define rrd_free(var, size) _rrd_free(var, size)
#define rrl_free(var, size) _rrl_free(var, size)
#endif


#if (VERBOSE > 1)
#define dbv2printf(...) fprintf(stderr, __VA_ARGS__)
#define dbv2print_vector(...) print_vector(__VA_ARGS__)
#define dbv2write(var) _write(var, "VERB_"__func__ "_" __FILE__ ":"  STR(__LINE__) "_" #var "_")
#else
#define dbv2printf(...)
#define dbv2print_vector(...)
#define dbv2write(var)
#endif

// General Purpose functions

template <typename T>
std::ostream &array_to_stream(const T *in, size_t size,
                              std::string varname = "",
                              const std::string sep = "\n",
                              std::ostream &stream = std::cout)
{
	stream << varname << sep;
	for (size_t i = 0; i < size; ++i)
		stream << in[i] << sep;
	return stream;
}


template <typename T>
void print_vector(std::string vname, size_t size, const  T *vect, std::ostream &stream = std::cout)
{
	#ifdef VERBOSE
	stream << "Node: " << nanos6_get_cluster_node_id() << ": ";
	#endif
	stream << vname
	       << "(" << vect << ":" << size * sizeof(T) << ")"
	       << "[" << size << "]={";
	for (size_t i = 0; i < size; ++i) {
		if (i > 0)
			stream << "; ";

		stream << vect[i];
	}
	stream << " }";
}

template <typename _Key, typename _Tp>
std::ostream &operator<< (std::ostream& os, const std::pair<_Key, _Tp> &in)
{
	os << "<" << in.first << ";" << in.second << ">";
	return os;
}

template <typename T>
inline void write_all(std::string &filename, const T *in_array, size_t numboxes)
{
	for (size_t id = 0; id < numboxes; ++id)
		write_task(filename, in_array[id], id);
}

template <typename T>
inline void reduce_sum_task(T *vout, const T *vin, size_t size)
{
	#pragma oss task in(vin[0; size]) out(vout[0; 1])
	{
		*vout = 0;
		for (size_t i = 0; i < size; ++i)
			*vout += vin[i];

		dbvprint_vector("Reducing: ", size, vin);
		#ifdef VERBOSE
		std::cout << " = " << *vout << std::endl;
		#endif
	}
}


inline void ompss_memcpy_task(void *pout, const void *pin, size_t size)
{
	char *tin = (char *) pin;
	char *tout = (char *) pout;

	assert(tin != nullptr);
	assert(tout != nullptr);

	#pragma oss task in(tin[0; size]) out(tout[0; size])
	{
		dbv2printf("Copy %ld bytes from %p -> %p\n", size, pin, pout);
		memcpy(tout, tin, size);
	}
}

inline void *_rrd_malloc(size_t size, const char info[] = "")
{
	void *ret = dmalloc(size);
	dbv2printf("%s %s(%lu) = [%p -> %p]\n", info, __func__, size,  ret, (char*)ret + size);

	assert(size == 0 || ret != NULL);
	return ret;
}

inline void _rrd_free(void *in, size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%p, %lu)\n", info, __func__, in, size);
	dfree(in, size);
}

inline void *_rrl_malloc(size_t size, const char info[] = "")
{
	void *ret = lmalloc(size);
	dbv2printf("%s %s(%lu) = [%p -> %p]\n", info, __func__, size, ret, (char*)ret + size);
	assert(size == 0 || ret != NULL);

	return ret;
}

inline void _rrl_free(void *in, size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%p, %lu)\n", info, __func__, in, size);
	lfree(in, size);
}

template<typename T, typename Container>
T *stl_to_local(size_t &copied, const Container &vin)
{
	copied = vin.size();

	T *ret = (T *) rrl_malloc(copied * sizeof(T));

	dbvprintf("STL copy %ld elements -> %p\n", copied, (void *) ret);

	// Copy from container to local memory, this is initialization any way.
	size_t i = 0;
	for (const T &a : vin)
		ret[i++] = a;

	return ret;
}

template<typename T>
void write(const T *A, std::string prefix)
{
	std::string filename = prefix + std::to_string(A->id) + ".verb";
	std::ofstream stream(filename);
	assert(stream.is_open());
	A->write(stream);
	stream.close();
}



template<typename T, typename R>
void get_offsets(R &count, R *indices, R *offsets, size_t numboxes,
                 const T *array, std::function<R(const T &)> f)
{
	count = 0;
	for (size_t id = 0; id < numboxes; ++id) {
		const R tmp = f(array[id]);

		indices[id] = tmp;
		offsets[id] = count;

		count += tmp;
	}
}

#endif // OMPSS_UTILS_HPP
