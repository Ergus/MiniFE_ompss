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

#ifdef NANOS6
#include "nanos6.h"
#else
#define nanos6_get_cluster_node_id() 0
#define nanos6_get_cluster_nodes() 1
#endif


// General use macros Here
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

#define REGISTER_ELAPSED_TIME(time_inc, time_total)			\
	{								\
		time_inc = mytimer() - time_inc;			\
		time_total += time_inc;					\
	}

// Debug conditional macros here.
#ifndef NDEBUG
#define dbprintf(...) fprintf(stderr, __VA_ARGS__)
#else
#define dbprintf(...)
#endif

#if !defined NDEBUG && defined VERBOSE
#define dbvprintf(ARG, ...) fprintf(stderr, "Node: %d: " ARG, nanos6_get_cluster_node_id(), ##__VA_ARGS__)
#if VERBOSE == 2
#define dbv2printf(...) fprintf(stderr, __VA_ARGS__)
#else
#define dbv2printf(...)
#endif
#else
#define dbvprintf(...)
#define dbv2printf(...)
#endif

#if !defined NDEBUG && defined VERBOSE
#define rrd_malloc(size) _rrd_malloc(size, __FILE__ ":"  STR(__LINE__))
#define rrl_malloc(size) _rrl_malloc(size, __FILE__ ":"  STR(__LINE__))
#define rrd_free(var, size) _rrd_free(var, size, __FILE__ ":"  STR(__LINE__) "(" #var ")")
#define rrl_free(var, size) _rrl_free(var, size, __FILE__ ":"  STR(__LINE__) "(" #var ")")
#else
#define rrd_malloc(size) _rrd_malloc(size)
#define rrl_malloc(size) _rrl_malloc(size)
#define rrd_free(var, size) _rrd_free(var, size)
#define rrl_free(var, size) _rrl_free(var, size)
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

#pragma oss task in(vin[0; size]) out(vout[0; 1])
inline void reduce_sum_task(double *vout, double *vin, size_t size)
{
	double tmp = 0.0;
	for (size_t i = 0; i < size; ++i)
		tmp += vin[i];

	vout[0] = tmp;

	#ifdef VERBOSE
	print_vector("Reducing: ", size, vin);
	std::cout << " = " << *vout << std::endl;
	#endif
}

#pragma oss task in(vin[0; size]) out(vout[0; 1])
inline void reduce_sum_task(int *vout, const int *vin, size_t size)
{
	*vout = 0;
	for (size_t i = 0; i < size; ++i)
		*vout += vin[i];

	#ifndef NDEBUG
	print_vector("Reducing: ", size, vin, std::cout);
	std::cout << *vout << std::endl;
	#endif
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


#ifdef NANOS6 // ===============================================================

// Distributed Memory
inline void *_rrd_malloc(size_t size, const char info[] = "")
{
	void *ret = nanos6_dmalloc(size, nanos6_equpart_distribution, 0, NULL);
	dbv2printf("%s %s(%lu) = [%p -> %p]\n", info, __func__, size,  ret, (char*)ret + size);

	assert(size == 0 || ret != NULL);
	return ret;
}

inline void _rrd_free(void *in, size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%p, %lu)\n", info, __PRETTY_FUNCTION__, in, size);
	nanos6_dfree(in, size);
}

inline void *_rrl_malloc(size_t size, const char info[] = "")
{
	void *ret = nanos6_lmalloc(size);
	dbv2printf("%s %s(%lu) = [%p -> %p]\n",
	          info, __func__, size, ret, (char*)ret + size);
	assert(size == 0 || ret != NULL);

	return ret;
}

inline void _rrl_free(void *in, size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%p, %lu)\n", info, __func__, in, size);
	nanos6_lfree(in, size);
}

template<typename T, typename Container>
size_t stl_to_global_task(T *vout, const Container &vin)
{
	const size_t sz = vin.size();

	dbvprintf("STL copy %ld elements -> %p\n", sz, (void *) vout);
	if (sz == 0)
		return 0;

	T *tmp = (T *) rrl_malloc(sz * sizeof(T));

	dbvprintf("STL copy %ld elements\n", sz);
	// Copy from container to local memory, this is initialization any way.
	size_t i = 0;
	for (const T &a : vin)
		tmp[i++] = a;

	ompss_memcpy_task(vout, tmp, sz * sizeof(T));

	#pragma oss taskwait
	rrl_free(tmp, sz * sizeof(T));

	return sz;
}

#else // NANOS6 ================================================================

static inline void *_rrd_malloc(size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%ld) = ", info, __func__, size);
	void *ret = malloc(size);
	dbv2printf("[%p -> %p]\n", ret, (char*)ret + size);

	return ret;
}

static inline void _rrd_free(void *in, size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%p, %ld)\n", info, __func__, in, size);
	free(in);
}

static inline void *_rrl_malloc(size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%ld) = ", info, __func__, size);
	void *ret = malloc(size);
	assert(size == 0 || ret != NULL);
	dbv2printf("Using libc lmalloc [%p -> %p] size %ld\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void _rrl_free(void *in, size_t size, const char info[] = "")
{
	dbv2printf("%s %s(%p, %ld)\n", info, __func__, in, size);
	free(in);
}

template<typename T, typename Container>
size_t stl_to_global_task(T *vout, const Container &vin)
{
	const size_t sz = vin.size();

	// Copy from container to local memory, this is initialization any way.
	size_t i = 0;
	for (const T &a : vin)
		vout[i++] = a;

	dbvprintf("STL copy %ld elements\n", sz);

	return sz;
}

#endif // NANOS6 ===============================================================

#endif // OMPSS_UTILS_HPP
