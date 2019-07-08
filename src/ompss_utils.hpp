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
#endif

// General use macros Here
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
#else
#define dbvprintf(...)
#endif

#ifndef NDEBUG
#define dbarray_to_stream(ARRAY, size, ...) \
	array_to_stream(ARRAY, size, #ARRAY, ##__VA_ARGS__)
#else
#define dbarray_to_stream(...)
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
void print_vector(std::string vname, size_t size,const  T *vect ,std::ostream &stream)
{
	stream << vname << "["  << size << "]={";
	for (size_t i = 0; i < size; ++i) {
		if (i > 0)
			stream << "; ";

		stream << vect[i];
	}
	stream << " }\n";
}

template <typename _Key, typename _Tp>
std::ostream& operator<< (std::ostream& os, const std::pair<_Key, _Tp> &in)
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


#ifdef NANOS6 // ===============================================================

// Nanos6 defined (this can go in a file)
#include "nanos6.h"

// Distributed Memory
static inline void *rrd_malloc(size_t size)
{
	dbvprintf("Using nanos6_dmalloc ");
	void *ret = nanos6_dmalloc(size, nanos6_equpart_distribution, 0, NULL);
	assert(size == 0 || ret != NULL);

	dbvprintf("[%p -> %p] size %ld\n", ret, (char*)ret + size, size);
	return ret;
}

static inline void rrd_free(void *in, size_t size)
{
	dbvprintf("Using nanos6_dfree(%p)\n", in);
	if (in)
		nanos6_dfree(in, size);
}

// Local Memory
static inline void *rrl_malloc(size_t size)
{
	void *ret = nanos6_lmalloc(size);
	assert(size == 0 || ret != NULL);
	dbvprintf("Using nanos6_lmalloc [%p -> %p] size %ld\n",
	          ret, (char*)ret + size, size);

	return ret;
}

static inline void *rrl_calloc(size_t nmemb, size_t size)
{
	const int bytes = nmemb * size;
	void *ret = rrl_malloc(bytes);
	assert(size == 0 || ret != NULL);
	dbvprintf("Using nanos6_lmalloc (calloc) [%p -> %p] size %ld\n",
		 ret, (char*)ret + size, size);
	memset(ret, 0, bytes);

	return ret;
}


static inline void rrl_free(void *in, size_t size)
{
	dbvprintf("Using nanos6_lfree(%p)\n", in);
	if (in)
		nanos6_lfree(in, size);
}

inline void ompss_memcpy_task(void *pout, const void *pin, size_t size)
{
	char *tin = (char *) pin;
	char *tout = (char *) pout;

	dbvprintf("Copy %ld bytes from %p -> %p\n", size, pin, pout);
	#pragma oss task in(tin[0; size]) out(tout[0; size])
	{
		memcpy(tout, tin, size);
	}
}

template <typename T>
void reduce_sum_task(T *vout, const T *vin, size_t size)
{
	#pragma oss task in(vin[0; size]) out(*vout)
	{
		*vout = 0;
		for (size_t i = 0; i < size; ++i)
			*vout += vin[i];
	}
}

template<typename T, typename Container>
size_t stl_to_global_task(T *vout, const Container &vin)
{
	const size_t sz = vin.size();

	dbvprintf("STL copy %ld elements -> %p\n", sz, (void *) vout);
	if (sz == 0)
		return 0;

	T *tmp = (T *) rrl_malloc(sz * sizeof(T));

	dbvprintf("STL allocated %ld elements -> %p\n");
	// Copy from container to local memory, this is initialization any way.
	size_t i = 0;
	for (const T &a : vin)
		tmp[i++] = a;

	ompss_memcpy_task(vout, tmp, sz * sizeof(T));

	#pragma oss taskwait
	rrl_free(tmp, sz * sizeof(T));

	return sz;
}

#define get_node_id() nanos6_get_cluster_node_id()
#define get_nodes_nr() nanos6_get_cluster_nodes()

#else // NANOS6 ================================================================

// use libc functions

static inline void *rrd_malloc(size_t size)
{
	void *ret = malloc(size);
	dbvprintf("Using libc dmalloc [%p -> %p] size %ld\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void rrd_free(void *in, size_t size  __attribute__((unused)))
{
	dbvprintf("Using libc dfree(%p)\n", in);
	free(in);
}

static inline void *rrl_malloc(size_t size)
{
	void *ret = malloc(size);
	assert(size == 0 || ret != NULL);
	dbvprintf("Using libc lmalloc [%p -> %p] size %ld\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void *rrl_calloc(size_t nmemb, size_t size)
{
	void *ret = calloc(nmemb, size);
	assert(size == 0 || ret != NULL);
	dbvprintf("Using libc lcalloc [%p -> %p] size %ld\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void rrl_free(void *in, size_t size  __attribute__((unused)))
{
	dbvprintf("Using libc lfree(%p)\n", in);
	free(in);
}

inline void ompss_memcpy_task(void *pout, const void *pin, size_t size)
{
	memcpy(pout, pin, size);
}

template <typename T>
void reduce_sum_task(T *out, const T *in, size_t size)
{
	*out = 0;
	for (size_t i = 0; i < size; ++i)
		*out += in[i];
}

template<typename T, typename Container>
size_t stl_to_global_task(T *vout, const Container &vin)
{
	const size_t sz = vin.size();

	// Copy from container to local memory, this is initialization any way.
	size_t i = 0;
	for (const T &a : vin)
		vout[i++] = a;

	dbvprintf("Copy %ld bytes -> %p\n", sz, (void *) vout);

	return sz;
}

#define get_node_id() 0
#define get_nodes_nr() 1

#endif // NANOS6 ===============================================================

#endif // OMPSS_UTILS_HPP
