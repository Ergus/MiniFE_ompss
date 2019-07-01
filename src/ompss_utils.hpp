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

#include <cstdio>
#include <cstdlib>
#include <cassert>

// General use macros Here
#define REGISTER_ELAPSED_TIME(time_inc, time_total)			\
	{								\
		time_inc = mytimer() - time_inc;			\
		time_total += time_inc;					\
	}

// Debug conditional macros here.
#ifdef NDEBUG
#define dbprintf(...)
#else
#define dbprintf(...) fprintf(stderr, __VA_ARGS__)
#endif

template <typename T>
void reduce_sum(T &out, const T *in, size_t size)
{
	out = 0;
	for (size_t i = 0; i < size; ++i)
		out += in[i];
}

#ifdef NANOS6 // ===============================================================

// Nanos6 defined (this can go in a file)
#include "nanos6.h"

// Distributed Memory
static inline void *rrd_malloc(size_t size)
{
	void *ret = nanos6_dmalloc(size, nanos6_equpart_distribution, 1, NULL);
	assert(ret != NULL);
	dbprintf("Using nanos6_dmalloc [%p -> %p] size %d\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void rrd_free(void *in, size_t size)
{
	dbprintf("Using nanos6_dfree(%p)\n", in);
	nanos6_dfree(in, size);
}

// Local Memory
static inline void *rrl_malloc(size_t size)
{
	void *ret = nanos6_lmalloc(size);
	assert(ret != NULL);
	dbprintf("Using nanos6_lmalloc [%p -> %p] size %d\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void *rrl_calloc(size_t nmemb, size_t size)
{
	const int bytes = nmemb * size;

	void *ret = rrl_malloc(bytes);
	memset(ret, 0, bytes);

	return ret;
}


static inline void rrl_free(void *in, size_t size)
{
	dbprintf("Using nanos6_lfree(%p)\n", in);
	nanos6_lfree(in, size);
}

template<typename T>
void copy_local_to_global_task(T *vout, const T *vin, size_t size)
{
	const size_t nbytes = size * sizeof(T);

	// TODO: Task here
	// in vin[0; size]
	// out vout[0; size]
	{
		memcpy(vout, vin, nbytes);
	}

}

#define get_node_id() nanos6_get_cluster_node_id()
#define get_nodes_nr() nanos6_get_cluster_nodes()

#else // NANOS6 ================================================================

// use libc functions

static inline void *rrd_malloc(size_t size)
{
	void *ret = malloc(size);
	assert(ret != NULL);
	dbprintf("Using libc dmalloc [%p -> %p] size %d\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void rrd_free(void *in, size_t size  __attribute__((unused)))
{
	dbprintf("Using libc dfree(%p)\n", in);
	free(in);
}

static inline void *rrl_malloc(size_t size)
{
	void *ret = malloc(size);
	assert(ret != NULL);
	dbprintf("Using libc lmalloc [%p -> %p] size %d\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void *rrl_calloc(size_t nmemb, size_t size)
{
	void *ret = calloc(nmemb, size);
	assert(ret != NULL);
	dbprintf("Using libc lcalloc [%p -> %p] size %d\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void rrl_free(void *in, size_t size  __attribute__((unused)))
{
	dbprintf("Using libc lfree(%p)\n", in);
	free(in);
}

template<typename T>
void copy_local_to_global_task(T *vout, const T *vin, size_t size)
{
	memcpy(vout, vin, size * sizeof(T));
}


#define get_node_id() 0
#define get_nodes_nr() 1

#endif // NANOS6 ===============================================================

#endif // OMPSS_UTILS_HPP
