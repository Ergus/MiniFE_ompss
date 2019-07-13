#ifndef SINGLETON_HPP
#define SINGLETON_HPP

#include "ompss_utils.hpp"

#define MAXVECTORS 10

class singleton {
public:
	size_t numboxes;

	// Info recv
	int global_nrecv_neighbors;
	int global_nexternals;

	int *nrecv_neighbors;   // neighbors to receive from / process
	int *recv_neighbors;    // list of neighbors to receive from
	double **recv_ptr;      // List of remote pointers to copy from
	int *recv_length;       // List of lengths to copy / neighbors

	int *nexternals;        // number of external elements/ process
	int *external_index;

	// Info sent
	int global_nsend_neighbors; // total neighbors processes.
	int global_nelements_to_send; // total elements

	int *nsend_neighbors;       // total neighbors processes / process
	int *send_neighbors;        // full list of neighbors to send
	int *send_length;           // full list of lengths / neighbors / process

  	int *nelements_to_send;     // total elements to send / process
	int *elements_to_send;      // full list of local indices
	double *send_buffer;        // send buffer

	size_t nvectors;
	int vector_global_size[MAXVECTORS]; // TODO: this is totally arbitrary
	double *vector_coefs[MAXVECTORS];

	singleton(size_t _numboxes) :
		numboxes(_numboxes),
		global_nrecv_neighbors(-1),
		global_nexternals(-1),

		nrecv_neighbors((int *) rrl_malloc(_numboxes * sizeof(int))),
		recv_neighbors(nullptr),
		recv_ptr(nullptr),
		recv_length(nullptr),

		nexternals((int *) rrl_malloc(_numboxes * sizeof(int))),
		external_index(nullptr),

		global_nsend_neighbors(-1),
		global_nelements_to_send(-1),

		nsend_neighbors((int *) rrl_malloc(_numboxes * sizeof(int))),
		send_neighbors(nullptr),
		send_length(nullptr),

		nelements_to_send(((int *) rrl_malloc(_numboxes * sizeof(int)))),
		elements_to_send(nullptr),
		send_buffer(nullptr),

		nvectors(0)
	{
		for (size_t i = 0 ; i < MAXVECTORS; ++i) {
			vector_global_size[i] = -1;
			vector_coefs[i] = nullptr;
		}
	}

	~singleton()
	{
		// recv
		rrl_free(nrecv_neighbors, numboxes * sizeof(int));
		rrl_free(nexternals, numboxes * sizeof(int));

		if (global_nrecv_neighbors > 0) {
			rrd_free(recv_neighbors, global_nrecv_neighbors * sizeof(int));
			rrd_free(recv_ptr, global_nrecv_neighbors * sizeof(double *));
			rrd_free(recv_length, global_nrecv_neighbors * sizeof(int));
		}

		if (global_nexternals)
			rrd_free(external_index, global_nexternals * sizeof(int));

		// send
		rrl_free(nsend_neighbors, numboxes * sizeof(int));
		rrl_free(nelements_to_send, numboxes * sizeof(int));

		if (global_nsend_neighbors) {
			rrd_free(send_neighbors, global_nsend_neighbors * sizeof(int));
			rrd_free(send_length, global_nsend_neighbors * sizeof(int));
		}

		if (global_nelements_to_send > 0) {
			rrd_free(elements_to_send, global_nelements_to_send * sizeof(int));
			rrd_free(send_buffer, global_nelements_to_send * sizeof(double));
		}

		for (size_t i = 0; i < nvectors; ++i) {
			if (vector_coefs[i] != nullptr) {
				assert(vector_global_size[i] >= 0);
				rrd_free(vector_coefs[i], vector_global_size[nvectors] * sizeof(double));
			}
		}

	}


	bool allocate_recv(int _global_nrecv_neighbors, int _global_nexternals)
	{
		assert (!(recv_neighbors || recv_ptr || recv_length || external_index));

		global_nrecv_neighbors = _global_nrecv_neighbors;
		global_nexternals = _global_nexternals;

		recv_neighbors = (int *) rrd_malloc(global_nrecv_neighbors * sizeof(int));
		recv_length = (int *) rrd_malloc(global_nrecv_neighbors * sizeof(int));
		recv_ptr = (double **) rrd_malloc(global_nrecv_neighbors * sizeof(double *));

		external_index = (int *) rrd_malloc(global_nexternals * sizeof(int));

		return (recv_neighbors && recv_ptr && recv_length && external_index);
	}


	bool allocate_send(int _global_nsend_neighbors, int _global_nelements_to_send)
	{
		assert(!(send_neighbors || send_length || elements_to_send || send_buffer));

		global_nsend_neighbors = _global_nsend_neighbors;
		global_nelements_to_send = _global_nelements_to_send;

		send_neighbors = (int *) rrd_malloc(global_nsend_neighbors * sizeof(int));
		send_length = (int *) rrd_malloc(global_nsend_neighbors * sizeof(int));

		elements_to_send = (int *) rrd_malloc(global_nelements_to_send * sizeof(int));
		send_buffer = (double *) rrd_malloc(global_nelements_to_send * sizeof(double));

		return (send_neighbors && send_length && elements_to_send && send_buffer);
	}

	double *allocate_vectors(size_t total_size)
	{
		if (nvectors == MAXVECTORS)
			return nullptr;

		vector_global_size[nvectors] = total_size;
		vector_coefs[nvectors] = (double *) rrd_malloc(total_size * sizeof(double));

		return vector_coefs[nvectors++];
	}

	void* operator new(size_t sz)
	{
		void * const tmp = rrl_malloc(sz);
		dbvprintf("Calling: %s, size: %lu\n", __PRETTY_FUNCTION__, sz);
		return tmp;
	}

	static void operator delete(void* ptr, std::size_t sz)
	{
		dbvprintf("Calling: %s, address %p size: %lu\n", __PRETTY_FUNCTION__, ptr, sz);
		return rrl_free(ptr, sz);
	}

};


#endif
