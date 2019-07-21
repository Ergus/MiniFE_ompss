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

	int global_nrows;
	int *rows;
	int global_nrow_offsets;
	int *row_offsets;
	int global_nrow_coords;
	int *row_coords;

	int global_nnz;
	int *packed_cols;
	double *packed_coefs;

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

		nvectors(0),

		global_nrows(-1),
		rows(nullptr),
		global_nrow_offsets(-1),
		row_offsets(nullptr),
		global_nrow_coords(-1),
		row_coords(nullptr),

		global_nnz(-1),
		packed_cols(nullptr),
		packed_coefs(nullptr)
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

		if (global_nrows > 0)
			rrd_free(rows, global_nrows * sizeof(int));

		if (global_nrow_offsets > 0)
			rrd_free(row_offsets, global_nrow_offsets * sizeof(int));

		if(global_nrow_coords > 0)
			rrd_free(row_coords, global_nrow_coords * sizeof(int));

		if (global_nnz > 0) {
			rrd_free(packed_cols, global_nnz * sizeof(int));
			rrd_free(packed_coefs, global_nnz * sizeof(double));
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

	template <typename T>
	void allocate_rows(T *A_array)
	{
		allocate_lambda<T, int>(global_nrows, rows,
		                A_array,
		                [](const T &A){return A.nrows;},
		                [](T &A, int *val){A.rows = val;});

		allocate_lambda<T, int>(global_nrow_offsets, row_offsets,
		                A_array,
		                [](const T &A){return A.nrows + 1;},
		                [](T &A, int *val){A.row_offsets = val;});

		allocate_lambda<T, int>(global_nrow_coords, row_coords,
		                A_array,
		                [](const T &A){return 3 * A.nrows;},
		                [](T &A, int *val){A.row_coords = val;});
	}

	template <typename T>
	void allocate_packed(T *A_array, size_t numboxes)
	{
		allocate_lambda<T, int>(global_nnz, packed_cols,
		                        A_array,
		                        [](const T &A){return A.nnz;},
		                        [](T &A, int *val){A.packed_cols = val;});

		allocate_lambda<T, double>(global_nnz, packed_coefs,
		                           A_array,
		                           [](const T &A){return A.nnz;},
		                           [](T &A, double *val){A.packed_coefs = val;});

	}

private:
	template <typename T, typename R>
	void allocate_lambda(int &total_size, R *&buffer,
	                     T *A_array,
	                     std::function<int(const T &)> fun1,
	                     std::function<void(T &, R*)> fun2)
	{
		int *indices = (int *) alloca(numboxes * sizeof(int));
		int *offsets = (int *) alloca(numboxes * sizeof(int));

		get_offsets<T, R>(total_size, indices, offsets, numboxes, A_array, fun1);

		buffer = (R *) rrd_malloc(total_size * sizeof(R));

		for (size_t id = 0; id < numboxes; ++id)
			fun2(A_array[id], &buffer[offsets[id]]);
	}
};


#endif
