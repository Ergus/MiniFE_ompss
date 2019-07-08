#-----------------------------------------------------------------------

MINIFE_TYPES =  \
        -DMINIFE_SCALAR=double   \
        -DMINIFE_LOCAL_ORDINAL=int      \
        -DMINIFE_GLOBAL_ORDINAL=int

MINIFE_MATRIX_TYPE = -DMINIFE_CSR_MATRIX
# MINIFE_MATRIX_TYPE = -DMINIFE_ELL_MATRIX

#-----------------------------------------------------------------------

MINIFEFLAGS := $(MINIFE_TYPES) $(MINIFE_MATRIX_TYPE)

.PHONY: all debug verbose

all: CFLAGS := -O3 -DNDEBUG -std=c99
all: CXXFLAGS := -O3 -DNDEBUG -std=c++11

debug: CFLAGS := -O0 -g -std=c99 -Wall
debug: CXXFLAGS := -O0 -g -std=c++11 -Wall

verbose: CFLAGS := -O0 -g -std=c99 -Wall -DVERBOSE
verbose: CXXFLAGS := -O0 -g -std=c++11 -Wall -DVERBOSE

nanos: CXX := mcxx
nanos: CFLAGS := -O0 -g -std=c99 -Wall --ompss-2 -DNANOS6 -DVERBOSE
nanos: CXXFLAGS := -O0 -g -std=c++11 -Wall --ompss-2 -DNANOS6 -DVERBOSE

include make_targets
