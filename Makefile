CFLAGS=-Wall -Werror

include cpu.mk

ifdef AVX
    CFLAGS+= -march=corei7-avx -mavx
else
    ifdef SSE2
        CFLAGS+= -msse2
    endif
endif

LDFLAGS=-lpthread

ifdef DEBUG
    CFLAGS+= -g
else
    CFLAGS+= -g -O3 -DNDEBUG -ffast-math
endif

ifdef CAP_BENCH
    CFLAGS+= -DCAP_BENCH
endif

EXECUTABLES=speed_sparse channel_matrix analyse capacity analyse_mat \
            mult stride extract_plot
ifdef DEBUG
EXECUTABLES+= test_sparse
endif

SPARSE_OBJS= sparse.o

ifdef AVX
SPARSE_OBJS+= sparse_avx.o
EXECUTABLES+= mvec 
endif

default: ${EXECUTABLES}

cpu.mk:
	./cpu_probe.sh

sparse.o: sparse.c sparse.h
testlib.o: testlib.c testlib.h
channel_algorithms.o: channel_algorithms.c channel_algorithms.h

test_sparse: test_sparse.o ${SPARSE_OBJS} testlib.o
test_sparse: LDFLAGS += -lrt

speed_sparse: speed_sparse.o ${SPARSE_OBJS} testlib.o
speed_sparse: LDFLAGS += -lrt

channel_matrix: channel_matrix.o ${SPARSE_OBJS}

analyse: analyse.o ${SPARSE_OBJS}

analyse_mat: analyse_mat.o ${SPARSE_OBJS}

capacity: capacity.o channel_algorithms.o ${SPARSE_OBJS} log.o fastexp.o
capacity: LDFLAGS += -lm -lrt -lpthread

mvec: mvec.o
mvec: LDFLAGS+= -lrt

mult: mult.o ${SPARSE_OBJS}
mult: LDFLAGS+= -lrt

log_bench: log_bench.o log.o
log_bench: LDFLAGS+= -lrt -lm

sparse_linear.o: sparse_linear.c sparse_linear.h

mult_lin: mult_lin.o ${SPARSE_OBJS} sparse_linear.o
mult_lin: LDFLAGS+= -lrt

extract_plot: extract_plot.o ${SPARSE_OBJS}

stride: stride.o ${SPARSE_OBJS}

.PHONY: clean

clean:
	rm -f *.o ${EXECUTABLES} cpu.mk
