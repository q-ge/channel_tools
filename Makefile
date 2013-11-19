CFLAGS=-Wall# -Werror
CFLAGS+=-DDSFMT_MEXP=521

CFLAGS+=${EXTRA_CFLAGS}

dSFMT_SRC=dSFMT-src-2.2.1

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

DEBUG_EXECUTABLES= test_sparse test_hist

EXECUTABLES=speed_sparse channel_matrix analyse capacity analyse_mat \
            mult stride extract_plot sample_error channel_hist \
            summarise filter_samples drop_samples row_average \
            smooth
ifdef DEBUG
EXECUTABLES+= $(DEBUG_EXECUTABLES)
endif

SPARSE_OBJS= sparse.o

ifdef AVX
SPARSE_OBJS+= sparse_avx.o
EXECUTABLES+= mvec 
endif

TESTS=test/matrix1

.PHONY: default test clean

default: ${EXECUTABLES}

cpu.mk:
	./cpu_probe.sh

sparse.o: sparse.c sparse.h
testlib.o: testlib.c testlib.h
channel_algorithms.o: channel_algorithms.c channel_algorithms.h

test_sparse: test_sparse.o ${SPARSE_OBJS} testlib.o
test_sparse: LDFLAGS += -lrt
test_sparse.o: CFLAGS= -Wall -g -DDEBUG

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

sample_error: LDFLAGS+= -lm -lrt -lpthread
sample_error: sample_error.o ${SPARSE_OBJS} channel_algorithms.o \
              log.o fastexp.o ${dSFMT_SRC}/dSFMT.o

channel_hist: channel_hist.o ${SPARSE_OBJS}

test_hist: test_hist.o ${SPARSE_OBJS}
test_hist.o: CFLAGS= -Wall -g -DDEBUG

row_average: row_average.o ${SPARSE_OBJS}

### Tests

TEST_MATRICES=matrix1

TEST_TARGETS= \
    $(patsubst %,test/%.plot,${TEST_MATRICES}) \
    $(patsubst %,test/%.capacity,${TEST_MATRICES}) \
    $(patsubst %,test/%.sim,${TEST_MATRICES})

%.cm: %.samples.xz channel_matrix
	xzcat $< | ./channel_matrix $@

%.plot: %.cm extract_plot
	./extract_plot $< 1024 1024 > $@

%.capacity: %.cm capacity
	./capacity $< 1e-3 > $@

%.sim: %.cm sample_error
	./sample_error $< 100 10 1e-3 1 0 1 > $@

test: ${TEST_TARGETS} test_sparse
	./test_sparse

###

clean:
	rm -f *.o ${dSFMT_SRC}/*.o ${EXECUTABLES} ${DEBUG_EXECUTABLES} \
              ${TESTS} cpu.mk ${TEST_TARGETS}
