#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sparse.h"
#include "dSFMT-src-2.2.1/dSFMT.h"

#if 0
int
sample(dv_t *cp, dsfmt_t *rng) {
    double x= dsfmt_genrand_close_open(rng);
    int s= 0, e= cp->length;

    while(s+1 < e) {
        int m= (s+e)/2;

        if(x < cp->entries[m]) e= m;
        else                   s= m;
    }

    return s;
}

#define N 16

void
sample_v(dv_t *cp, dsfmt_t *rng, int *r) {
    double x[N];
    int i;

    dsfmt_fill_array_close_open(rng, x, N);

    for(i= 0; i < N; i++) {
        int s= 0, e= cp->length;

        while(s+1 < e) {
            int m= (s+e)/2;

            if(x[i] < cp->entries[m]) e= m;
            else                      s= m;
        }

        r[i]= s;
    }
}

csc_mat_t *
sampled_matrix(dv_t *cp, int cols, int samples, dsfmt_t *rng) {
    bsc_hist_t *H;
    csc_mat_t *M;
    int c;

    H= bsc_hist_new();
    for(c= 0; c < cols; c++) {
        int i;

#if 0
        for(i= 0; i < samples; i++) {
            int r= sample(cp, rng);
            bsc_hist_count(H, c, r, 1);
        }
#else
        for(i= 0; i < samples/N; i++) {
            int r[N];
            int j;

            sample_v(cp, rng, r);

            for(j= 0; j < N; j++)
                bsc_hist_count(H, c, r[j], 1);
        }
#endif
    }
    M= bsc_normalise(H);
    bsc_hist_destroy(H);
    return M;
}
#endif

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    dv_t *cum_prob;
    dsfmt_t rng;
    int samples;
    struct timespec start, end;
    int ncol, firstrow, lastrow, nrow;
    int *counts;
    int count_total;
    int c;

    if(argc < 3) {
        fprintf(stderr,
                "Usage: %s <channel hist> <samples per column>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    samples= atoi(argv[2]);

    if(fscanf(in, "%d %d %d\n", &ncol, &firstrow, &lastrow) != 3) {
        fprintf(stderr, "Malformed histogram\n");
        exit(EXIT_FAILURE);
    }

    if(ncol < 1 || firstrow < 0 || lastrow < 0 || lastrow < firstrow) {
        fprintf(stderr, "Malformed histogram\n");
        exit(EXIT_FAILURE);
    }

    nrow= lastrow - firstrow + 1;

    counts= calloc(nrow, sizeof(int));
    if(!counts) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    printf("Averaging matrix columns\n");

    {
        int clabel, r, count;
        count_total= 0;
        while(fscanf(in, "%d %d %d\n", &clabel, &r, &count) == 3) {
            if(r < firstrow || r > lastrow || c < 0) {
                fprintf(stderr, "Malformed histogram\n");
                exit(EXIT_FAILURE);
            }
            counts[r - firstrow]+= count;
            count_total+= count;
        }
    }

    printf("%ld\n", count_total);

    fclose(in);

#if 0
    sfmt_init_gen_rand(&rng, time(NULL));

    {
        int i;
        double t;

        clock_gettime(CLOCK_MONOTONIC, &start);
        for(i= 0; i < 100; i++) {
            sampled_matrix(cum_prob, M->ncol, samples, &rng);
        }
        clock_gettime(CLOCK_MONOTONIC, &end);
        t= (end.tv_sec + end.tv_nsec*1e-9)
         - (start.tv_sec + start.tv_nsec*1e-9);
        printf("%.3fms\n", (t / 100) * 1000);
    }

    dv_destroy(cum_prob);
    csc_mat_destroy(M);
#endif

    free(counts);

    return 0;
}
