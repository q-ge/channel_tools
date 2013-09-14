#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sparse.h"
#include "dSFMT-src-2.2.1/dSFMT.h"

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

csc_mat_t *
sampled_matrix(dv_t *cp, int cols, int samples, dsfmt_t *rng) {
    bsc_hist_t *H;
    csc_mat_t *M;
    int c;

    H= bsc_hist_new();
    for(c= 0; c < cols; c++) {
        int i;

        for(i= 0; i < samples; i++) {
            int r= sample(cp, rng);
            bsc_hist_count(H, c, r, 1);
        }
    }
    M= bsc_normalise(H);
    bsc_hist_destroy(H);
    return M;
}

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    dv_t *cum_prob;
    dsfmt_t rng;
    int samples;

    if(argc < 3) {
        fprintf(stderr,
                "Usage: %s <channel matrix> <samples per column>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    samples= atoi(argv[2]);

    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);

    if(!csc_check(M, 1)) abort();
    csc_prune_cols(M);
    csc_stats(M);

    if(STRIDE_OF(M) > 1) {
        fprintf(stderr, "Can't handle a strided matrix\n");
        exit(EXIT_FAILURE);
    }

    printf("Averaging matrix columns\n");

    cum_prob= dv_new(M->nrow);
    if(!cum_prob) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }
    dv_zero(cum_prob);

    {
        int64_t i;
        float Pc, Ps;

        /* Tally all matrix elements by row. */
        for(i= 0; i < M->nnz; i++)
            cum_prob->entries[M->rows[i]]+= M->entries[i];

        /* Scale the rows down, and find cumulative prob. */
        Pc= 0.0;
        for(i= 0; i < cum_prob->length; i++) {
            Pc+= cum_prob->entries[i] / M->nrow;
            cum_prob->entries[i]= Pc;
        }

        for(i= 0; i < cum_prob->length; i++)
            cum_prob->entries[i]/= Pc;

    }

    dsfmt_init_gen_rand(&rng, time(NULL));

    sampled_matrix(cum_prob, M->ncol, samples, &rng);

    dv_destroy(cum_prob);
    csc_mat_destroy(M);

    return 0;
}
