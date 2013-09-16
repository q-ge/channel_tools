#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sparse.h"
#include "dSFMT-src-2.2.1/dSFMT.h"
#include "channel_algorithms.h"

int
sample(dv_t *cp, dsfmt_t *rng) {
    float x= dsfmt_genrand_close_open(rng);
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

dv_t *
average_cols(csc_mat_t *M) {
    dv_t *v;
    int64_t i;
    float Pc;

    v= dv_new(M->nrow);
    if(!v) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }
    dv_zero(v);

    /* Tally all matrix elements by row. */
    for(i= 0; i < M->nnz; i++)
        v->entries[M->rows[i]]+= M->entries[i];

    /* Scale the rows down, and find cumulative prob. */
    Pc= 0.0;
    for(i= 0; i < v->length; i++) {
        Pc+= v->entries[i] / M->nrow;
        v->entries[i]= Pc;
    }

    /* Rescale to ensure total prob = 1 */
    for(i= 0; i < v->length; i++)
        v->entries[i]/= Pc;

    return v;
}

void
noisy_matrices(dv_t *cp, float epsilon, int n, int cols,
        int samples, dsfmt_t *rng) {
    int i;

    for(i= 0; i < n; i++) {
        csc_mat_t *M= sampled_matrix(cp, cols, samples, rng);
        float I= blahut_arimoto_4(M, epsilon);
        csc_mat_destroy(M);
        printf("%.12e\n", I);
    }
}

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    dv_t *cum_prob;
    dsfmt_t rng;
    int samples, runs;
    float epsilon;

    if(argc < 5) {
        fprintf(stderr,
                "Usage: %s <channel matrix> <samples per column> "
                "<runs> <max error>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    samples= atoi(argv[2]);
    runs= atoi(argv[3]);
    epsilon= atof(argv[4]);

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

    fprintf(stderr, "Averaging matrix columns...");
    fflush(stderr);
    cum_prob= average_cols(M);
    fprintf(stderr, " done.\n");

    dsfmt_init_gen_rand(&rng, time(NULL));

    fprintf(stderr, "Generating noisy matrices...");
    fflush(stderr);
    noisy_matrices(cum_prob, epsilon, runs, M->ncol, samples, &rng);
    fprintf(stderr, " done.\n");

    dv_destroy(cum_prob);
    csc_mat_destroy(M);

    return 0;
}
