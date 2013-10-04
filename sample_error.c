#include <assert.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "log.h"
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
sampled_matrix(dv_t *cp, int rows, int samples, dsfmt_t *rng) {
    bsc_hist_t *H;
    csc_mat_t *M;
    int r;

    H= bsc_hist_new();
    for(r= 0; r < rows; r++) {
        int i;

        for(i= 0; i < samples; i++) {
            int c= sample(cp, rng);
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
    int c;
    float Ps= 0.0;

    v= dv_new(M->ncol);
    if(!v) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }
    dv_zero(v);

    /* Sum along every column. */
    for(c= 0; c < M->ncol; c++) {
        int64_t i;

        for(i= M->ci[c]; i < M->ci[c+1]; i++) {
            if(M->entries[i] == 0) continue;
            v->entries[c]+= M->entries[i];
        }
    }

    /* Rescale to ensure total prob = 1 */
    for(c= 0; c < v->length; c++) Ps+= v->entries[c];
    for(c= 0; c < v->length; c++) v->entries[c]/= Ps;

    return v;
}

dv_t *
accumulate_prob(dv_t *p) {
    dv_t *cp;
    float Pc= 0.0;
    int i;

    cp= dv_new(p->length);
    if(!cp) {
        perror("dv_new");
        exit(EXIT_FAILURE);
    }

    for(i= 0; i < p->length; i++) {
        Pc+= p->entries[i];
        cp->entries[i]= Pc;
    }

    return cp;
}

pthread_mutex_t output_lock= PTHREAD_MUTEX_INITIALIZER;

void
noisy_matrices(dv_t *cp, float epsilon, int n, int rows,
        int samples, dsfmt_t *rng) {
    int i;

    for(i= 0; i < n; i++) {
        csc_mat_t *M= sampled_matrix(cp, rows, samples, rng);
        float I, e;
        csc_check(M, 1);
        I= blahut_arimoto(M, epsilon, &e);
        csc_mat_destroy(M);
        pthread_mutex_lock(&output_lock);
        printf("%.12e %.12e\n", I, e);
        pthread_mutex_unlock(&output_lock);
    }
}

struct job {
    dv_t *cp;
    float epsilon;
    int n;
    int rows;
    int samples;
    dsfmt_t rng;
    pthread_t thread;
};

static void *
worker(void *arg) {
    struct job *j= (struct job *)arg;
    noisy_matrices(j->cp, j->epsilon, j->n, j->rows, j->samples, &j->rng);
    return (void *)0;
}

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    dv_t *avg_prob, *cum_prob;
    int samples, runs, runs_per_job;
    float epsilon;
    int nthreads;
    struct job *jobs;
    int seed;
    int quiet= 0;
    int i;

    srandom(time(NULL));

    if(argc < 5) {
        fprintf(stderr,
                "Usage: %s <channel matrix> <samples per column> "
                "<runs> <max error> [-q]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    samples= atoi(argv[2]);
    runs= atoi(argv[3]);
    epsilon= atof(argv[4]);

    if(argc > 5 && !strcmp(argv[5], "-q"))
        quiet= 1;

    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);

    if(!csc_check(M, 1)) abort();
    csc_prune_cols(M);
    if(!quiet) csc_stats(M);

    if(STRIDE_OF(M) > 1) {
        fprintf(stderr, "Can't handle a strided matrix\n");
        exit(EXIT_FAILURE);
    }

    if(!quiet) fprintf(stderr, "Averaging matrix columns...");
    if(!quiet) fflush(stderr);
    avg_prob= average_cols(M);
    cum_prob= accumulate_prob(avg_prob);
    if(!quiet) fprintf(stderr, " done.\n");

    nthreads= sysconf(_SC_NPROCESSORS_ONLN);
    if(nthreads < 1) nthreads= 1;

    write_log_table();

    jobs= malloc(nthreads * sizeof(struct job));
    if(!jobs) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    if(!quiet) fprintf(stderr, "Generating noisy matrices...");
    if(!quiet) fflush(stderr);
    runs_per_job= runs / nthreads;
    for(i= 0; i < nthreads; i++) {
        jobs[i].cp= cum_prob;
        jobs[i].epsilon= epsilon;
        jobs[i].n= runs_per_job;
        jobs[i].rows= M->nrow;
        jobs[i].samples= samples;
        seed= random();
        dsfmt_init_gen_rand(&jobs[i].rng, seed);
        runs-= runs_per_job;
    }
    jobs[nthreads-1].n+= runs;
    for(i= 0; i < nthreads; i++) {
        if(pthread_create(&jobs[i].thread, NULL, worker, &jobs[i])) {
            perror("pthread_create");
            exit(EXIT_FAILURE);
        }
    }
    for(i= 0; i < nthreads; i++) {
        if(pthread_join(jobs[i].thread, NULL)) {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }
    if(!quiet) fprintf(stderr, " done.\n");

    dv_destroy(cum_prob);
    csc_mat_destroy(M);

    return 0;
}
