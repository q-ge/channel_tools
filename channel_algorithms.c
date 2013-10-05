#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "fastexp.h"
#include "log.h"
#include "sparse.h"

#define DEBUG_CAPACITY

#ifdef DEBUG_CAPACITY
#define D(fmt,...) fprintf(stderr,fmt,## __VA_ARGS__)
#else
#define D(fmt,...)
#endif

#ifdef CAP_BENCH
int iterations;
#endif

void
dv_normalise(dv_t *v) {
    float Ps;
    int i;

    Ps= 0;
    for(i= 0; i < v->length; i++) Ps+= v->entries[i];
    for(i= 0; i < v->length; i++) v->entries[i]/= Ps;
}

float
blahut_arimoto_precise(csc_mat_t *Q, float epsilon, float *e_obs) {
    double *p, *q, *c;
    double Il, Iu, e;
    double Il_best, e_best= INFINITY;
    int i, col, row;
    double *logQ, *logq;

    D("\n");

    p= malloc(Q->nrow * sizeof(double));
    q= malloc(Q->ncol * sizeof(double));
    c= malloc(Q->nrow * sizeof(double));
    if(!p || !q || !c) { perror("malloc"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(double));
    if(!logQ) { perror("malloc"); abort(); }
    for(i= 0; i < Q->nnz; i++) logQ[i]= log2((double)Q->entries[i]);
    D(" done.\n");

    logq= malloc(Q->ncol * sizeof(double));
    if(!logq) { perror("malloc"); abort(); }

    /* Start with a uniform input distribution. */
    for(i= 0; i < Q->nrow; i++) p[i]= 1.0 / Q->nrow;

#ifdef CAP_BENCH
    iterations= 0;
#endif

    while(1) {
#ifdef CAP_BENCH
        iterations++;
#endif

        /* Calculate q - output distribution. */
        {
            int c;

            for(c= 0; c < Q->ncol; c++) {
                int i;

                q[c]= 0.0;
                for(i= Q->ci[c]; i < Q->ci[c+1]; i++)
                    q[c]+= (double)Q->entries[i] * p[Q->rows[i]];
            }
        }

        /* Precalculate log(q). */
        for(i= 0; i < Q->ncol; i++) logq[i]= log2(q[i]);

        /* Calculate log(c). */
        bzero(c, Q->nrow * sizeof(double));
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                c[row]+= (double)Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= -INFINITY;
        for(i= 0; i < Q->nrow; i++) {
            if(c[i] > Iu) Iu= c[i];
        }

        {
            double tmp= 0.0;

            /* Update Il.  By taking out the largest element of c, the
             * rest are all less than 1. */
            for(row= 0; row < Q->nrow; row++) {
                tmp+= p[row] * exp2(c[row] - Iu);
            }
            Il= Iu + log2(tmp);
        }

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6e Iu= %.6e e= %.6e\n", Il, Iu, e);
        if(e < epsilon) break;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < Q->nrow; i++) {
            p[i]*= exp2(c[i] - Il);
        }
    } 

    if(e_obs) *e_obs= (float)e_best;
    return (float)Il_best;
}

/* Fourth attempt. Use log tables. */
float
blahut_arimoto(csc_mat_t *Q, float epsilon, float *e_obs) {
    dv_t *p, *q, *c;
    float Il, Iu, e;
    float Il_last= 0, Iu_last= INFINITY;
    float Il_best= 0, e_best= INFINITY;
    int i, col, row;
    float *logQ, *logq;

    D("\n");

    p= dv_new(Q->nrow);
    q= dv_new(Q->ncol);
    c= dv_new(Q->nrow);
    if(!p || !q || !c) { perror("dv_new"); abort(); }

    D("Precalculating logs...");
    logQ= malloc(Q->nnz * sizeof(float));
    if(!logQ) { perror("malloc"); abort(); }
    for(i= 0; i < Q->nnz; i++) logQ[i]= log2f_table(Q->entries[i]);
    D(" done.\n");

    logq= malloc(q->length * sizeof(float));
    if(!logq) { perror("malloc"); abort(); }

    /* Start with a uniform input distribution. */
    dv_uniform(p, 1.0 / Q->nrow);

#ifdef CAP_BENCH
    iterations= 0;
#endif

    while(1) {
#ifdef CAP_BENCH
        iterations++;
#endif

        /* Calculate q - output distribution. */
        mult_csc_dv(q, p, Q);

        dv_normalise(q);

        /* Precalculate log(q). */
        for(i= 0; i < q->length; i++) logq[i]= log2f_table(q->entries[i]);

        /* Calculate log(c). */
        dv_zero(c);
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);

                /* Use the log tables. */
                c->entries[row]+= Q->entries[i] * (logQ[i] - logq[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= dv_max(c);

        {
            float tmp= 0.0;

            /* Update Il.  By taking out the largest element of c, the
             * rest are all less than 1. */
            for(row= 0; row < p->length; row++) {
                tmp+= p->entries[row] * exp2f_fast(c->entries[row] - Iu);
            }
            Il= Iu + log2f_table(tmp);
        }

        e= Iu - Il;

        if(e < e_best) {
            Il_best= Il; e_best= e;
        }

        D("Il= %.6e Iu= %.6e e= %.6e\n", Il, Iu, e);
        if(e < epsilon || (Il==Il_last && Iu==Iu_last)) break;

        Il_last= Il; Iu_last= Iu;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < p->length; i++) {
            p->entries[i]*= exp2f_fast(c->entries[i] - Il);
        }

        dv_normalise(p);
    } 

    if(e_obs) *e_obs= e_best;
    return Il_best;
}
