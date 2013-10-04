#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>

#include "fastexp.h"
#include "log.h"
#include "sparse.h"

#undef DEBUG_CAPACITY

#ifdef DEBUG_CAPACITY
#define D(fmt,...) fprintf(stderr,fmt,## __VA_ARGS__)
#else
#define D(fmt,...)
#endif

#ifdef CAP_BENCH
int iterations;
#endif

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
    } 

    if(e_obs) *e_obs= e_best;
    return Il_best;
}
