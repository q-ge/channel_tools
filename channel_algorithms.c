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

/* A Naïve implementation of the original algorithm.  Blows up on large
   channel matrices, as the exponentials get enormous. */
float
blahut_arimoto_1(csc_mat_t *Q, float epsilon) {
    dv_t *p, *q, *c;
    float Il, Iu;
    int i, col, row;

    p= dv_new(Q->nrow);
    q= dv_new(Q->ncol);
    c= dv_new(Q->nrow);
    if(!p || !q || !c) { perror("dv_new"); abort(); }

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

        /* Calculate c */
        dv_uniform(c, 1.0);
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                assert(Q->entries[i] > 0);
                assert(q->entries[col] > 0);

                c->entries[row]*=
                    exp2f(Q->entries[i] *
                          log2f(Q->entries[i] /
                                q->entries[col]));
            }
        }

        /* Update Il and Iu. */
        Il= log2f(dv_dot(p, c));
        Iu= log2f(dv_max(c));

        D("Il= %.3f Iu= %.3f\n", Il, Iu);
        if(Iu - Il < epsilon) break;

        /* Update input distribution. */
        for(i= 0; i < p->length; i++) p->entries[i]*= c->entries[i] /
                                                      exp2f(Il);
    } 

    return Il;
}

/* Second attempt. Keep the terms in log form as long as possible, then
 * prescale them. */
float
blahut_arimoto_2(csc_mat_t *Q, float epsilon) {
    dv_t *p, *q, *c;
    float Il, Iu;
    int i, col, row;

    p= dv_new(Q->nrow);
    q= dv_new(Q->ncol);
    c= dv_new(Q->nrow);
    if(!p || !q || !c) { perror("dv_new"); abort(); }

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

        /* Calculate log(c). */
        dv_zero(c);
        for(col= 0; col < Q->ncol; col++) {
            for(i= Q->ci[col]; i < Q->ci[col+1]; i++) {
                row= Q->rows[i];

                /* Excessively small values here will blow up
                   under division, but have very little effect on
                   the outcome. */
                if(q->entries[col] < 1e-40) continue;

                assert(Q->entries[i] > 0);

                c->entries[row]+= Q->entries[i] *
                                  log2f(Q->entries[i] /
                                        q->entries[col]);

            }
        }

        /* Find the largest element of c.  This is Iu, as log is monotonic. */
        Iu= dv_max(c);

        {
            float tmp= 0.0;

            /* Update Il.  By taking out the largest element of c, the
             * rest are all less than 1. */
            for(row= 0; row < p->length; row++) {
                tmp+= p->entries[row] * exp2f(c->entries[row] - Iu);
            }
            Il= Iu + log2f(tmp);
        }

        D("Il= %.3f Iu= %.3f\n", Il, Iu);
        if(Iu - Il < epsilon) break;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < p->length; i++) {
            p->entries[i]*= exp2f(c->entries[i] - Il);
        }
    } 

    return Il;
}

/* Third attempt. Precalculate logs. */
float
blahut_arimoto_3(csc_mat_t *Q, float epsilon) {
    dv_t *p, *q, *c;
    float Il, Iu;
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
    for(i= 0; i < Q->nnz; i++) logQ[i]= log2f(Q->entries[i]);
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
        for(i= 0; i < q->length; i++) logq[i]= log2f(q->entries[i]);

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
                tmp+= p->entries[row] * exp2f(c->entries[row] - Iu);
            }
            Il= Iu + log2f(tmp);
        }

        D("Il= %.3f Iu= %.3f\n", Il, Iu);
        if(Iu - Il < epsilon) break;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < p->length; i++) {
            p->entries[i]*= exp2f(c->entries[i] - Il);
        }
    } 

    return Il;
}

/* Fourth attempt. Use log tables. */
float
blahut_arimoto_4(csc_mat_t *Q, float epsilon) {
    dv_t *p, *q, *c;
    float Il, Iu;
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

        D("Il= %.3f Iu= %.3f\n", Il, Iu);
        if(Iu - Il < epsilon) break;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < p->length; i++) {
            p->entries[i]*= exp2f_fast(c->entries[i] - Il);
        }
    } 

    return Il;
}

#define NTHREADS 2

struct thread_data {
    pthread_t id;

    pthread_cond_t start;
    pthread_mutex_t m_start;

    pthread_cond_t finished;
    pthread_mutex_t m_finished;

    int scol, ecol;

    /* Read-only, shared. */
    csc_mat_t *Q;
    dv_t *p, *q;
    float *logQ, *logq;

    /* Read-write, private. */
    dv_t *c;
};

void
start(struct thread_data *td, int i) {
    if(pthread_mutex_lock(&td[i].m_start))
        { perror("pthread_mutex_lock"); abort(); }

    if(pthread_cond_broadcast(&td[i].start))
        { perror("pthread_cond_broadcast"); abort(); }

    if(pthread_mutex_unlock(&td[i].m_start))
        { perror("pthread_mutex_unlock"); abort(); }
}

void
finished(struct thread_data *td) {
    if(pthread_mutex_lock(&td->m_finished))
        { perror("pthread_mutex_lock"); abort(); }

    if(pthread_cond_broadcast(&td->finished))
        { perror("pthread_cond_broadcast"); abort(); }

    if(pthread_mutex_unlock(&td->m_finished))
        { perror("pthread_mutex_unlock"); abort(); }
}

void wait_for(struct thread_data *td, int i) {
    if(pthread_cond_wait(&td[i].finished, &td[i].m_finished))
        { perror("pthread_cond_wait"); abort(); }
}

static void *
worker(void *_arg) {
    struct thread_data *td= (struct thread_data *)_arg;

    /* First, take the listener mutex for 'start'. */
    if(pthread_mutex_lock(&td->m_start))
        { perror("pthread_mutex_lock"); abort(); }

    /* Then, let the supervisor know that it's safe to signal us. */
    finished(td);

    while(1) {
        int col;

        /* Wait to be woken up. */
        if(pthread_cond_wait(&td->start, &td->m_start))
            { perror("pthread_cond_wait"); abort(); }

        dv_zero(td->c);
        for(col= td->scol; col < td->ecol; col++) {
            int i;
            for(i= td->Q->ci[col]; i < td->Q->ci[col+1]; i++) {
                int row= td->Q->rows[i];

                assert(td->Q->entries[i] > 0);

                /* Use the log tables. */
                td->c->entries[row]+=
                    td->Q->entries[i] * (td->logQ[i] - td->logq[col]);

            }
        }

        /* Signal completion. */
        finished(td);
    }
}

/* Fifth attempt. Threading. */
float
blahut_arimoto_5(csc_mat_t *Q, float epsilon) {
    dv_t *p, *q, *c;
    float Il, Iu;
    int i, row;
    float *logQ, *logq;
    struct thread_data *td;

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

    D("Creating worker threads...");
    td= malloc(NTHREADS * sizeof(struct thread_data));
    if(!td) { perror("malloc"); abort(); }
    for(i= 0; i < NTHREADS; i++) {
        td[i].c= dv_new(Q->nrow);
        if(!td[i].c) { perror("dv_new"); abort(); }

        td[i].scol= (Q->ncol / NTHREADS) * i;
        td[i].ecol= (Q->ncol / NTHREADS) * (i + 1);
        if(td[i].ecol > Q->ncol) td[i].ecol= Q->ncol;

        td[i].Q= Q;
        td[i].p= p;
        td[i].q= q;
        td[i].logQ= logQ;
        td[i].logq= logq;

        if(pthread_cond_init(&td[i].start, NULL))
            { perror("pthread_cond_init"); abort(); }
        pthread_mutex_init(&td[i].m_start, NULL);

        if(pthread_cond_init(&td[i].finished, NULL))
            { perror("pthread_cond_init"); abort(); }
        pthread_mutex_init(&td[i].m_finished, NULL);

        /* Make sure the thread doesn't signal us until
           we're ready. */
        if(pthread_mutex_lock(&td[i].m_finished))
            { perror("pthread_mutex_lock"); abort(); }

        /* Start it. */
        if(pthread_create(&td[i].id, NULL, worker, &td[i]))
            { perror("pthread_create"); abort(); }

        /* Wait for it to signal that it has the 'start' lock. */
        wait_for(td, i);
    }
    D(" done.\n");

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
        for(i= 0; i < NTHREADS; i++) start(td, i);
        for(i= 0; i < NTHREADS; i++) wait_for(td, i);

        for(i= 0; i < c->length; i++) {
            int j;

            c->entries[i]= 0.0;
            for(j= 0; j < NTHREADS; j++) {
                c->entries[i]+= td[j].c->entries[i];
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

        D("Il= %.3f Iu= %.3f\n", Il, Iu);
        if(Iu - Il < epsilon) break;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < p->length; i++) {
            p->entries[i]*= exp2f_fast(c->entries[i] - Il);
        }
    } 

    return Il;
}

#if 1
/* Sixth attempt. Striping. */
#define STRIP_LEN 8

/* Strips are contiguous blocks of entries, to optimise
   vectorised operation.  A strip is always complete. */
struct strip {
    float entries[STRIP_LEN];
    int32_t rows[STRIP_LEN];
};

/* A stripe is all matrix entries between two rows, blocked
   into strips. */
struct stripe {
    int *sci;
    struct strip *strips;
};

struct striped {
    int stripe_width; /* In strips, the last stripe may have fewer. */
    int nstripes;

    struct stripe *stripes;
};

struct striped *
stripe_mat(csc_mat_t *M, int stripe_width) {
    struct striped *SM;
    int s;
    int *stripe_col_starts;

    D("stripes of %d * %d\n", stripe_width, STRIP_LEN);

    SM= malloc(sizeof(struct striped));
    if(!SM) return NULL;

    SM->stripe_width= stripe_width;

    D("%d rows\n", M->nrow);

    /* How many stripes? (The last may be incomplete) */
    SM->nstripes= ROUND_UP(M->nrow, stripe_width * STRIP_LEN)
                / (stripe_width * STRIP_LEN);
    D("%d stripes\n", SM->nstripes);

    SM->stripes= malloc(SM->nstripes * sizeof(struct stripe));
    if(!SM->stripes) return NULL;

    /* We advance these after building each stripe, so we don't
       have to search for the start of the next. */
    stripe_col_starts= calloc(M->ncol, sizeof(int));
    if(!stripe_col_starts) return NULL;
    {
        int c;

        for(c= 0; c < M->ncol; c++) stripe_col_starts[c]= M->ci[c];
    }

    /* Build the stripes one at a time. */
    for(s= 0; s < SM->nstripes; s++) {
        int c;
        int end_row= (s+1) * stripe_width * STRIP_LEN;
        int nstrips= 0;
        int t; /* Strip number. */

        /* There's an index for each column, even if it's empty. */
        SM->stripes[s].sci= calloc(M->ncol + 1, sizeof(int));
        if(!SM->stripes[s].sci) return NULL;

        /* First pass to count strips. */
        for(c= 0; c < M->ncol; c++) {
            int i;

            //assert(M->rows[stripe_col_starts[c]] >= start_row);

            /* Stay within both the column, and the (row) span
               of the stripe.  */
            for(i= stripe_col_starts[c];
                i < M->ci[c+1] && M->rows[i] < end_row; i++);

            /* How many new strips? */
            nstrips+= ROUND_UP(i - stripe_col_starts[c], STRIP_LEN)
                    / STRIP_LEN;
        }
        D("Stripe %d has %d strips\n", s, nstrips);

        /* Allocate strips. */
        SM->stripes[s].strips= malloc(nstrips * sizeof(struct strip));

        /* Second pass to fill strips. */
        t= 0;
        for(c= 0; c < M->ncol; c++) {
            int i, j;
            int ss= stripe_col_starts[0]; /* Strip start (row). */

            SM->stripes[s].sci[c]= t;

            /* Walk the column up to the stripe end row. */
            for(i= stripe_col_starts[c];
                i < M->ci[c+1] && M->rows[i] < end_row; i++) {

                /* Start a new strip if required. */
                if(i - ss >= STRIP_LEN) { t++; ss= i; }

                /* Copy entry and row number. */
                SM->stripes[s].strips[t].entries[i - ss]= M->entries[i];
                SM->stripes[s].strips[t].rows[i - ss]= M->rows[i];
            }

            /* Fill the last strip with zeroes. */
            for(j= i % STRIP_LEN; j % STRIP_LEN != 0; j++) {
                SM->stripes[s].strips[t].entries[j]= 0;
                SM->stripes[s].strips[t].rows[j]= M->rows[i-1];
            }

            /* Update the column start, for the next stripe. */
            stripe_col_starts[c]= i;
        }
        SM->stripes[s].sci[M->ncol]= t; /* Sentinel. */
    }

    free(stripe_col_starts);

    return SM;
}

void
csc_striped_mult(dv_t *y, dv_t *x, struct striped *M) {
    int s;

    dv_zero(y);

    for(s= 0; s < M->nstripes; s++) {
        int c;
        struct stripe *stripe= &M->stripes[s];

        for(c= 0; c < x->length; c++) {
            int i;

            for(i= stripe->sci[c]; i < stripe->sci[c+1]; i++) {
                struct strip *strip= &stripe->strips[i];
                int j;

                for(j= 0; j < STRIP_LEN; j++) {
                    int row= strip->rows[j];

                    y->entries[c]+= strip->entries[j] * x->entries[row];
                }
            }
        }
    }
}

float
blahut_arimoto_6(csc_mat_t *Q, float epsilon) {
    dv_t *p, *q, *c;
    float Il, Iu;
    int i, col, row;
    float *logQ, *logq;
    struct striped *SQ;

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

    D("Striping matrix...");
    SQ= stripe_mat(Q, 4096 / STRIP_LEN); /* XXX */
    D(" done.\n");

    logq= malloc(q->length * sizeof(float));
    if(!logq) { perror("malloc"); abort(); }

    /* Start with a uniform input distribution. */
    dv_uniform(p, 1.0 / Q->nrow);

#ifdef CAP_BENCH
    iterations= 0;
#endif

    while(1) {
        int s;
#ifdef CAP_BENCH
        iterations++;
#endif

        /* Calculate q - output distribution. */
        csc_striped_mult(q, p, SQ);

        /* Precalculate log(q). */
        for(i= 0; i < q->length; i++) logq[i]= log2f_table(q->entries[i]);

        /* Calculate log(c). */
        dv_zero(c);
        for(s= 0; s < SQ->nstripes; s++) {
            for(col= 0; col < Q->ncol; col++) {
                for(i= SQ->stripes[s].sci[col];
                    i < SQ->stripes[s].sci[col+1]; i++) {
                    int j;

                    for(j= 0; j < STRIP_LEN; j++) {
                        row= SQ->stripes[s].strips[i].rows[j];

                        assert(SQ->stripes[s].strips[i].entries[j] > 0);

                        /* Use the log tables. */
                        c->entries[row]+=
                            SQ->stripes[s].strips[i].entries[j]
                         * (log2f_table(SQ->stripes[s].strips[i].entries[j])
                            - logq[col]);
                    }
                }
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

        D("Il= %.3f Iu= %.3f\n", Il, Iu);
        if(Iu - Il < epsilon) break;

        /* Update input distribution.  Here again, we avoid overflow. */
        for(i= 0; i < p->length; i++) {
            p->entries[i]*= exp2f_fast(c->entries[i] - Il);
        }
    } 

    return Il;
}
#endif
