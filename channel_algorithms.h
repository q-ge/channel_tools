#ifndef __CHANNEL_ALGORITHMS_H
#define __CHANNEL_ALGORITHMS_H

#include "sparse.h"

#ifdef CAP_BENCH
extern int iterations;
#endif
float blahut_arimoto(csc_mat_t *Q, float epsilon, float *e_obs);
float blahut_arimoto_precise(csc_mat_t *Q, float epsilon, float *e_obs);
float blahut_arimoto_precise_squeezed(csc_mat_t *Q, float epsilon, float *e_obs);
float blahut_arimoto_ld_squeezed(csc_mat_t *Q, float epsilon, float *e_obs);
float ba_phased(csc_mat_t *Q, float epsilon, float *e_obs);

#endif /* __CHANNEL_ALGORITHMS_H */
