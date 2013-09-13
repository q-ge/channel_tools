#ifndef __CHANNEL_ALGORITHMS_H
#define __CHANNEL_ALGORITHMS_H

#include "sparse.h"

#ifdef CAP_BENCH
extern int iterations;
#endif

float blahut_arimoto_1(csc_mat_t *Q, float epsilon);
float blahut_arimoto_2(csc_mat_t *Q, float epsilon);
float blahut_arimoto_3(csc_mat_t *Q, float epsilon);
float blahut_arimoto_4(csc_mat_t *Q, float epsilon);
float blahut_arimoto_5(csc_mat_t *Q, float epsilon);
float blahut_arimoto_6(csc_mat_t *Q, float epsilon);

#endif /* __CHANNEL_ALGORITHMS_H */
