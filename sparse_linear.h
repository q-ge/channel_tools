#ifndef __SPARSE_LINEAR_H
#define __SPARSE_LINEAR_H

#include "sparse.h"

void lin_mult_4(dv_t *y, dv_t *x, csc_mat_t *A);
void linearise(csc_mat_t *M);

#endif /* __SPARSE_LINEAR_H */
