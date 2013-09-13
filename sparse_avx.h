#ifndef __SPARSE_AVX_H
#define __SPARSE_AVX_H

void mult_csc_dv_avx(dv_t *y, dv_t *x, csc_mat_t *A, int start, int end);
void csc_str_mult_nv_4(dv_t *y, dv_t *x, csc_mat_t *A);
void csc_str_mult_nv_8(dv_t *y, dv_t *x, csc_mat_t *A);
void csc_mult_cf_4_4(dv_t *y, dv_t *x, csc_mat_t *A);

#endif /* __SPARSE_AVX_H */
