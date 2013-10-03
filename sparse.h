#ifndef __SPARSE_H
#define __SPARSE_H

#include <stdint.h>
#include <stdio.h>

#define ROUND_UP(x,n) (((x)+(n)-1)-(((x)+(n)-1)%n))
#define ROUND_DOWN(x,n) ((x)-((x)%(n)))

/* Quantum for allocating column space. */
#define ROW_BLOCK 16

/* Block-Sparse-Column histogram. */
typedef struct bsc_hist {
    int end_col;     /* First unallocated column. */
    int end_row;     /* First row unallocated in all columns. */
    int *start_rows; /* First non-empty row by column block. */
    int *end_rows;   /* First unallocated row by column block. */
    int **entries;   /* Entries indexed by block and index, then row. */
    int *row_total;  /* Totals by (allocated) row. */
    int total;       /* Total of all counts. */
    int64_t nalloc;  /* Number of allocated entries. */
    int64_t nnz;     /* Number of non-zero entries. */
} bsc_hist_t;

#define CSC_MAX_STRIDE 7
#define CSC_MAX_SPAN 7
#define CSC_VERSION 1
#define CSC_MAGIC "CSC_MATRIXLE"

#define CSC_M_STRIDERADIX 0x7
#define CSC_F_CFREE       0x8
#define CSC_I_ROWSPAN     4
#define CSC_M_ROWSPAN     (0x7 << CSC_I_ROWSPAN)
#define CSC_M_ALLFLAGS (CSC_M_ROWSPAN | CSC_F_CFREE | CSC_M_STRIDERADIX)

#define STRIDE_OF(M) (1 << ((M)->flags & CSC_M_STRIDERADIX))
#define SPAN_OF(M) (1 << (((M)->flags & CSC_M_ROWSPAN) >> CSC_I_ROWSPAN))

#define CSC_ROW_INVALID 0xff

typedef struct csc_matrix {
    /* We use explicit sizes to allow serialisation. */
    int32_t ver;    /* Version. */
    int32_t flags;  /* Bits 2:0 -> log2(stride length), 0 means no stride. */
    int32_t nrow;   /* Number of Rows */
    int32_t ncol;   /* Number of Columns */
    int64_t nnz;    /* Number of Non-Zero entries */

    /* Only one of the {ci,si} is used, depending on whether there is a
     * stride.  The other is NULL. */
    int32_t *ci;    /* start-of-Column Indicies in entries
                       (size ncol+1) */
    int32_t *si;    /* start-of-Stride Indicies in entries
                       (size (ncol/STRIDE)+1) */

    /* Only used with a stride, otherwise NULL. */
    uint8_t *sc;    /* Column number within the stride (size nnz) */

    int32_t *rows;  /* Row numbers of entries (size nnz, or
                       nnz/stride if collision-free) */

    uint8_t *row_offsets;
                    /* Offset from the row number of the containing block.
                       (size nnz), NULL if not collision-free. */

    float *entries; /* Entry values (size nnz) */
} csc_mat_t;

typedef struct dense_vec {
    int length;
    float *entries;
} dv_t;

typedef enum {
    E_CSC_SUCCESS  = 0,
    E_CSC_TRUNC    = 1,
    E_CSC_BADMAGIC = 2,
    E_CSC_ERRNO    = 3,
    E_CSC_COUNT
} csc_errno_t;

#ifdef __AVX__
#include "sparse_avx.h"
#endif

bsc_hist_t *bsc_hist_new(void);
void bsc_hist_destroy(bsc_hist_t *M);
void bsc_extend(bsc_hist_t *M, int c);
void bsc_extend_col(bsc_hist_t *M, int c, int r);
void bsc_hist_count(bsc_hist_t *M, int c, int r, int n);
csc_mat_t *bsc_normalise(bsc_hist_t *H);
uint64_t bsc_size(bsc_hist_t *H);
void bsc_stats(bsc_hist_t *H);
int bsc_check(bsc_hist_t *H, int verbose);

uint64_t csc_size(csc_mat_t *M);
void csc_mat_destroy(csc_mat_t *M);
csc_errno_t csc_store_binary(csc_mat_t *M, FILE *f);
csc_mat_t *csc_load_binary(FILE *f, csc_errno_t *e);
void mult_csc_dv(dv_t *y, dv_t *x, csc_mat_t *A);
void csc_str_mult_nv(dv_t *y, dv_t *x, csc_mat_t *A);
void csc_stats(csc_mat_t *M);
void csc_prune_cols(csc_mat_t *M);
void csc_stride(csc_mat_t *M, int r_stride);
void csc_perror(csc_errno_t e, const char *s);
int csc_check(csc_mat_t *M, int verbose);
int estimate_cfree(csc_mat_t *M, int row_span);
void csc_make_cfree(csc_mat_t *M, int row_span);
void csc_mult_cf(dv_t *y, dv_t *x, csc_mat_t *A);
void mult_csc_dv_p(dv_t *y, dv_t *x, csc_mat_t *A, int n);
void csc_align(csc_mat_t *M, int n);

dv_t *dv_new(int length);
void dv_destroy(dv_t *v);
void dv_uniform(dv_t *v, float x);
void dv_zero(dv_t *v);
float dv_dot(dv_t *u, dv_t *v);
float dv_max(dv_t *u);

#endif /* __SPARSE_H */
