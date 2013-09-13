/* AVX-optimised csc matrix operations. */
#include <assert.h>

#include <immintrin.h>

#include "sparse.h"

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

#if 0
void
mult_csc_dv_avx(dv_t *y, dv_t *x, csc_mat_t *A, int start, int end) {
    int c, i;

    assert(x); assert(A); assert(y);
    assert(STRIDE_OF(A) == 1);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(0 <= start && start <= end && end <= A->ncol);

    for(c= start; c < end; c++) {
        float tmp= 0.0;
        __v8sf tmp_v= {0,0,0,0,0,0,0,0};
        int len= A->ci[c+1] - A->ci[c];
        __v8sf ones= {1,1,1,1,1,1,1,1};
        __v8sf dp;

        for(i= 0; i < ROUND_DOWN(len, 8); i+= 8) {
            int j= A->ci[c] + i;
            __v8sf A_e= __builtin_ia32_loadups256(A->entries + j);
            __v8sf x_e= {0,0,0,0,0,0,0,0};
            __v4sf tmp_x;

            tmp_x[0]= x->entries[A->rows[j]];
            tmp_x[1]= x->entries[A->rows[j+1]];
            tmp_x[2]= x->entries[A->rows[j+2]];
            tmp_x[3]= x->entries[A->rows[j+3]];
            x_e= __builtin_ia32_vinsertf128_ps256(x_e, tmp_x, 0);

            tmp_x[0]= x->entries[A->rows[j+4]];
            tmp_x[1]= x->entries[A->rows[j+5]];
            tmp_x[2]= x->entries[A->rows[j+6]];
            tmp_x[3]= x->entries[A->rows[j+7]];
            x_e= __builtin_ia32_vinsertf128_ps256(x_e, tmp_x, 1);

            tmp_v+= A_e * x_e;
        }
        dp= __builtin_ia32_dpps256(tmp_v, ones, 0xff);
        tmp+= dp[4] + dp[0];

        for(; i < len; i++) {
            int j= A->ci[c] + i;
            tmp+= A->entries[j] * x->entries[A->rows[j]];
        }

        y->entries[c]= tmp;
    }
}
#endif

void
mult_csc_dv_avx(dv_t *y, dv_t *x, csc_mat_t *A, int start, int end) {
    int c, i;
    //float stuff= 0;

    assert(x); assert(A); assert(y);
    assert(STRIDE_OF(A) == 1);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(0 <= start && start <= end && end <= A->ncol);

#if 1
    for(c= start; c < end; c++) {
        __v4sf tmp_v= {0,0,0,0};
        __v4sf ones= {1,1,1,1};
        __v4sf dp;

        for(i= A->ci[c]; i < A->ci[c+1]; i+= 4) {
            __v4sf A_e;
            __v4sf x_e= {0,0,0,0};

#if 0
            __builtin_prefetch(&x->entries[A->rows[i+4]]);
            __builtin_prefetch(&x->entries[A->rows[i+5]]);
            __builtin_prefetch(&x->entries[A->rows[i+6]]);
            __builtin_prefetch(&x->entries[A->rows[i+7]]);
#endif

            asm volatile("movntdqa %1, %0" : "=x"(A_e) : "m"(A->entries[i]));

            x_e[0]= x->entries[A->rows[i+0]];
            x_e[1]= x->entries[A->rows[i+1]];
            x_e[2]= x->entries[A->rows[i+2]];
            x_e[3]= x->entries[A->rows[i+3]];

            tmp_v+= A_e * x_e;
        }
        dp= __builtin_ia32_dpps(tmp_v, ones, 0xff);

        y->entries[c]= dp[0];
    }
#else
    for(c= start; c < ROUND_DOWN(end,4); c+=1) {
        __v4sf tmp_v= {0,0,0,0};
        __v4sf ones= {1,1,1,1};
        //__v4sf sums= {0,0,0,0};
        __v4sf dp;

        for(i= A->ci[c]; i < A->ci[c+1]; i+= 4) {
            //__v4sf A_e;
            __v4sf A_e= __builtin_ia32_loadups(A->entries + i);
            __v4sf x_e= {0,0,0,0};

            //asm volatile("movntdqa %1, %0" : "=x"(A_e) : "m"(A->entries[i]));

            x_e[0]= x->entries[A->rows[i+0]];
            x_e[1]= x->entries[A->rows[i+1]];
            x_e[2]= x->entries[A->rows[i+2]];
            x_e[3]= x->entries[A->rows[i+3]];
#if 0
            x_e[0]= x->entries[i - A->ci[c]];
            x_e[1]= x->entries[i - A->ci[c] + 1];
            x_e[2]= x->entries[i - A->ci[c] + 2];
            x_e[3]= x->entries[i - A->ci[c] + 3];
#endif

            tmp_v+= A_e * x_e;
        }
        dp= __builtin_ia32_dpps(tmp_v, ones, 0xff);

        //if(i%4 == 0) y->entries[c]= dp[0];
        y->entries[c]= dp[0];
        //stuff+= dp[0];
    }
#endif

    //y->entries[0]= stuff;
}

void
csc_str_mult_nv_4(dv_t *y, dv_t *x, csc_mat_t *A) {
    int s;

    assert(x); assert(A); assert(y);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(A->ncol % STRIDE_OF(A) == 0);
    assert(STRIDE_OF(A) == 4);

    for(s= 0; s < A->ncol/4; s++) {
        int i;
        __v4sf acc= {0,0,0,0};

        for(i= A->si[s]; i < ROUND_DOWN(A->si[s+1],2); i+=4) {
            __v4sf tmp1= {0,0,0,0};
            __v4sf tmp2= {0,0,0,0};
            __v4sf tmp3= {0,0,0,0};
            __v4sf tmp4= {0,0,0,0};
            int c1= A->sc[i];
            int c2= A->sc[i+1];
            int c3= A->sc[i+2];
            int c4= A->sc[i+3];

            tmp1[0]= A->entries[i] * x->entries[A->rows[i]];
            tmp2[0]= A->entries[i+1] * x->entries[A->rows[i+1]];
            tmp3[0]= A->entries[i+2] * x->entries[A->rows[i+2]];
            tmp3[0]= A->entries[i+3] * x->entries[A->rows[i+3]];
            tmp1= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp1, c1*4);
            tmp2= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp2, c2*4);
            tmp3= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp3, c3*4);
            tmp4= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp4, c4*4);
            acc+= tmp1 + tmp2 + tmp3 + tmp4;
        }

        for(;i < A->si[s+1]; i++) {
            __v4sf tmp= {0,0,0,0};
            int c= A->sc[i];
            assert(c < 4);

            tmp[0]= A->entries[i] * x->entries[A->rows[i]];
            tmp= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp, c*4);
            acc+= tmp;
        }

        __builtin_ia32_storeups(y->entries + s*4, acc);
    }
}

void
csc_str_mult_nv_8(dv_t *y, dv_t *x, csc_mat_t *A) {
    int s;

    assert(x); assert(A); assert(y);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(A->ncol % STRIDE_OF(A) == 0);
    assert(STRIDE_OF(A) == 8);

    for(s= 0; s < A->ncol/8; s++) {
        int i;
        __v8sf acc= {0,0,0,0,0,0,0,0};

        for(i= A->si[s]; i < A->si[s+1]; i++) {
            __v8sf tmp= {0,0,0,0,0,0,0,0};
            assert(A->sc[i] < 8);

            tmp[A->sc[i]]= A->entries[i] * x->entries[A->rows[i]];
            acc+= tmp;
        }

        __builtin_ia32_storeups256(y->entries + s*8, acc);
    }
}

void
csc_mult_cf_4_4(dv_t *y, dv_t *x, csc_mat_t *A) {
    __v4sf col_sum= {0, 0, 0, 0};
    int s;

    assert(x); assert(A); assert(y);
    assert(A->nrow == x->length);
    assert(A->ncol == y->length);
    assert(A->ncol % STRIDE_OF(A) == 0);
    assert(STRIDE_OF(A) == 4);
    assert(A->flags & CSC_F_CFREE);
    assert(SPAN_OF(A) == 4);

    for(s= 0; s < A->ncol/4; s++) {
        int i;

        assert(A->si[s] % 4 == 0);
        assert(A->si[s+1] % 4 == 0);

        /* Zero the sum. */
        col_sum= __builtin_ia32_xorps(col_sum, col_sum);

        for(i= A->si[s]; i < A->si[s+1]; i+= 4) {
            __v4sf A_entries, x_entries;
            int row_base;
            __v4si row_offsets;

            asm volatile("foo:\n");

            /* Get the first row index. */
            row_base= A->rows[i/4];

            /* Load and unpack the row offsets. */
            row_offsets[0]= *((uint32_t *)(&A->row_offsets[i]));
            row_offsets=
                __builtin_ia32_pmovzxbd128((__v16qi)row_offsets);

            /* Load the x (row-indexed) values. */
            x_entries= __builtin_ia32_loadups(x->entries + row_base);

            /* Load the matrix entries. */
            A_entries= __builtin_ia32_loadups(A->entries + i);

            /* Shuffle the x values. */
            x_entries= __builtin_shuffle(x_entries, row_offsets);

            col_sum+= A_entries * x_entries;

            asm volatile("bar:\n");
        }

        __builtin_ia32_storeups(y->entries + s*4, col_sum);
    }
}
