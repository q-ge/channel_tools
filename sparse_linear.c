#include <assert.h>
#include <stdlib.h>

#include <immintrin.h>

#include "sparse.h"

void
linearise(csc_mat_t *M) {
    int s;
    void *new_entries;
    uint64_t i;

    if(!csc_check(M, 1)) abort();

    assert(M->nrow < 1<<24);
    assert(STRIDE_OF(M) < 1<<8);

    new_entries= calloc(2 * M->nnz, sizeof(float));
    if(!new_entries) { perror("malloc"); exit(EXIT_FAILURE); }

    for(s= 0; s < M->ncol / STRIDE_OF(M); s++) {
        for(i= M->si[s]; i < M->si[s+1]; i++) {
            ((float *)new_entries)[2*i]= M->entries[i];
            ((uint32_t *)new_entries)[2*i + 1]=
                (M->sc[i] & 0xff) | ((M->rows[i] & 0xffffff) << 24);
        }
    }

    free(M->entries);
    M->entries= new_entries;
}

void
lin_mult_4(dv_t *y, dv_t *x, csc_mat_t *A) {
    int s;

    assert(STRIDE_OF(A) == 4);

    for(s= 0; s < A->ncol / STRIDE_OF(A); s++) {
        int i;
        __v4sf acc= {0,0,0,0};
        //int len= A->si[s+1] - A->si[s];

#if 0
        for(i= 0; i < ROUND_DOWN(len,2); i+= 2) {
            uint64_t pack= ((uint64_t *)A->entries)[A->si[s] + i];
            uint64_t pack2= ((uint64_t *)A->entries)[A->si[s] + i+1];
            int r= pack >> 56;
            int c= (pack >> 32) & 0xff;
            int r2= pack2 >> 56;
            int c2= (pack2 >> 32) & 0xff;
            __v4sf A_e= {0,0,0,0};
            __v4sf A_e2= {0,0,0,0};
            __v4sf x_e= {0,0,0,0};
            __v4sf x_e2= {0,0,0,0};
            __v4sf tmp;

            A_e[0]= pack & 0xffffffff;
            A_e2[0]= pack2 & 0xffffffff;
            x_e[0]= x->entries[r];
            x_e2[0]= x->entries[r2];
            tmp= A_e * x_e;
            tmp= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp, c*4);
            acc+= tmp;
            tmp= A_e2 * x_e2;
            tmp= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp, c2*4);
            acc+= tmp;
        }
#endif

        for(i= A->si[s]; i < A->si[s+1]; i++) {
            //uint64_t pack= ((uint64_t *)A->entries)[A->si[s] + i];
            uint64_t pack= ((uint64_t *)A->entries)[i];
            int r= pack >> 56;
            int c= (pack >> 32) & 0xff;
            __v4sf A_e= {0,0,0,0};
            __v4sf x_e= {0,0,0,0};
            __v4sf tmp;

            //A_e[0]= pack & 0xffffffff;
            A_e[0]= pack & 0xffffffff;
            x_e[0]= x->entries[r];
            tmp= A_e * x_e;
            tmp= (__v4sf)__builtin_ia32_pslldi128((__v4si)tmp, c*4);
            acc+= tmp;
        }
        __builtin_ia32_storeups(y->entries + s*4, acc);
    }
}
