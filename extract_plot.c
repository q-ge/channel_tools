#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "sparse.h"

int
main(int argc, char *argv[]) {
    FILE *in;
    csc_mat_t *M;
    csc_errno_t e;
    int nrow, ncol, rowbin, colskip;
    int cmin, cmax, ccount;
    int c;
    float *bins;

    if(argc < 4) {
        fprintf(stderr, "Usage: %s <matrix_filename> <nrow> <ncol>\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }

    in= fopen(argv[1], "rb");
    if(!in) { perror("fopen"); exit(EXIT_FAILURE); }

    nrow= atoi(argv[2]);
    if(0 >= nrow) {
        fprintf(stderr, "nrow <= 0\n");
        exit(EXIT_FAILURE);
    }

    ncol= atoi(argv[3]);
    if(0 >= ncol) {
        fprintf(stderr, "ncol <= 0\n");
        exit(EXIT_FAILURE);
    }

    M= csc_load_binary(in, &e);
    if(!M) { csc_perror(e, "csc_load_binary"); exit(EXIT_FAILURE); }
    fclose(in);

    if(!csc_check(M, 1)) abort();
    printf("# ");
    csc_stats(M);

    if(nrow > M->nrow) {
        fprintf(stderr, "nrow too large\n");
        exit(EXIT_FAILURE);
    }

    if(ncol > M->ncol) {
        fprintf(stderr, "ncol too large\n");
        exit(EXIT_FAILURE);
    }

    if(STRIDE_OF(M) != 1) {
        fprintf(stderr, "Too lazy to handle strided matrices\n");
        exit(EXIT_FAILURE);
    }

    cmax= 0;
    cmin= INT_MAX;

    for(c= 0; c < M->ncol; c++) {
        int is= M->ci[c], ie= M->ci[c+1];

        if(is < ie) {
            if(c < cmin) cmin= c;
            if(c > cmax) cmax= c;
        }
    }
    ccount= cmax - cmin + 1;

    if(ccount < ncol) ncol= ccount;

    printf("# %d rows %d cols\n", nrow, ncol);

    rowbin= M->nrow / nrow;

    colskip= M->ncol / ncol;

    bins= calloc(nrow, sizeof(float));
    if(!bins) {
        perror("calloc");
        exit(EXIT_FAILURE);
    }

    for(c= 0; c < M->ncol; c+=colskip) {
        int i, b;

        bzero(bins, nrow * sizeof(float));

        for(i= M->ci[c]; i < M->ci[c+1]; i++) {
            int r= M->rows[i];
            b= (r * nrow) / M->nrow;
            bins[b]+= M->entries[i];
        }

        for(b= 0; b < nrow; b++) {
            printf("%f %d %.12e\n",
                ((float)(b * rowbin) + ((b+1) * rowbin - 1))/2,
                c, bins[b] / rowbin);
        }
    }

    free(bins);
    csc_mat_destroy(M);

    return 0;
}
