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
    int nrow, ncol, colbin;
    int rmin, rmax, rcount;
    int r, c;
    float *bins;
    int *row_bin_numbers;

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

    rmax= 0;
    rmin= INT_MAX;

    for(c= 0; c < M->ncol; c++) {
        int is= M->ci[c], ie= M->ci[c+1];

        if(is < ie) {
            if(M->rows[is] < rmin) rmin= M->rows[is];
            if(M->rows[ie-1] > rmax) rmax= M->rows[ie-1];
        }
    }
    rcount= rmax - rmin + 1;

    if(rcount < nrow) nrow= rcount;

    printf("# %d rows %d cols\n", nrow, ncol);

    colbin= M->ncol / ncol;

    bins= malloc(nrow * sizeof(float));
    row_bin_numbers= malloc(nrow * sizeof(int));
    if(!bins || !row_bin_numbers) { perror("malloc"); exit(EXIT_FAILURE); }

    for(r= 0; r < rcount; r++) {
        int index= (r * nrow) / rcount;
        row_bin_numbers[index]= r;
    }

    for(c= 0; c < M->ncol; c+=colbin) {
        int d;
        float col_bin_prob= 0.0;

        bzero(bins, nrow * sizeof(float));

        for(d= 0; d < colbin && c + d < M->ncol; d++) {
            int i;

            for(i= M->ci[c+d]; i < M->ci[c+d+1]; i++) {
                int index;
                r= M->rows[i];
                index= ((r-rmin) * nrow) / rcount;
                bins[index]+= M->entries[i];
                col_bin_prob+= M->entries[i];
            }
        }

        for(r= 0; r < nrow; r++) {
            printf("%d %d %f\n", c + colbin/2, row_bin_numbers[r],
                    bins[r] / col_bin_prob);
        }
    }

    free(bins);
    csc_mat_destroy(M);

    return 0;
}
