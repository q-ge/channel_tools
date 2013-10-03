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
    int cmin= INT_MAX, cmax= 0, ccount;
    int nrow, ncol, rowbin, colskip;
    int xmin= INT_MAX, xmax= 0, ymin= INT_MAX, ymax= 0;
    int c, x, y;
    float *plot;

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

    if(nrow > M->nrow) nrow= M->nrow;

    if(ncol > M->ncol) ncol= M->ncol;

    if(STRIDE_OF(M) != 1) {
        fprintf(stderr, "Too lazy to handle strided matrices\n");
        exit(EXIT_FAILURE);
    }

    for(c= 0; c < M->ncol; c++) {
        if(M->ci[c] < M->ci[c+1]) {
            if(c < cmin) cmin= c;
            if(c > cmax) cmax= c;
        }
    }
    ccount= cmax - cmin + 1;
    if(ncol > ccount) ncol= ccount;

    plot= calloc(nrow * ncol, sizeof(float));

    for(c= 0; c < M->ncol; c++) {
        int i;

        for(i= M->ci[c]; i < M->ci[c+1]; i++) {
            if(M->entries[i] > 0) {
                int r= M->rows[i];
                int x= (r * nrow) / M->nrow;
                int y= ((c-cmin) * ncol) / ccount;

                plot[y * nrow + x]+= M->entries[i];

                if(x < xmin) xmin= x;
                if(x > xmax) xmax= x;
                if(y < ymin) ymin= y;
                if(y > ymax) ymax= y;
            }
        }
    }

    rowbin= M->nrow / nrow;
    colskip= ccount / ncol;
    for(x= xmin; x <= xmax; x++) {
        int r= ((x * rowbin) + ((x+1) * rowbin -1))/2;

        for(y= ymin; y <= ymax; y++) {
            int c= y * colskip + cmin;

            printf("%d %d %.12e\n", r, c, plot[y * nrow + x] / rowbin);
        }
    }

    free(plot);
    csc_mat_destroy(M);

    return 0;
}
