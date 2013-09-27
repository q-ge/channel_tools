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
    int nrow, ncol, rowskip, colbin;
    int maxrow;
    int r, c;
    float *bins;
    int mask= 0;

    if(argc < 4) {
        fprintf(stderr, "Usage: %s <matrix_filename> <nrow> <ncol> "
                " [<max row>]\n",
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

    if(argc >= 5) maxrow= atoi(argv[4]);
    else maxrow= INT_MAX;
    //if(argc >= 5 && !strncmp(argv[4], "-m", 2)) mask= 1;

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

    if(maxrow > M->nrow) maxrow= M->nrow;

    rowskip= maxrow / nrow;
    if(maxrow % nrow != 0) rowskip++;
    colbin=  M->ncol / ncol;
    if(M->ncol % ncol != 0) colbin++;
    printf("# %d rows %d cols\n", nrow, ncol);

    bins= malloc(nrow * sizeof(float));
    if(!bins) { perror("malloc"); exit(EXIT_FAILURE); }

    for(c= 0; c < M->ncol; c+=colbin) {
        int d;
        float col_bin_prob= 0.0;

        bzero(bins, nrow * sizeof(float));

        for(d= 0; d < colbin && c + d < M->ncol; d++) {
            int i;

            for(i= M->ci[c+d]; i < M->ci[c+d+1]; i++) {
                r= M->rows[i];
                if(r > maxrow) continue;
                assert(0 <= r / rowskip &&
                       r / rowskip < nrow);
                bins[r / rowskip]+= M->entries[i];
                col_bin_prob+= M->entries[i];
            }
        }

        for(r= 0; r < nrow; r++) {
            if(mask) {
                if(bins[r] != 0)
                    printf("%d %d %f\n", c + colbin/2, r * rowskip, 1.0);
                else
                    printf("%d %d %f\n", c + colbin/2, r * rowskip, 0.0);
            }
            else
                printf("%d %d %f\n", c + colbin/2, r * rowskip,
                        bins[r] / col_bin_prob);
        }
    }

    free(bins);
    csc_mat_destroy(M);

    return 0;
}
