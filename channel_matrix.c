#include <stdio.h>

#include "sparse.h"

#define BLOCK_SIZE 16

int
main(int argc, char *argv[]) {
    bsc_hist_t *H;
    csc_mat_t *M;
    csc_errno_t e;
    int c, r;
    FILE *out;

    if(argc < 2) {
        printf("Usage: %s <output_filename>\n", argv[0]);
        return 1;
    }

    out= fopen(argv[1], "wb");
    if(!out) {
        perror("fopen");
        return 1;
    }

    printf("Building histogram...");
    fflush(stdout);
    H= bsc_hist_new();
    while(scanf("%d %d\n", &r, &c) == 2) bsc_hist_count(H, c, r, 1);
    printf(" done.\n");
    bsc_stats(H);

    fprintf(stderr, "Building matrix\n");
    M= bsc_normalise(H);
    csc_stats(M);
    bsc_hist_destroy(H);

    fprintf(stderr, "Writing matrix\n");
    e= csc_store_binary(M, out);
    if(e != E_CSC_SUCCESS) {
        csc_perror(e, "csc_write_binary");
        return 1;
    }

    return 0;
}
