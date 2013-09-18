#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "sparse.h"

#define MAXLINE 1024

int
main(int argc, char *argv[]) {
    bsc_hist_t *H;
    csc_mat_t *M;
    csc_errno_t e;
    int c, r;
    FILE *out;
    int cmin= INT_MIN, cmax= INT_MAX;
    int climit= -1;
    int *counts;
    int discard= 0;
    size_t malformed= 0, out_of_range= 0;
    char buf[MAXLINE];

    if(argc != 2 && argc != 4 && argc != 5 && argc != 6) {
        printf("Usage: %s <output_filename> [<col. min> <col. max> "
               "[<count limit> [<discard>]]]\n",
                argv[0]);
        return 1;
    }

    out= fopen(argv[1], "wb");
    if(!out) {
        perror("fopen");
        return 1;
    }

    if(argc >= 4) {
        cmin= atoi(argv[2]);
        cmax= atoi(argv[3]);
    }

    if(argc >= 5) {
        climit= atoi(argv[4]);
        counts= calloc(cmax - cmin + 1, sizeof(int));
        if(!counts) {
            perror("calloc");
            exit(EXIT_FAILURE);
        }
    }

    if(argc >= 6) {
        discard= atoi(argv[5]);
    }

    printf("Building histogram...");
    fflush(stdout);
    H= bsc_hist_new();
    while(fgets(buf, MAXLINE, stdin)) {
        int n;

        n= sscanf(buf, "%d %d\n", &c, &r);

        if(n != 2) {
            malformed++;
            continue;
        }

        if(c < cmin || cmax < c) {
            out_of_range++;
            continue;
        }

        if(climit >= 0) {
            int i= c - cmin;

            if(counts[i] >= climit + discard)
                continue;

            counts[i]++;

            if(counts[i] <= discard)
                continue;
        }

        bsc_hist_count(H, c, r, 1);
    }
    printf(" done.\n");
    bsc_stats(H);

    printf("Building matrix\n");
    M= bsc_normalise(H);
    csc_stats(M);
    bsc_hist_destroy(H);

    printf("Writing matrix\n");
    e= csc_store_binary(M, out);
    if(e != E_CSC_SUCCESS) {
        csc_perror(e, "csc_write_binary");
        return 1;
    }

    fclose(out);

    printf("%lu malformed entries, %lu columns out of range\n",
            malformed, out_of_range);

    return 0;
}
