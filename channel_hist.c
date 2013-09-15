#include <assert.h>
#include <limits.h>
#include <stdio.h>

#include "sparse.h"

#define BLOCK_SIZE 16

int
main(int argc, char *argv[]) {
    bsc_hist_t *H;
    int c, r;
    FILE *out;
    int ne_cols= 0;
    int min_row= INT_MAX;
    int max_row= 0;

    if(argc < 2) {
        printf("Usage: %s <output_filename>\n", argv[0]);
        return 1;
    }

    out= fopen(argv[1], "w");
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

    for(c= 0; c < H->end_col; c++) {
        if(H->start_rows[c] < H->end_rows[c]) ne_cols++;
        if(H->start_rows[c] < H->end_rows[c]) {
            if(H->start_rows[c] < min_row)
                min_row= H->start_rows[c];
            if(max_row < H->end_rows[c])
                max_row= H->end_rows[c];
        }
    }
    assert(min_row <= max_row);

    fprintf(out, "%d %d %d\n", ne_cols, min_row, max_row);

    for(c= 0; c < H->end_col; c++) {
        for(r= H->start_rows[c]; r < H->end_rows[c]; r++) {
            int ri= r - H->start_rows[c];
            if(H->entries[c][ri] > 0) {
                fprintf(out, "%d %d %d\n", c, r, H->entries[c][ri]);
            }
        }
    }

    fclose(out);

    return 0;
}

