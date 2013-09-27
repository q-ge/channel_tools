#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char *argv[]) {
    int cmin, cmax, rmin, rmax;
    int c, r;

    if(argc < 5) {
        fprintf(stderr, "Usage: %s <cmin> <cmax> <rmin> <rmax>\n", argv[0]);
        return 1;
    }

    cmin= atoi(argv[1]);
    cmax= atoi(argv[2]);
    rmin= atoi(argv[3]);
    rmax= atoi(argv[4]);

    while(scanf("%d %d\n", &c, &r) == 2) {
        if(cmin <= c && c <= cmax &&
           rmin <= r && r <= rmax) {
            printf("%d %d\n", c, r);
        }
    }

    return 0;
}
