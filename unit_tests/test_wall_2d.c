/**
 * @file test_wall_2d.c
 * @brief Test program for 2D wall collision functions
 */
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../wall/wall_2d.h"

int main(int argc, char** argv) {
    real wall_r[] = {4.7, 6.1, 6.1, 8.4, 8.4, 7.6, 6.1, 4.8, 4.0, 4.0, 4.7,
        4.7};
    real wall_z[] = {-4.25, -4.25, -3.5, -1.1, 1.75, 2.5, 3.6, 3.6, 2.4, -1.5,
        -3.25, -4.25};

    wall_2d_data wdata;
    wdata.n = 12;
    wdata.wall_r = wall_r;
    wdata.wall_z = wall_z;

    srand(0);

    real rmin = 0;
    real rmax = 10;
    real zmin = -6;
    real zmax = 6;

    int i;
    for(i = 0; i < 1000; i++) {
        real r = ((real)rand()/(real)RAND_MAX)*(rmax-rmin) + rmin;
        real z = ((real)rand()/(real)RAND_MAX)*(zmax-zmin) + zmin;
        int inside = wall_2d_inside(r, z, &wdata);
        printf("%lf, %lf, %d\n", r, z, inside);
    }

    return 0;
}
