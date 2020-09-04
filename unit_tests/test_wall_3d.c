/**
 * @file test_wall_3d.c
 * @brief Test program for 3D wall collision functions
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../wall/wall_3d.h"

#define N 10000 /**< Number of repetitions in each test */

#define TETRA_N_TRIANGLES 4

/**
 *  Generate a simple 3D wall to test against.
 *
 */
void tetra_wall(wall_3d_offload_data *offload_data,
        real **offload_array) {



	int iTri,iedge;
	const int nTri = TETRA_N_TRIANGLES;

	offload_data->n = nTri;
	offload_data->offload_array_length = nTri*3*3;

    /* The data in the offload array is to be in the format
     *  [x1 y1 z1 x2 y2 z2 x3 y3 z3; ... ]
     */
	*offload_array = (real*) malloc(9 * offload_data->n * sizeof(real));
	iTri=-1;

	iTri++;;iedge=0; /*o,x,z*/
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 1.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 1.0; /* z */

	iTri++;;iedge=0; /*o,x,y*/
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 1.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 1.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */

	iTri++;;iedge=0; /*o,z,y*/
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 1.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 1.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */

	iTri++;;iedge=0; /*x,y,z*/
	(*offload_array)[iTri*9 + iedge*3 + 0] = 1.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 1.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 0.0; /* z */
	iedge++;
	(*offload_array)[iTri*9 + iedge*3 + 0] = 0.0; /* x */
	(*offload_array)[iTri*9 + iedge*3 + 1] = 0.0; /* y */
	(*offload_array)[iTri*9 + iedge*3 + 2] = 1.0; /* z */

}



/**
 * Guess two random points and see if the line between them intersects a wall.
 */
void test_wall_hit(wall_3d_data* wdata) {
    srand(0);
    real Rmin = 2;
    real Rmax = 9;
    real phimin = 0;
    real phimax = 2*3.1415926;
    real zmin = -6;
    real zmax = 6;
    int i;
    for(i = 0; i < N; i++) {
        real q1[3], q2[3];
        q1[0] = ((real)rand()/(real)RAND_MAX)*(Rmax-Rmin) + Rmin;
        q2[0] = ((real)rand()/(real)RAND_MAX)*(Rmax-Rmin) + Rmin;
        q1[1] = ((real)rand()/(real)RAND_MAX)*(phimax-phimin) + phimin;
        q2[1] = ((real)rand()/(real)RAND_MAX)*(phimax-phimin) + phimin;
        q1[2] = ((real)rand()/(real)RAND_MAX)*(zmax-zmin) + zmin;
        q2[2] = ((real)rand()/(real)RAND_MAX)*(zmax-zmin) + zmin;

        int w = wall_3d_hit_wall(q1[0], q1[1], q1[2], q2[0], q2[1], q2[2], wdata);

        printf("%lf %lf %lf %lf %lf %lf %d\n", q1[0], q1[1], q1[2],
               q2[0], q2[1], q2[2], w);
    }
}


void test_collisions(wall_3d_data wdata, real* offload_array) {
    real xmin = 2;
    real xmax = 9;
    real ymin = -1;
    real ymax = 1;
    real zmin = -6;
    real zmax = 6;
    int i;
    for(i = 0; i < N; i++) {
        real q1[3], q2[3];
        q1[0] = ((real)rand()/(real)RAND_MAX)*(xmax-xmin) + xmin;
        q2[0] = ((real)rand()/(real)RAND_MAX)*(xmax-xmin) + xmin;
        q1[1] = ((real)rand()/(real)RAND_MAX)*(ymax-ymin) + ymin;
        q2[1] = ((real)rand()/(real)RAND_MAX)*(ymax-ymin) + ymin;
        q1[2] = ((real)rand()/(real)RAND_MAX)*(zmax-zmin) + zmin;
        q2[2] = ((real)rand()/(real)RAND_MAX)*(zmax-zmin) + zmin;

        real w;
        int j;
        for(j = 0; j < wdata.n; j++) {
            w = wall_3d_tri_collision(q1, q2, &offload_array[9*j],
                    &offload_array[9*j+3], &offload_array[9*j+6]);
            if(w >= 0) {
                break;
            }
        }

        if(j == wdata.n) {
            j = -1;
        }
        printf("%lf %lf %lf %lf %lf %lf %lf %d\n", q1[0], q1[1], q1[2],
               q2[0], q2[1], q2[2], w, j);
    }
}

void test_tree(wall_3d_data* wdata, real* offload_array) {
    int i;
    int ncell = wdata->ngrid*wdata->ngrid*wdata->ngrid;
    printf("%d\n", ncell);
    for(i = 0; i < ncell; i++) {
        int ntris = wdata->tree_array[wdata->tree_array[i]];
        printf("%d\n", ntris);
        int j;
        for(j = 0; j < ntris; j++) {
            printf(" %d\n", wdata->tree_array[wdata->tree_array[i]+1+j]);
        }
    }
}

void test_tri_in_cube(void) {
    real t1[3] = {-1,-1,0.5};
    real t2[3] = { 2, 2,0.5};
    real t3[3] = { 2,-1,0.5};
    real bb1[3] = {0,0,0};
    real bb2[3] = {1,1,1};

    printf("%d\n", wall_3d_tri_in_cube(t1, t2, t3, bb1, bb2));
}

int test_all_tris(real *offload_array, real q1[3], real q2[3], int *shouldhit){
	real w;
    int iTri = 0;
    const int nTri =TETRA_N_TRIANGLES;

    printf("Ray (%f,%f,%f) ---> (%f,%f,%f)\n",
            q1[0],q1[1],q1[2],
            q2[0],q2[1],q2[2]);
    for(iTri = 0; iTri < nTri; iTri++){
        w = wall_3d_tri_collision(q1, q2, &offload_array[9*iTri],
                &offload_array[9*iTri+3], &offload_array[9*iTri+6]);
        if(  (w >= 0.0 && w<= 1.0)  ){
            printf("Tri %d          intersects. ", iTri);
            if ( shouldhit[iTri] ){
                printf("ok!\n");
            }else{
                printf("fail!\n");
                return 1;
            }
        } else {
            printf("Tri %d does not intersect.  ", iTri);
            if ( shouldhit[iTri] ){
                printf("fail!\n");
                return 1;
            }else{
                printf("ok!\n");
            }

        }

    }

	return 0;
}


int test_fixed_rays(real *offload_array ){

	printf("Triangles are: 0:o,x,z 1:o,x,y 2:o,z,y 3:x,y,z\n");

	int ray=-1;

    real q1[3], q2[3];
    const int nTri =TETRA_N_TRIANGLES;
    int shouldhit[nTri]; /* For each tetra triangle */

    ray=0; printf("Ray %d\n",ray);
    q1[0]=0.1; q1[1]=0.1; q1[2]=0.1;
    q2[0]=7.1; q2[1]=8.1; q2[2]=9.1;
    shouldhit[0]=0; shouldhit[1]=0; shouldhit[2]=0; shouldhit[3]=1;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;

    ray++; printf("Ray %d\n",ray);
    q1[0]=0.1; q1[1]=0.1; q1[2]=0.1;
    q2[0]=1.1; q2[1]=0.1; q2[2]=0.1;
    shouldhit[0]=0; shouldhit[1]=0; shouldhit[2]=0; shouldhit[3]=1;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;

    ray=2; printf("Ray %d\n",ray);
    q1[0]=0.1;  q1[1]=0.1; q1[2]=0.1;
    q2[0]=-1.1; q2[1]=0.1; q2[2]=0.1;
    shouldhit[0]=0; shouldhit[1]=0; shouldhit[2]=1; shouldhit[3]=0;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;

    ray++; printf("Ray %d\n",ray);
    q1[0]=-0.1; q1[1]=0.1; q1[2]=0.1;
    q2[0]=-1.1; q2[1]=0.1; q2[2]=0.1;
    shouldhit[0]=0; shouldhit[1]=0; shouldhit[2]=0; shouldhit[3]=0;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;

    ray++; printf("Ray %d\n",ray);
    q1[0]= 0.0; q1[1]=0.0; q1[2]=0.1;
    q2[0]= 0.0; q2[1]=0.0; q2[2]=0.2;
    shouldhit[0]=1; shouldhit[1]=0; shouldhit[2]=1; shouldhit[3]=0;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;

    ray=4; printf("Ray %d\n",ray);
    q1[0]= 0.1; q1[1]=0.0; q1[2]=0.1;
    q2[0]= 0.1; q2[1]=0.0; q2[2]=1.2;
    shouldhit[0]=1; shouldhit[1]=0; shouldhit[2]=0; shouldhit[3]=1;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;

    ray=5; printf("Ray %d\n",ray);
    q1[0]= 0.0; q1[1]=0.0; q1[2]=1.1;
    q2[0]= 0.0; q2[1]=0.0; q2[2]=1.2;
    shouldhit[0]=0; shouldhit[1]=0; shouldhit[2]=0; shouldhit[3]=0;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;

    ray++; printf("Ray %d\n",ray);
    q1[0]= 0.5; q1[1]=0.0; q1[2]=0.5;
    q2[0]= 0.0; q2[1]=0.5; q2[2]=0.5;
    shouldhit[0]=1; shouldhit[1]=0; shouldhit[2]=1; shouldhit[3]=1;
    if(test_all_tris(offload_array, q1, q2, shouldhit))
        return 1;


    return 0;
}



int main(int argc, char** argv) {
    wall_3d_offload_data offload_data;
    real* offload_array;

    /* Get a sample wall */
    tetra_wall(&offload_data, &offload_array);

    wall_3d_init_offload(&offload_data, &offload_array);

    wall_3d_data wdata;
    wall_3d_init(&wdata, &offload_data, offload_array);

    //test_wall_hit(&wdata);
    //test_collisions(wdata, offload_array);
    //test_tree(&wdata, offload_array);
    //test_tri_in_cube();

    if (test_fixed_rays(offload_array)) {
        return 1;
    }

    return 0;
}
