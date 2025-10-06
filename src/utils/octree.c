/**
 * Implements octree.h.
 */
#include "octree.h"
#include "defines.h"
#include "list.h"
#include "mathlib.h"
#include <stdlib.h>

void octree_create(
    octree_node **node, real x_min, real x_max, real y_min, real y_max,
    real z_min, real z_max, int depth)
{
    octree_node *n = (octree_node *)malloc(sizeof(octree_node));
    (*node) = n;
    real epsilon = 1e-6;
    n->bb1[0] = x_min - epsilon;
    n->bb2[0] = x_max + epsilon;
    n->bb1[1] = y_min - epsilon;
    n->bb2[1] = y_max + epsilon;
    n->bb1[2] = z_min - epsilon;
    n->bb2[2] = z_max + epsilon;
    real x = (x_max + x_min) / 2;
    real y = (y_max + y_min) / 2;
    real z = (z_max + z_min) / 2;
    if (depth > 1)
    {
        octree_create(&(n->n000), x_min, x, y_min, y, z_min, z, depth - 1);
        octree_create(&(n->n100), x, x_max, y_min, y, z_min, z, depth - 1);
        octree_create(&(n->n010), x_min, x, y, y_max, z_min, z, depth - 1);
        octree_create(&(n->n110), x, x_max, y, y_max, z_min, z, depth - 1);
        octree_create(&(n->n001), x_min, x, y_min, y, z, z_max, depth - 1);
        octree_create(&(n->n101), x, x_max, y_min, y, z, z_max, depth - 1);
        octree_create(&(n->n011), x_min, x, y, y_max, z, z_max, depth - 1);
        octree_create(&(n->n111), x, x_max, y, y_max, z, z_max, depth - 1);
        n->list = NULL;
    }
    else
    {
        n->n000 = NULL;
        n->n100 = NULL;
        n->n010 = NULL;
        n->n110 = NULL;
        n->n001 = NULL;
        n->n101 = NULL;
        n->n011 = NULL;
        n->n111 = NULL;
        list_int_create(&(n->list));
    }
}

void octree_free(octree_node **node)
{
    if ((*node)->n000 != NULL)
    {
        octree_free(&((*node)->n000));
        octree_free(&((*node)->n100));
        octree_free(&((*node)->n010));
        octree_free(&((*node)->n110));
        octree_free(&((*node)->n001));
        octree_free(&((*node)->n101));
        octree_free(&((*node)->n011));
        octree_free(&((*node)->n111));
    }
    else
    {
        list_int_free(&((*node)->list));
    }
    free(*node);
    *node = NULL;
}

void octree_add(octree_node *node, real t1[3], real t2[3], real t3[3], int id)
{
    if (node->n000 == NULL)
    {
        /* leaf node, add the triangle to the list */
        list_int_add(node->list, id);
    }
    else
    {
        if (octree_tri_in_cube(
                t1, t2, t3, node->n000->bb1, node->n000->bb2) > 0)
            octree_add(node->n000, t1, t2, t3, id);
        if (octree_tri_in_cube(
                t1, t2, t3, node->n100->bb1, node->n100->bb2) > 0)
            octree_add(node->n100, t1, t2, t3, id);
        if (octree_tri_in_cube(
                t1, t2, t3, node->n010->bb1, node->n010->bb2) > 0)
            octree_add(node->n010, t1, t2, t3, id);
        if (octree_tri_in_cube(
                t1, t2, t3, node->n110->bb1, node->n110->bb2) > 0)
            octree_add(node->n110, t1, t2, t3, id);
        if (octree_tri_in_cube(
                t1, t2, t3, node->n001->bb1, node->n001->bb2) > 0)
            octree_add(node->n001, t1, t2, t3, id);
        if (octree_tri_in_cube(
                t1, t2, t3, node->n101->bb1, node->n101->bb2) > 0)
            octree_add(node->n101, t1, t2, t3, id);
        if (octree_tri_in_cube(
                t1, t2, t3, node->n011->bb1, node->n011->bb2) > 0)
            octree_add(node->n011, t1, t2, t3, id);
        if (octree_tri_in_cube(
                t1, t2, t3, node->n111->bb1, node->n111->bb2) > 0)
            octree_add(node->n111, t1, t2, t3, id);
    }
}


list_int_node *octree_get(octree_node *node, real p[3])
{
    if (node->n000 == NULL)
    {
        return node->list;
    }
    else
    {
        real x = (node->bb1[0] + node->bb2[0]) / 2;
        real y = (node->bb1[1] + node->bb2[1]) / 2;
        real z = (node->bb1[2] + node->bb2[2]) / 2;

        if (p[0] < x && p[1] < y && p[2] < z)
            return octree_get(node->n000, p);
        if (p[0] > x && p[1] < y && p[2] < z)
            return octree_get(node->n100, p);
        if (p[0] < x && p[1] > y && p[2] < z)
            return octree_get(node->n010, p);
        if (p[0] > x && p[1] > y && p[2] < z)
            return octree_get(node->n110, p);
        if (p[0] < x && p[1] < y && p[2] > z)
            return octree_get(node->n001, p);
        if (p[0] > x && p[1] < y && p[2] > z)
            return octree_get(node->n101, p);
        if (p[0] < x && p[1] > y && p[2] > z)
            return octree_get(node->n011, p);
        if (p[0] > x && p[1] > y && p[2] > z)
            return octree_get(node->n111, p);
    }

    /* We should not ever get here but this
       supresses compiler warnings. */
    return NULL;
}

int octree_tri_in_cube(
    real t1[3], real t2[3], real t3[3], real bb1[3], real bb2[3]
) {
    /* check if any point is inside the cube */
    if(bb1[0] <= t1[0] && t1[0] <= bb2[0]
        && bb1[1] <= t1[1] && t1[1] <= bb2[1]
        && bb1[2] <= t1[2] && t1[2] <= bb2[2])
        return 1;
    if(bb1[0] < t2[0] && t2[0] <= bb2[0]
        && bb1[1] <= t2[1] && t2[1] <= bb2[1]
        && bb1[2] <= t2[2] && t2[2] <= bb2[2])
        return 1;
    if(bb1[0] <= t3[0] && t3[0] <= bb2[0]
        && bb1[1] <= t3[1] && t3[1] <= bb2[1]
        && bb1[2] <= t3[2] && t3[2] <= bb2[2])
        return 1;

    /* no such luck; check if any of the cube edges intersects the triangle */
    real c000[3]; c000[0] = bb1[0]; c000[1] = bb1[1]; c000[2] = bb1[2];
    real c100[3]; c100[0] = bb2[0]; c100[1] = bb1[1]; c100[2] = bb1[2];
    real c010[3]; c010[0] = bb1[0]; c010[1] = bb2[1]; c010[2] = bb1[2];
    real c110[3]; c110[0] = bb2[0]; c110[1] = bb2[1]; c110[2] = bb1[2];
    real c001[3]; c001[0] = bb1[0]; c001[1] = bb1[1]; c001[2] = bb2[2];
    real c101[3]; c101[0] = bb2[0]; c101[1] = bb1[1]; c101[2] = bb2[2];
    real c011[3]; c011[0] = bb1[0]; c011[1] = bb2[1]; c011[2] = bb2[2];
    real c111[3]; c111[0] = bb2[0]; c111[1] = bb2[1]; c111[2] = bb2[2];

    if(   octree_tri_collision(c000, c100, t1, t2, t3) >= 0
       || octree_tri_collision(c000, c010, t1, t2, t3) >= 0
       || octree_tri_collision(c110, c010, t1, t2, t3) >= 0
       || octree_tri_collision(c110, c100, t1, t2, t3) >= 0
       || octree_tri_collision(c000, c001, t1, t2, t3) >= 0
       || octree_tri_collision(c010, c011, t1, t2, t3) >= 0
       || octree_tri_collision(c100, c101, t1, t2, t3) >= 0
       || octree_tri_collision(c110, c111, t1, t2, t3) >= 0
       || octree_tri_collision(c001, c101, t1, t2, t3) >= 0
       || octree_tri_collision(c001, c011, t1, t2, t3) >= 0
       || octree_tri_collision(c111, c011, t1, t2, t3) >= 0
       || octree_tri_collision(c111, c101, t1, t2, t3) >= 0)
        return 1;

    /* check for triangle edges intersecting cube quads */
    if(   octree_quad_collision(t1, t2, c000, c100, c110, c010) == 1
       || octree_quad_collision(t1, t2, c000, c010, c011, c001) == 1
       || octree_quad_collision(t1, t2, c000, c100, c101, c001) == 1
       || octree_quad_collision(t1, t2, c010, c110, c111, c011) == 1
       || octree_quad_collision(t1, t2, c100, c101, c111, c110) == 1
       || octree_quad_collision(t1, t2, c001, c101, c111, c011) == 1)
        return 1;
    if(   octree_quad_collision(t3, t2, c000, c100, c110, c010) == 1
       || octree_quad_collision(t3, t2, c000, c010, c011, c001) == 1
       || octree_quad_collision(t3, t2, c000, c100, c101, c001) == 1
       || octree_quad_collision(t3, t2, c010, c110, c111, c011) == 1
       || octree_quad_collision(t3, t2, c100, c101, c111, c110) == 1
       || octree_quad_collision(t3, t2, c001, c101, c111, c011) == 1)
        return 1;
    if(   octree_quad_collision(t1, t3, c000, c100, c110, c010) == 1
       || octree_quad_collision(t1, t3, c000, c010, c011, c001) == 1
       || octree_quad_collision(t1, t3, c000, c100, c101, c001) == 1
       || octree_quad_collision(t1, t3, c010, c110, c111, c011) == 1
       || octree_quad_collision(t1, t3, c100, c101, c111, c110) == 1
       || octree_quad_collision(t1, t3, c001, c101, c111, c011) == 1)
        return 1;
    return 0;
}


int octree_quad_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3], real t4[3]
) {
    if(octree_tri_collision(q1, q2, t1, t2, t3) >= 0
       || octree_tri_collision(q1, q2, t1, t3, t4) >= 0)
        return 1;
    else
        return 0;
}


double octree_tri_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3]
) {
    real q12[3], Q12[3];
    Q12[0] = q2[0] - q1[0];
    Q12[1] = q2[1] - q1[1];
    Q12[2] = q2[2] - q1[2];
    math_unit(Q12, q12);

    real edge12[3];
    edge12[0] = t2[0] - t1[0];
    edge12[1] = t2[1] - t1[1];
    edge12[2] = t2[2] - t1[2];

    real edge13[3];
    edge13[0] = t3[0] - t1[0];
    edge13[1] = t3[1] - t1[1];
    edge13[2] = t3[2] - t1[2];

    real h[3];
    math_cross(q12, edge13, h);
    real det = math_dot(h, edge12);

    /* Check that the triangle has non-zero area */
    real normal[3], area;
    math_cross(edge12, edge13, normal);
    area = math_norm(normal);

    real w = -1.0;
    if( area > WALL_EPSILON ) {
        /* If ray is parallel to the triangle, nudge it a little bit so we don't
           have to handle annoying special cases */
        if( fabs(det) < WALL_EPSILON ) {
            Q12[0] = q2[0] - q1[0] + 2 * WALL_EPSILON * normal[0] / area;
            Q12[1] = q2[1] - q1[1] + 2 * WALL_EPSILON * normal[1] / area;
            Q12[2] = q2[2] - q1[2] + 2 * WALL_EPSILON * normal[2] / area;
            math_unit(Q12, q12);
            math_cross(q12, edge13, h);
            det = math_dot(h, edge12);
        }

        real tq11[3];
        tq11[0] = q1[0] - t1[0];
        tq11[1] = q1[1] - t1[1];
        tq11[2] = q1[2] - t1[2];

        real n[3];
        math_cross(tq11, edge12, n);

        real u = math_dot(h, tq11) / det;
        real v = math_dot(q12, n) / det;

        if( ( u >= 0.0 && u <= 1.0 ) && ( v >= 0.0 && u + v <= 1.0 )  ) {
            w = ( math_dot(n, edge13) / det ) / math_norm(Q12);
            if( w > 1.0 ) {
                w = -1.0;
            }
        }
    }

    return w;

}
