import pytest
import ctypes
import numpy as np

from a5py.libascot import LIBASCOT, init_fun


# typedef struct Octree
# {
#     struct Octree *n000; /**< The [xmin, ymin, zmin] of child nodes [m].      */
#     struct Octree *n100; /**< The [xmax, ymin, zmin] of child nodes [m].      */
#     struct Octree *n010; /**< The [xmin, ymax, zmin] of child nodes [m].      */
#     struct Octree *n110; /**< The [xmax, ymax, zmin] of child nodes [m].      */
#     struct Octree *n001; /**< The [xmin, ymin, zmax] of child nodes [m].      */
#     struct Octree *n101; /**< The [xmax, ymin, zmax] of child nodes [m].      */
#     struct Octree *n011; /**< The [xmin, ymax, zmax] of child nodes [m].      */
#     struct Octree *n111; /**< The [xmax, ymax, zmax] of child nodes [m].      */
#     float bb1[3];        /**< Bounding box xyz minimum limit [m].             */
#     float bb2[3];        /**< Bounding box xyz maximum limit [m].             */
#     list_int_node *list; /**< Linked list for storing triangle IDs.           */
# } Octree;

# Octree_create(
#     Octree **node, float x_min, float x_max, float y_min, float y_max,
#     float z_min, float z_max, size_t depth);

# void Octree_free(Octree **node);

# Octree_add(Octree *node, float t1[3], float t2[3], float t3[3], size_t id);

init_fun(
    "Octree_tri_collision",
    *(2*[ctypes.POINTER(ctypes.c_double)]),
    *(3*[ctypes.POINTER(ctypes.c_float)]),
    restype=ctypes.c_float,
)

init_fun(
    "Octree_tri_in_cube",
    *(5*[ctypes.POINTER(ctypes.c_float)]),
    restype=ctypes.c_int,
)


@pytest.mark.parametrize(
    "p1, p2, v1, v2, k_expected",
    [
        ([0.1, 0.1], [0.1, 0.1], [0.0, 1.0], [1.0, 0.0], 0.5),
        ([0.1, 0.1], [0.1, 0.1], [1.0, 0.0], [0.0, 1.0], 0.5),
        ([-0.1, -0.1], [-0.1, -0.1], [0.0, 1.0], [1.0, 0.0], -1.0),
    ],
)
@pytest.mark.parametrize(
    "triangle_axes, ray_axis",
    [
        ((0, 1), 2),
        ((0, 2), 1),
        ((1, 2), 0),
    ],
)
def test_octree_tri_collision(
    p1, p2, v1, v2, k_expected, triangle_axes, ray_axis
):
    q1 = np.zeros(3, dtype=np.float64)
    q2 = np.zeros(3, dtype=np.float64)

    q1[ray_axis-2] = p1[0]
    q1[ray_axis-1] = p1[1]
    q2[ray_axis-2] = p2[0]
    q2[ray_axis-1] = p2[1]
    q1[ray_axis] = 1.0
    q2[ray_axis] = -1.0

    v1_3d = np.zeros(3, dtype=np.float32)
    v2_3d = np.zeros(3, dtype=np.float32)
    v3 = np.zeros(3, dtype=np.float32)

    a, b = triangle_axes

    v1_3d[a], v1_3d[b] = v1
    v2_3d[a], v2_3d[b] = v2

    k = LIBASCOT.Octree_tri_collision(
        q1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        q2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        v1_3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        v2_3d.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        v3.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    )
    assert k == k_expected

@pytest.mark.parametrize(
    "t1, t2, t3, expected",
    [
        ((0.1, 0.1, 0.1), (0.2, 0.2, 0.2), (0.3, 0.3, 0.3), True),
        ((-0.1, -0.1, -0.1), (-0.2, -0.2, -0.2), (-0.3, -0.3, -0.3), False),
        ((0.1, 0.1, 0.1), (-0.2, -0.2, -0.2), (-0.3, -0.3, -0.3), True),
        ((-10.5, 0.5, 0.5), (10.5, 0.5, 0.5), (0.5, 10.5, 0.5), True),
        ((-0.5, 10., 0.5), (-0.5, -10., 0.5), (100., 0.5, 0.5), True),
    ],
)
def test_octree_tri_in_cube(t1, t2, t3, expected):
    t1 = np.array(t1, dtype=np.float32)
    t2 = np.array(t2, dtype=np.float32)
    t3 = np.array(t3, dtype=np.float32)
    bb1 = np.array([0., 0., 0.], dtype=np.float32)
    bb2 = np.array([1., 1., 1.], dtype=np.float32)
    isin = LIBASCOT.Octree_tri_in_cube(
        t1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        t2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        t3.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        bb1.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        bb2.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    )
    assert isin == expected
