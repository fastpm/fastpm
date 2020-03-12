#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <alloca.h>
#include <mpi.h>
#include <math.h>

#include "../libfastpm/spherebox.h"

/*
[ -0.986404 -1.30458 -1.51315 ] - [ 513.05 513.066 513.21]
[ 0000078.8033 ]: usmesh intersection from 0.9335 to 1.0000 with 8 tiles.
[ 0000078.8033 ]: usmesh: tile 0 bounding box intersects shell r = (0 20.9997)
[ 0000078.9478 ]: usmesh: tile 1 bounding box intersects shell r = (0 20.9997)
[ 0000079.0941 ]: usmesh: tile 2 bounding box intersects shell r = (0 20.9997)
[ 0000079.2404 ]: usmesh: tile 3 bounding box intersects shell r = (0 20.9997)
[ 0000079.3866 ]: usmesh: tile 4 bounding box intersects shell r = (0 20.9997)
[ 0000079.5330 ]: usmesh: tile 5 bounding box intersects shell r = (0 20.9997)
[ 0000079.6792 ]: usmesh: tile 6 bounding box intersects shell r = (0 20.9997)
[ 0000079.8253 ]: usmesh: tile 7 bounding box intersects shell r = (0 20.9997)
*/

void test_sphere_inside_box(double r1, int expected) {
    double xmin[3] = {-4, -4, -4};
    double xmax[3] = { 20, 4, 4};
    box b = make_box(xmin, xmax);
    sphere s = make_sphere0(r1);
    int got = sphere_inside_box(s, b);
    if (got != expected) {
        printf("%s r1 = %g got = %d expected = %d\n", __func__, r1, got, expected);
    }
}

void test_sphere_intersects_or_contains_box(double r1, int expected) {
    double xmin[3] = {-4, -4, -4};
    double xmax[3] = { 20, 4, 4};
    box b = make_box(xmin, xmax);
    sphere s = make_sphere0(r1);
    int got = sphere_intersects_or_contains_box(s, b);
    if (got != expected) {
        printf("%s r1 = %g got = %d expected = %d\n", __func__, r1, got, expected);
    }
}

void test_box_inside_sphere(double r1, int expected) {
    double xmin[3] = {-4, -4, -4};
    double xmax[3] = { 20, 4, 4};
    box b = make_box(xmin, xmax);
    sphere s = make_sphere0(r1);
    int got = box_inside_sphere(b, s);
    if (got != expected) {
        printf("%s r1 = %g got = %d expected = %d\n", __func__, r1, got, expected);
    }
}



int main(int argc, char * argv[]) {

    test_sphere_inside_box(3.9, 1);
    test_sphere_inside_box(4.1, 0);
    test_sphere_inside_box(20, 0);

    test_sphere_intersects_or_contains_box(3.9, 0);
    test_sphere_intersects_or_contains_box(4.1, 1);
    test_sphere_intersects_or_contains_box(30, 1);

    test_box_inside_sphere(30, 1);
    test_box_inside_sphere(4.1, 0);
    test_box_inside_sphere(3.9, 0);

    return 0;
}
