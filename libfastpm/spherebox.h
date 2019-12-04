/* Sphere Box intersection.
 *
 * Retrieved from http://theorangeduck.com/page/correct-box-sphere-intersection
 *
 * rainwoodman:
 * The original blog post is listed mostly as is. Except the following changes in the source code:
 *  - double -> real
 *  - bool -> int
 *  - Added the missing type definitions.
 *  - box type also remembers its corners.
 *  - Add box_in_sphere (uses the corners).
 *  - rename sphere_intersects_box to sphere_intersects_or_contains_box.
 *
 * Text starting with ">" are quoted from the original blog post.
 *
 * -----------------------------
 * > Created on April 2, 2013, 1:21 p.m.
 *
 * > Edit [31/03/2019]:
 * > The code presented here (embarrassingly given the name of this article I
 * know) before an edit in 2019 actually contained an error and produced false
 * positives in some cases where the sphere could be closely nestled against
 * the box but not completely intersecting - in this case it would say an
 * intersection was taking place when it was actually not.
 * The code you see here has since been updated to fix this error.  I guess the
 * lesson is don't trust code you find on the internet, in particular if it has
 * been written by me!
 *
 * > Testing for intersection between a box and sphere is not as trivial as it may first seem.
 *
 * > I found many resources on the web explaining how to do it,
 * but after struggling for a while I realized all of them listed an incorrect algorithm.
 * What they stated was that you had to try for intersection against all face planes of the box.
 * Some did mentioned briefly that this is flawed, but most didn't seem to even be aware.
 * With a little bit of thought the reason this is incorrect becomes obvious.
 *
 * > Box Sphere Intersection
 *
 * > If you simply check for intersection with the face planes of the box you
 * get a bunch of false positives where the planes extend past the ends of the
 * box. While perhaps a few false positives may not be a big deal in some
 * cases, when there are large spheres and acceleration structures such as quad
 * trees, false positives can mean many, many more operations than usual.
 * Actually doing the correct checks are not much more expensive and overall
 * improve performance.

 * > To start we need to write three tests for checking if a sphere is inside,
 * outside or intersecting a plane. This can be done by taking the signed
 * distance from the plane and comparing to the sphere radius.

 * > If the distance is negative and greater than the radius we know it is
 * inside. If the distance is positive and greater than the radius we know it
 * is outside. Finally if the absolute distance is less than or equal to the
 * radius we know there is an intersection.
 */

#include <math.h>

#ifndef real
#define real double
#endif

typedef union {
    struct {
        real x;
        real y;
        real z;
    };
    real v[3];
} vec3;

typedef struct {
    vec3 position;
    vec3 direction;
} plane;

typedef struct {
    union {
    struct {
    plane front;
    plane back;
    plane top;
    plane bottom;
    plane right;
    plane left;
    };
    plane planes[6];
    };
    vec3 corners[8];
} box;

typedef struct {
    vec3 center;
    real radius;
} sphere;

box make_box(double xmin[3], double xmax[3])
{
    box bbox = {0};

    double *p = xmin;
    double *q = xmax;
    for(int d = 0; d < 3; d ++ ) {
        bbox.back.position.v[d] = p[d];
        bbox.left.position.v[d] = p[d];
        bbox.bottom.position.v[d] = p[d];

        bbox.front.position.v[d] = q[d];
        bbox.right.position.v[d] = q[d];
        bbox.top.position.v[d] = q[d];
    }

    /* Pointing outward. Other components are zeroed at var init.*/
    bbox.back.direction.x = -1;
    bbox.left.direction.y = -1;
    bbox.bottom.direction.z = -1;

    bbox.front.direction.x = 1;
    bbox.right.direction.y = 1;
    bbox.top.direction.z = 1;
    for(int d = 0; d < 3; d ++) {
        bbox.corners[0].v[d] = p[d];
        bbox.corners[1].v[d] = p[d];
        bbox.corners[2].v[d] = p[d];
        bbox.corners[3].v[d] = p[d];

        bbox.corners[4].v[d] = q[d];
        bbox.corners[5].v[d] = q[d];
        bbox.corners[6].v[d] = q[d];
        bbox.corners[7].v[d] = q[d];
    }
    bbox.corners[1].x = q[0];
    bbox.corners[2].y = q[1];
    bbox.corners[3].z = q[2];

    bbox.corners[5].x = p[0];
    bbox.corners[6].y = p[1];
    bbox.corners[7].z = p[2];
    return bbox;
}

sphere make_sphere0(double r)
{
    sphere s;
    for(int d = 0; d < 3; d ++) {
        s.center.v[d] = 0;
    }
    s.radius = r;
    return s;
}

real vec3_dist(vec3 v1, vec3 v2) {
  return sqrt((v1.x - v2.x) * (v1.x - v2.x) +
         (v1.y - v2.y) * (v1.y - v2.y) +
         (v1.z - v2.z) * (v1.z - v2.z));
}

real vec3_dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 vec3_sub(vec3 a, vec3 b) {
    return (vec3){.x=a.x - b.x, .y=a.y - b.y, .z=a.z - b.z};
}

vec3 vec3_mul(vec3 a, real f) {
    return (vec3){.x=a.x * f, .y=a.y *f, .z=a.z *f};
}

real plane_distance(plane p, vec3 point) {
  return vec3_dot(vec3_sub(point, p.position), p.direction);
}

int point_inside_sphere(sphere s, vec3 point) {
  return vec3_dist(s.center, point) < s.radius;
}

int sphere_inside_plane(sphere s, plane p) {
  return -plane_distance(p, s.center) > s.radius;
}

int sphere_outside_plane(sphere s, plane p) {
  return plane_distance(p, s.center) > s.radius;
}

int sphere_intersects_plane(sphere s, plane p) {
  return fabs(plane_distance(p, s.center)) <= s.radius;
}

/* > Using these we can test for if a sphere is inside, outside or intersecting
 * a box (specified as 6 planes). The reasoning behind this representation of 6
 * planes is that we can use it with non axis-aligned boxes, and boxes of odd
 * shapes such as view frustums (which was my original application). When
 * constructing these planes it is important to ensure all the plane normals
 * face outward.
 *
 * rainwoodman: the algebra definitions of inside, intersect and outside are
 *   - inside : s * b = s
 *   - intersects_or_contains : 0 != s * b < s  (this includes b inside s)
 *   - outside : s * b = 0
 *
 * >The test for if a sphere is inside a box is easy. Just check that the sphere
 * is inside all of the face planes.
 */

int sphere_inside_box(sphere s, box b) {

  if (!sphere_inside_plane(s, b.front))  { return 0; }
  if (!sphere_inside_plane(s, b.back))   { return 0; }
  if (!sphere_inside_plane(s, b.top))    { return 0; }
  if (!sphere_inside_plane(s, b.bottom)) { return 0; }
  if (!sphere_inside_plane(s, b.left))   { return 0; }
  if (!sphere_inside_plane(s, b.right))  { return 0; }

  return 1;

}

/* > The intersection test, our original goal, is a little harder, and for this
 * we will need to know not just if a sphere is inside a plane, but also the
 * closest point on the place where the sphere is intersecting, and also the
 * radius of the circle which is produced by the intersection.
 */
int sphere_intersects_plane_point(sphere s, plane p, vec3* point, real* radius) {
  real d = plane_distance(p, s.center);
  vec3 proj = vec3_mul(p.direction, d);
  *point = vec3_sub(s.center, proj);
  *radius = sqrt(fmax(s.radius * s.radius - d * d, 0));
  return fabs(d) <= s.radius; 
}

/* > Then what we need to do is this: for each face check if there is an
 * intersection, and if there is an intersection then ensure that the distance
 * between the closest point on the intersecting plane is within the radius of
 * the circle of intersection.
 */

/* rainwoodman: this returns true if the box is inside the sphere;
 * so not exactly "intersects" box; renamed the function to avoid
 * confusion.
 */
int sphere_intersects_or_contains_box(sphere s, box b) {

  vec3 point;
  real radius;

  if (sphere_intersects_plane_point(s, b.top, &point, &radius)) {
   if (plane_distance(b.left,  point) <= radius &&
        plane_distance(b.right, point) <= radius &&
        plane_distance(b.front, point) <= radius &&
        plane_distance(b.back,  point) <= radius) {
      return 1;
    }

  }

  if (sphere_intersects_plane_point(s, b.bottom, &point, &radius)) {

    if (plane_distance(b.left,  point) <= radius &&
        plane_distance(b.right, point) <= radius &&
        plane_distance(b.front, point) <= radius &&
        plane_distance(b.back,  point) <= radius) {
      return 1;
    }

  }

  if (sphere_intersects_plane_point(s, b.left, &point, &radius)) {

    if (plane_distance(b.top,    point) <= radius &&
        plane_distance(b.bottom, point) <= radius &&
        plane_distance(b.front,  point) <= radius &&
        plane_distance(b.back,   point) <= radius) {
      return 1;
    }

  }

  if (sphere_intersects_plane_point(s, b.right, &point, &radius)) {

    if (plane_distance(b.top,    point) <= radius &&
        plane_distance(b.bottom, point) <= radius &&
        plane_distance(b.front,  point) <= radius &&
        plane_distance(b.back,   point) <= radius) {
      return 1;
    }

  }

  if (sphere_intersects_plane_point(s, b.front, &point, &radius)) {

    if (plane_distance(b.top,    point) <= radius &&
        plane_distance(b.bottom, point) <= radius &&
        plane_distance(b.left,   point) <= radius &&
        plane_distance(b.right,  point) <= radius) {
      return 1;
    }

  }

  if (sphere_intersects_plane_point(s, b.back, &point, &radius)) {

    if (plane_distance(b.top,    point) <= radius &&
        plane_distance(b.bottom, point) <= radius &&
        plane_distance(b.left,   point) <= radius &&
        plane_distance(b.right,  point) <= radius) {
      return 1;
    }

  }

  return 0;

}
/* > Once we have these tests, checking if a box is outside is as simple as
 * ensuring it is not inside or intersecting.
 */
int sphere_outside_box(sphere s, box b) {
  return !(sphere_inside_box(s, b) || sphere_intersects_or_contains_box(s, b));
}

/* > It may seem more expensive but really it is only a handful more dot
 * products and conditionals which the compiler can hopefully optimize away
 * nicely. I feel like capturing the logic and making it clear is more
 * important. Hope this helps any of you who were in my situation!
 */

/* rainwoodman:
 * We actually also need to test if a box is inside a sphere.
 *
 **/
int box_inside_sphere(box b, sphere s) {
  for (int i = 0; i < 8; i ++) {
    if (!point_inside_sphere(s, b.corners[i])) return 0;
  }
  return 1;
}
