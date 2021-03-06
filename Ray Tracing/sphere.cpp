#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // Part 1, Task 4:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
    double a = dot(r.d, r.d);
    double b = dot(r.d, r.o-o) * 2;
    double c = dot(r.o-o, r.o-o) - r2;
    double disc = b*b - 4*a*c;
    double sqdisc = sqrt(disc)/a;
    double t;
    if (disc>0) {
        t1 = (-b + sqdisc)/(2*a);
        t2 = (-b - sqdisc)/(2*a);
        if (t1>t2) {
            double temp = t1;
            t1 = t2;
            t2 = temp;
        }
        return true;
    }
    if (disc==0) {
        t1 = -b/(2*a);
        t2 = t1;
        return true;
    }
    return false;
}

bool Sphere::intersect(const Ray& r) const {

  // Part 1, Task 4:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
    double t1, t2;
    if (test(r, t1, t2)) { 
        if ((t1 < r.max_t) && (t1 > r.min_t)) {
            return true;
        }
    }
    return false;

}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // Part 1, Task 4:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
    double t1, t2;
    if (test(r, t1, t2)) {
        if ((t1 < r.max_t) && (t1 > r.min_t)) {
            i->t = t1;
            i->primitive = this;
            i->n = normal(r.o + r.d*t1);
            i->bsdf = get_bsdf();
            r.max_t = t1;
            return true;
        }
    }
    return false;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
