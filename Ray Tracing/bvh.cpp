#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // Part 2, Task 1:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  BVHNode *node = new BVHNode(bbox);

  if (prims.size() <= max_leaf_size) {
    node->prims = new vector<Primitive *>(prims);
    return node;
  }
  vector<Primitive *> right;
  vector<Primitive *> left;
  
  int index = 0;
  if (bbox.extent[0] <= bbox.extent[1])
      index = 1;
  if (bbox.extent[index] <= bbox.extent[2])
      index = 2;
  
  for (Primitive *p : prims) {
      Vector3D mid = p->get_bbox().centroid();
      if (mid[index] < centroid_box.centroid()[index])
          right.push_back(p);
      else
          left.push_back(p);
  }
  vector<Primitive *> temp;
  int size;
  if (left.empty()) {
      size = right.size();
      for (int i=0; i<size; i++) {
          if (i < size/2) {
              left.push_back(right[i]);
              continue;
          }
          temp.push_back(left[i]);
      }
      right = temp;
  }
  if (right.empty()) {
      size = left.size();
      for (int i=0; i<size; i++) {
          if (i < size/2) {
              right.push_back(left[i]);
              continue;
          }
          temp.push_back(right[i]);
      }
      left = temp;
  }
  //cout << "Going Right" << endl;
  node->r = construct_bvh(right, max_leaf_size);
  //cout << "Going Left" << endl;
  node->l = construct_bvh(left, max_leaf_size);
  //cout << "AAAAND we're back" << endl;
  return node;
}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  // Part 2, task 3: replace this.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.

  double x, y;
  bool hit = false;
  //cout << "Start intersect" << endl;
  if (!(node->bb).intersect(ray, x, y))
      return false;
  if ((x <= ray.min_t) || (x >= ray.max_t))
      return false;
  //cout << "Mid inter " << endl;
  if (node->isLeaf()) {
      for (Primitive *p : *(node->prims)) {
          total_isects++;
          if (p->intersect(ray))
              return true;
      }
  }
  else {
      if (intersect(ray, node->r))
          hit = true;
      if (intersect(ray, node->l))
          hit = true;
  }
  //cout << "Hit 111 = " << hit << endl;
  return hit;

}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {
    
  // Part 2, task 3: replace this
//cout << "Start1" << endl;
  if (!intersect(ray, node))
      return false;
  
  //cout << "Start2" << endl;
  if (node->isLeaf()) {
      for (Primitive *p : *(node->prims)) {
          total_isects++;
          if (p->intersect(ray, i))
              return true;
      }
  }
  //cout << "MID2" << endl;
  bool hit = false;
  if (intersect(ray, i, node->r))
      hit = true;
  if (intersect(ray, i, node->l))
      hit = true;
  //cout << "Hit2 = " << hit << endl;
  return hit;
}

}  // namespace StaticScene
}  // namespace CGL
