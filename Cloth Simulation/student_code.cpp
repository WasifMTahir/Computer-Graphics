#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{
  void BezierCurve::evaluateStep()
  {
    // TODO Part 1.
    // Perform one step of the Bezier curve's evaluation at t using de Casteljau's algorithm for subdivision.
    // Store all of the intermediate control points into the 2D vector evaluatedLevels.
      vector<Vector2D> step;
      int recent = evaluatedLevels.size() - 1;
      int points = evaluatedLevels[recent].size();
      for (int i=1; i<points; i++) {
          Vector2D temp1;
          Vector2D temp2;
          temp1 = evaluatedLevels[recent][i];
          temp2 = evaluatedLevels[recent][i-1];
          temp1.x = t*temp1.x + (1-t)*temp2.x;
          temp1.y = t*temp1.y + (1-t)*temp2.y;
          step.push_back(temp1);
      }
      if (step.size())
          evaluatedLevels.push_back(step);
  }


  Vector3D BezierPatch::evaluate(double u, double v) const
  {
    // TODO Part 2.
    // Evaluate the Bezier surface at parameters (u, v) through 2D de Casteljau subdivision.
    // (i.e. Unlike Part 1 where we performed one subdivision level per call to evaluateStep, this function
    // should apply de Casteljau's algorithm until it computes the final, evaluated point on the surface)
      vector<Vector3D> step;
      for (int i=0; i<4; i++) {
          vector<Vector3D> cPoints;
          for (int j=0; j<4; j++)
              cPoints.push_back(controlPoints[i][j]);
          step.push_back(evaluate1D(cPoints, u));
      }
      return evaluate1D(step, v);
  }

  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> points, double t) const
  {
    // TODO Part 2.
    // Optional helper function that you might find useful to implement as an abstraction when implementing BezierPatch::evaluate.
    // Given an array of 4 points that lie on a single curve, evaluates the Bezier curve at parameter t using 1D de Casteljau subdivision.
      int numPoints = points.size();
      if (!numPoints)
          return Vector3D();
      if (numPoints == 1)
          return points[0];
      
      vector<Vector3D> step;
      for (int i=1; i<numPoints; i++) {
          Vector3D temp;
          temp.x = t*points[i].x + (1-t)*points[i-1].x;
          temp.y = t*points[i].y + (1-t)*points[i-1].y;
          temp.z = t*points[i].z + (1-t)*points[i-1].z;
          step.push_back(temp);
      }
      return evaluate1D(step, t);
  }



  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // TODO Returns an approximate unit normal at this vertex, computed by
    // TODO taking the area-weighted average of the normals of neighboring
    // TODO triangles, then normalizing.
    Vector3D n;
    HalfedgeCIter h = halfedge();
    n = cross(h->vertex()->position - h->twin()->vertex()->position, h->next()->vertex()->position - h->next()->twin()->vertex()->position);
    h = h->twin()->next();
    while(true) {
        n = n + cross(h->vertex()->position - h->twin()->vertex()->position, h->next()->vertex()->position - h->next()->twin()->vertex()->position);
        h = h->twin()->next();
        if (h == halfedge())
            return n.unit();
    }
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // TODO This method should flip the given edge and return an iterator to the flipped edge.
    
    if (e0->isBoundary())
        return e0;
    EdgeIter e1 = e0->halfedge()->next()->edge();
    EdgeIter e2 = e0->halfedge()->next()->next()->edge();
    EdgeIter e3 = e0->halfedge()->twin()->next()->edge();
    EdgeIter e4 = e0->halfedge()->twin()->next()->next()->edge();
    
    VertexIter v1 = e0->halfedge()->vertex();
    VertexIter v2 = e0->halfedge()->twin()->vertex();
    VertexIter v3 = e0->halfedge()->next()->next()->vertex();
    VertexIter v4 = e0->halfedge()->twin()->next()->next()->vertex();
    
    HalfedgeIter h1 = e0->halfedge();
    HalfedgeIter h2 = e0->halfedge()->next();
    HalfedgeIter h3 = e0->halfedge()->next()->next();
    HalfedgeIter h4 = e0->halfedge()->twin();
    HalfedgeIter h5 = e0->halfedge()->twin()->next();
    HalfedgeIter h6 = e0->halfedge()->twin()->next()->next();
    HalfedgeIter h7 = e0->halfedge()->next()->twin();
    HalfedgeIter h8 = e0->halfedge()->next()->next()->twin();
    HalfedgeIter h9 = e0->halfedge()->twin()->next()->twin();
    HalfedgeIter h10 = e0->halfedge()->twin()->next()->next()->twin();
    
    FaceIter f1 = e0->halfedge()->face();
    FaceIter f2 = e0->halfedge()->twin()->face();
    FaceIter f3 = h7->face();
    FaceIter f4 = h8->face();
    FaceIter f5 = h9->face();
    FaceIter f6 = h10->face();
    
    h1->setNeighbors(h2, h4, v4, e0, f1);
    h2->setNeighbors(h3, h8, v3, e2, f1);
    h3->setNeighbors(h1, h9, v1, e3, f1);
    h4->setNeighbors(h5, h1, v3, e0, f2);
    h5->setNeighbors(h6, h10, v4, e4, f2);
    h6->setNeighbors(h4, h7, v2, e1, f2);
    h7->setNeighbors(h7->next(), h6, v3, e1, f3);
    h8->setNeighbors(h8->next(), h2, v1, e2, f4);
    h9->setNeighbors(h9->next(), h3, v4, e3, f5);
    h10->setNeighbors(h10->next(), h5, v2, e4, f6);
    
    v1->halfedge() = h3;
    v2->halfedge() = h6;
    v3->halfedge() = h4;
    v4->halfedge() = h1;
    
    e0->halfedge() = h1;
    e1->halfedge() = h6;
    e2->halfedge() = h2;
    e3->halfedge() = h3;
    e4->halfedge() = h5;
    
    f1->halfedge() = h1;
    f2->halfedge() = h4;
    
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // TODO This method should split the given edge and return an iterator to the newly inserted vertex.
    // TODO The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    
    if (e0->isBoundary())
        return VertexIter();

    EdgeIter e1 = e0->halfedge()->next()->edge();
    EdgeIter e2 = e0->halfedge()->next()->next()->edge();
    EdgeIter e3 = e0->halfedge()->twin()->next()->edge();
    EdgeIter e4 = e0->halfedge()->twin()->next()->next()->edge();
    EdgeIter e5 = newEdge();
    EdgeIter e6 = newEdge();
    EdgeIter e7 = newEdge();
    e5->isNew = true;
    e6->isNew = true;
    e7->isNew = true;
    
    VertexIter v1 = e0->halfedge()->vertex();
    VertexIter v2 = e0->halfedge()->twin()->vertex();
    VertexIter v3 = e0->halfedge()->next()->next()->vertex();
    VertexIter v4 = e0->halfedge()->twin()->next()->next()->vertex();
    VertexIter v5 = newVertex();
    v5->isNew = true;
    v5->position = (v1->position + v2->position) / 2;
    
    HalfedgeIter h1 = e0->halfedge();
    HalfedgeIter h2 = e0->halfedge()->next();
    HalfedgeIter h3 = e0->halfedge()->next()->next();
    HalfedgeIter h4 = e0->halfedge()->twin();
    HalfedgeIter h5 = e0->halfedge()->twin()->next();
    HalfedgeIter h6 = e0->halfedge()->twin()->next()->next();
    HalfedgeIter h7 = e0->halfedge()->next()->twin();
    HalfedgeIter h8 = e0->halfedge()->next()->next()->twin();
    HalfedgeIter h9 = e0->halfedge()->twin()->next()->twin();
    HalfedgeIter h10 = e0->halfedge()->twin()->next()->next()->twin();
    HalfedgeIter h11 = newHalfedge();
    HalfedgeIter h12 = newHalfedge();
    HalfedgeIter h13 = newHalfedge();
    HalfedgeIter h14 = newHalfedge();
    HalfedgeIter h15 = newHalfedge();
    HalfedgeIter h16 = newHalfedge();
    
    FaceIter f1 = e0->halfedge()->face();
    FaceIter f2 = e0->halfedge()->twin()->face();
    FaceIter f3 = h7->face();
    FaceIter f4 = h8->face();
    FaceIter f5 = h9->face();
    FaceIter f6 = h10->face();
    FaceIter f7 = newFace();
    FaceIter f8 = newFace();
    
    h1->setNeighbors(h14, h4, v1, e0, f1);
    h2->setNeighbors(h13, h7, v2, e1, f7);
    h3->setNeighbors(h1, h8, v3, e2, f1);
    h4->setNeighbors(h5, h1, v5, e0, f2);
    h5->setNeighbors(h16, h9, v1, e3, f2);
    h6->setNeighbors(h12, h10, v4, e4, f8);
    h7->setNeighbors(h7->next(), h2, v3, e1, f3);
    h8->setNeighbors(h8->next(), h3, v1, e2, f4);
    h9->setNeighbors(h9->next(), h5, v4, e3, f5);
    h10->setNeighbors(h10->next(), h6, v2, e4, f6);
    h11->setNeighbors(h2, h12, v5, e5, f7);
    h12->setNeighbors(h15, h11, v2, e5, f8);
    h13->setNeighbors(h11, h14, v3, e6, f7);
    h14->setNeighbors(h3, h13, v5, e6, f1);
    h15->setNeighbors(h6, h16, v5, e7, f8);
    h16->setNeighbors(h4, h15, v4, e7, f2);
    
    v1->halfedge() = h1;
    v2->halfedge() = h12;
    v3->halfedge() = h2;
    v4->halfedge() = h16;
    v5->halfedge() = h4;
    
    e0->halfedge() = h4;
    e1->halfedge() = h2;
    e2->halfedge() = h3;
    e3->halfedge() = h5;
    e4->halfedge() = h6;
    e5->halfedge() = h12;
    e6->halfedge() = h13;
    e7->halfedge() = h15;

    f1->halfedge() = h3;
    f2->halfedge() = h5;
    f7->halfedge() = h2;
    f8->halfedge() = h6;
    
    return v5;
    
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    VertexIter begin = mesh.verticesBegin();
    VertexIter end = mesh.verticesEnd();
    EdgeIter beginEdge = mesh.edgesBegin();
    EdgeIter endEdge = mesh.edgesEnd();
    EdgeIter last;

    for (VertexIter vert = begin; vert != end; vert++) {
        vert->isNew = false;
        int n = 0;
        double u = 0;
        HalfedgeIter half = vert->halfedge();
        Vector3D neighbor_position_sum = half->twin()->vertex()->position;
        half = half->twin()->next();
        while (half != vert->halfedge()) {
            neighbor_position_sum = neighbor_position_sum + half->twin()->vertex()->position;
            half = half->twin()->next();
            n++;
        }
        u = 3/(8*n);
        vert->newPosition = (1 - n*u) * vert->position + u * neighbor_position_sum;
    }
    for (EdgeIter edge = beginEdge; edge != endEdge; edge++) {
        Vector3D A = edge->halfedge()->vertex()->position;
        Vector3D B = edge->halfedge()->twin()->vertex()->position;
        Vector3D C = edge->halfedge()->next()->next()->vertex()->position;
        Vector3D D = edge->halfedge()->twin()->next()->next()->vertex()->position;
        edge->newPosition = 3.0/8.0 * (A+B) + 1.0/8.0 * (C+D);
        last = edge;
        edge->isNew = false;
    } 
    
    for (EdgeIter e = beginEdge; e != endEdge; e++) {
        last = e;
    }
    for (EdgeIter edge = beginEdge; edge != last; edge++) {
        mesh.splitEdge(edge)->halfedge()->vertex()->newPosition = edge->newPosition;
    }
    for (EdgeIter edge = beginEdge; edge != endEdge; edge++) {
        if (edge->isNew) {
            bool a = edge->halfedge()->vertex()->isNew;
            bool b = edge->halfedge()->twin()->vertex()->isNew;
            if ((a | b) & !(a & b))
                mesh.flipEdge(edge);
        }
    }
    for (VertexIter vert = begin; vert != end; vert++) {
        vert->position = vert->newPosition;
    }
}  
}
