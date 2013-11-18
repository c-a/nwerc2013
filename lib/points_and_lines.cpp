struct point;
struct vec {
  double x, y;
  vec(double _x, double _y): x(_x), y(_y) {}
  vec(const point& a, const point& b); // definition after point
  vec scale(double s) const {
    return vec(x*s, y*s);
  }
  double dot(const vec& b) const{
    return x*b.x + y*b.y;
  }
  double cross(const vec& b) {
    return x*b.y-y*b.x;
  }
  double norm_sq() const {
    return (x*x+y*y);
  }
};

struct point {
  double x, y;
  point(double _x, double _y): x(_x), y(_y) {}
  point(const point& other) {
    x = other.x; y = other.y;
  }
  bool operator<(const point& other) const {
    if (fabs(x-other.x) > EPS) return x < other.x; else return y < other.y;
  }
  bool operator==(const point& other) const {
    return (fabs(x-other.x) < EPS && fabs(y-other.y) < EPS);
  }
  point translate(const vec& v) const {
    return point(x + v.x, y + v.y);
  }
  point rotate(double theta) {
    double rad = theta * M_PI / 180;
    return point(x*cos(rad)-y*sin(rad), x*sin(rad)+y*cos(rad));
  }
};
vec::vec(const point& a, const point& b) { // needs to be here since it uses point
  x = b.x-a.x; y = b.y-a.y;
}

double dist(const point& p1, const point& p2) { // Euclidean distance
    return sqrt(vec(p1, p2).norm_sq());
}
double angle(const point& a, const point& o, const point& b) { // returns angle of aob in rad
  vec oa(o, a), ob(o, b);
  return acos(oa.dot(ob)) / sqrt(oa.norm_sq()*ob.norm_sq());
}
// returns true if point r is to the left side of line pq
bool ccw(const point& p, const point& q, const point& r) {
  return vec(p,q).cross(vec(p,r)) > 0;
}
// returns true if point r is to the left side of line pq
bool collinear(const point& p, const point& q, const point& r) {
  return fabs(vec(p,q).cross(vec(p,r))) < EPS;
}

double distToLine(const point& p, const point& a, const point& b, point& c) {
  // formula: c = a + u * ab;
  vec ap(a, p), ab(a, b);
  double u = ap.dot(ab) / ab.norm_sq();
  c = a.translate(ab.scale(u));
  return dist(p, c);
}

// returns the distance from p to the line segment ab defined by two points a and b. The closest point is stored in c.
double distToLineSegment(const point& p, const point& a, const point& b, point& c) {
  vec ap(a, p), ab(a, b);
  double u = ap.dot(ab) / ab.norm_sq();
  if (u < 0.0) { // closer to a
    c = a; return dist(p, a); // Euclidean distance between p and a
  }
  else if (u > 1.0) { // closer to b
    c = b; return dist(p, b); // Euclidean distance between p and b
  }
  return distToLine(p, a, b, c);
}

struct line {
  double a, b, c;
  line(double _a, double _b, double _c): a(_a), b(_b), c(_c) {}
  line(point& p1, point& p2) {
    if (fabs(p1.x-p2.x) < EPS) { // vertical line
      a = 1.0; b = 0.0; c = -p1.x; // default values
    } else {
      a = -(double)(p1.y-p2.y)/(p1.x-p2.x);
      b = 1.0; // IMPORTANT: we fix the value of b to 1.0
      c = -(double)(a*p1.x) - p1.y;
    }
  }
  bool parallel(const line& l) const {
    return (fabs(a-l.a) < EPS) && (fabs(b-l.b) < EPS);
  }
  bool intersect(const line& l, point& p) const {
    if (parallel(l)) return false;
    // solve system of 2 linear algebraic equations with 2 unknowns
    p.x = (l.b*c - b*l.c) / (l.a*b - a*l.b);
    if (fabs(b) > EPS) p.y = -(a*p.x + c); else p.y = -(l.a*p.x + l.c);
    return true;
  }
  bool operator==(const line& other) const {
    return parallel(other) && (fabs(c-other.c) < EPS);
  }
};
