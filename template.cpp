
typedef double T;
struct pt {
    T x, y;

    pt operator+(pt p) { return {x + p.x, y + p.y}; }

    pt operator-(pt p) { return {x - p.x, y - p.y}; }

    pt operator*(T d) { return {x * d, y * d}; }

    pt operator/(T d) { return {x / d, y / d}; } // only for floatingpoint
    pt translate(pt v, pt p) { return p + v; }

    pt rot(pt p, double a) {
        return {p.x * cos(a) - p.y * sin(a), p.x * sin(a) + p.y * cos(a)};
    }


};

T sq(pt p) { return p.x * p.x + p.y * p.y; }

double abs(pt p) { return sqrt(sq(p)); }

bool operator==(pt a, pt b) { return a.x == b.x && a.y == b.y; }

bool operator!=(pt a, pt b) { return !(a == b); }

T dot(pt v, pt w) { return v.x * w.x + v.y * w.y; }

T cross(pt v, pt w) { return v.x * w.y - v.y * w.x; }

T orient(pt a, pt b, pt c) { return cross(b - a, c - a); }

//left turn,collinear, right turn = orient(A, B, C) > 0, orient(A, B, C) = 0, orient(A, B,
pt perp(pt p) { return {-p.y, p.x}; }

double angle(pt v, pt w) {
    return acos(clamp(dot(v, w) / abs(v) / abs(w), -1.0, 1.0));
}

bool isConvex(vector<pt> p) {
    bool hasPos = false, hasNeg = false;
    for (int i = 0, n = p.size(); i < n; i++) {
        int o = orient(p[i], p[(i + 1) % n], p[(i + 2) % n]);
        if (o > 0) hasPos = true;
        if (o < 0) hasNeg = true;
    }
    return !(hasPos && hasNeg);
}

bool half(pt p) { // true if in blue half
    assert(p.x != 0 || p.y != 0); // the argument of (0,0) isundefined
    return p.y > 0 || (p.y == 0 && p.x < 0);
}

void polarSort(vector<pt> &v) {
    sort(v.begin(), v.end(), [](pt v, pt w) {
        return make_tuple(half(v), 0, sq(v)) <
               make_tuple(half(w), cross(v, w), sq(w));
    });
}

//We can perform a polar sort around some point O other than the
//        origin: we just have to subtract that point O from the vectors #»v and
//#»w when comparing them. This as if we translated the whole plane so
//that O is moved to (0, 0):
//void polarSortAround(p
struct line {
    pt v;
    T c;

// From direction vector v and offset c
    line(pt v, T c) : v(v), c(c) {}

// From equation ax+by=c
    line(T a, T b, T c) : v({b, -a}), c(c) {}

// From points P and Q
    line(pt p, pt q) : v(q - p), c(cross(v, p)) {}

// Will be defined later:
// - these work with T = int
    T side(pt p) { return cross(v, p) - c; }

    double dist(pt p) { return abs(side(p)) / abs(v); }

    line perpThrough(pt p) { return {p, p + perp(v)}; }

    bool cmpProj(pt p, pt q) {
        return dot(v, p) < dot(v, q);
    }

    line translate(pt t) { return {v, c + cross(v, t)}; }

    // - these require T = double

    line shiftLeft(double dist) { return {v, c + dist * abs(v)}; }
//    There is a unique intersection point between two lines l1 and l2 if and only
//    if # »vl1 × # »vl2 != 0
    bool inter(line l1, line l2, pt &out) {
        T d = cross(l1.v, l2.v);
        if (d == 0) return false;
        out = (l2.v*l1.c - l1.v*l2.c) / d; // requires floating-pointcoordinates
        return true;
    }
    pt proj(pt p) {return p - perp(v)*side(p)/sq(v);}
    pt refl(pt p) {return p - perp(v)*2*side(p)/sq(v);}

};
bool inDisk(pt a, pt b, pt p) {
    return dot(a-p, b-p) <= 0;
}
bool onSegment(pt a, pt b, pt p) {
    return orient(a,b,p) == 0 && inDisk(a,b,p);
}
bool properInter(pt a, pt b, pt c, pt d, pt &out) {
    double oa = orient(c,d,a),
            ob = orient(c,d,b),
            oc = orient(a,b,c),
            od = orient(a,b,d);
// Proper intersection exists iff opposite signs
    if (oa*ob < 0 && oc*od < 0) {
        out = (a*ob - b*oa) / (ob-oa);
        return true;
    }
    return false;
}
//    Segment-point distance
double segPoint(pt a, pt b, pt p) {
    if (a != b) {
        line l(a,b);
        if (l.cmpProj(a,p) && l.cmpProj(p,b)) // if closest toprojection
            return l.dist(p); // output distance toline
    }
    return min(abs(p-a), abs(p-b)); // otherwise distance to A or B
}
// Segment-segment distance
double segSeg(pt a, pt b, pt c, pt d) {
    pt dummy;
    if (properInter(a,b,c,d,dummy))
        return 0;
    return min({segPoint(a,b,c), segPoint(a,b,d),
                segPoint(c,d,a), segPoint(c,d,b)});
}
double areaPolygon(vector<pt> p) {
    double area = 0.0;
    for (int i = 0, n = p.size(); i < n; i++) {
        area += cross(p[i], p[(i+1)%n]); // wrap back to 0 if i == n-1
    }
    return abs(area) / 2.0;
}
bool above(pt a, pt p) {
    return p.y >= a.y;
}
// check if [PQ] crosses ray from A
bool crossesRay(pt a, pt p, pt q) {
    return (above(a,q) - above(a,p)) * orient(a,p,q) > 0;
}
bool inPolygon(vector<pt> p, pt a, bool strict = true) {
    int numCrossings = 0;
    for (int i = 0, n = p.size(); i < n; i++) {
        if (onSegment(p[i], p[(i+1)%n], a))
            return !strict;
        numCrossings += crossesRay(a, p[i], p[(i+1)%n]);
    }
    return numCrossings & 1; // inside if odd number of crossings
}
