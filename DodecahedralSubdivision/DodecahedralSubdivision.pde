import java.util.Map;
import java.util.Collections;
import java.util.Collection;
import java.util.Set;
int SubDivLevel = 1;
float speed = 1.0; 
float time = 0.0;
int NGONs1 = 5; //number of polygon sides
float RNG1 = 300;  //maximum radius of a polygon vertex from polygon center for the initializing polygon
float atime = 1.0; //(animation time in seconds)
float frac = .2; //fractional size (recommend that this is < .5)
PVector[] NG1pos = new PVector[NGONs1];
ArrayList<ArrayList<Polygon>> PFamily = new ArrayList<ArrayList<Polygon>>();
ArrayList<PShape> shapes = new ArrayList<PShape>();
class Polygon{
  ArrayList<PVector> vertices;
  PVector center;
  ArrayList<Edge> arcs; //these are animation arcs for a given subdivision
  ArrayList<Float> arcincs;   // animation scalar increment value
  ArrayList<PVector> spokenorms;  // normal vectors used in computing animation position... unit vectors in the direction of the spoke
  ArrayList<Float> sincvals;  //  current spoke scalar 'distance' values on animation cycle 
  ArrayList<PVector> intpolypoints; //interior polygon points for Pshape rendering...the same as p2s...no difference... deprecated at the moment...
  ArrayList<PVector> p2s;  //points originating on the smallest interior polygon boundaries
  ArrayList<Edge> edges;
  HashMap<ArrayList<PVector>,Integer> pointsToEdge; // A vertices to edges index lookup
  //points are keyed to a given winding order, ArrayList<PVector> is of point paired size.
  HashMap<Integer,Integer> ptToSubdivPt;
  ArrayList<PVector> subdivpts;
  HashMap<ArrayList<Integer>,Boolean> SubdivPtPairtoThick;
  HashMap<Integer,Integer> subdivPtToPt;
  Polygon(ArrayList<PVector> verts, PVector PCent){
    vertices = verts;
    center = PCent;
    arcs = new ArrayList<Edge>();
    arcincs = new ArrayList<Float>();
    spokenorms = new ArrayList<PVector>();
    sincvals = new ArrayList<Float>();

    intpolypoints = new ArrayList<PVector>();
    p2s = new ArrayList<PVector>();
    edges = new ArrayList<Edge>();
    pointsToEdge = new HashMap<ArrayList<PVector>,Integer>();
  }
  Polygon(){
    vertices = new ArrayList<PVector>();
    center = new PVector();
    arcs = new ArrayList<Edge>();
    arcincs = new ArrayList<Float>();
    spokenorms = new ArrayList<PVector>();
    sincvals = new ArrayList<Float>();

    intpolypoints = new ArrayList<PVector>();
    p2s = new ArrayList<PVector>();
    edges = new ArrayList<Edge>();
    pointsToEdge = new HashMap<ArrayList<PVector>,Integer>();
  }
  void subdivPtToPtbuild(){
  //assuming that ptToSubdivPt exists
    subdivPtToPt = new HashMap<Integer,Integer>();
    for (Map.Entry<Integer,Integer> me : ptToSubdivPt.entrySet()) {
      int value = me.getKey();
      int ikey = me.getValue();
      subdivPtToPt.put(ikey,value);
    }
  }
  
  Edge getEdgefromSubdivPt(int i){
    //subdivPtToPt must already be computed
    //i is a subdivpt indexed point.
    Set<Integer> k = subdivPtToPt.keySet();
    int minval = -1;
    for (int l : k){

      boolean t3 = i > l;
      boolean t4 = l > minval;
      if (t3  && t4){
        minval = l;
      }
    }
    int vptindex1 = subdivPtToPt.get(minval);
    ArrayList<Integer> npts = new ArrayList<Integer>();
    npts.add(vptindex1);
    return edges.get(vptindex1);
    
  }
}

class Edge{
  PVector p1;  // p1 is point 1 on the edge pair
  PVector p2;  // p2 is point 2 on the edge pair
  int p1index;  // polygon vertices index
  int p2index;  
  PVector genCenter;  //generator center.  This is the edge coordinate generating center
  // genCenter is the main generating circle for such edge
  float genR; // radius of genrating circle.
  float angle1; //polar angle1 generating p1
  float angle2; //polar angle2 generating p2
  boolean thick; //Thick edge is true otherwise thin edge is false
  boolean linear;  // Linear if true otherwise arc if false
  Circle circle;
  int interiornp1;  //this is only called for main edge to interior edges mappings
  int interiornp2; //this is only called for main edge to interior edges mappings 
  int pole1;  //poles are the original vertices on the polygon prior to subdivision
              //This is a subdivision index mapping.
  int pole2;  //used in special case edge mappings
  boolean pole;  //flagged edge that wraps around a pole
  Edge(PVector P1, PVector P2, boolean Thick, int P1index, int P2index){
    //This is a assumed linear edge constructor call...use this constructor call
    //only if the edge is linear
    p1 = P1;
    p2 = P2;
    linear = true;
    thick = Thick;
    p1index = P1index;
    p2index = P2index;
  }
  Edge(PVector P1, PVector P2, boolean Thick, PVector GenCenter,
       float GenR, float Angle1, float Angle2, int P1index, int P2index){
    //This is a assumed linear edge constructor call...use this constructor call
    //only if the edge is linear
    p1 = P1;
    p2 = P2;
    linear = true;
    thick = Thick;
    genCenter = GenCenter;
    genR = GenR;
    angle1 = Angle1;
    angle2 = Angle2;
    p1index = P1index;
    p2index = P2index;
    circle = new Circle(GenCenter,GenR);
  }
  Edge(PVector P1, PVector P2, PVector GenCenter,
       float GenR, int P1index, int P2index){
    //This is a special edge constructor call used for edge Circle to Interior
    // Circles mapping. That is in forming the subdivisions of an existing
    // polygon.
    p1 = P1;
    p2 = P2;
    linear = false;
    //thick = Thick;
    genCenter = GenCenter;
    genR = GenR;
    //angle1 = Angle1;
    //angle2 = Angle2;
    p1index = P1index;
    p2index = P2index;
    //interiornp1 = Interiornp1;
    //interiornp2 = Interiornp2;
    circle = new Circle(GenCenter,GenR);
  }
  Edge(Circle circlei){
    genCenter = circlei.center;
    genR = circlei.radius;
    circle = circlei;
    angle1 = 0.0;
    angle2 = 2.0*PI;
  }
  Edge(Circle circlei, float Angle1, float Angle2){
    genCenter = circlei.center;
    genR = circlei.radius;
    circle = circlei;
    angle1 = Angle1;
    angle2 = Angle2;
  }
  Edge(){
    p1 = new PVector(0.0,0.0,0.0);
    p2 = new PVector(0.0,0.0,0.0);  // p2 is point 2 on the edge pair
    p1index = 0;  // polygon vertices index
    p2index = 0;  
    genCenter = new PVector(0.0,0.0,0.0);  //generator center.  This is the edge coordinate generating center
  // genCenter is the main generating circle for such edge
    genR = 1.0; // radius of genrating circle.
    angle1 = 0.0; //polar angle1 generating p1
    angle2 = 0.0; //polar angle2 generating p2
    thick = false; //Thick edge is true otherwise thin edge is false
    linear = false;  // Linear if true otherwise arc if false
    circle = new Circle(genCenter,genR);
    interiornp1 = 0;  //this is only called for main edge to interior edges mappings
    interiornp2 = 0; //this is only called for main edge to interior edges mappings 
    pole1 = 0;  //poles are the original vertices on the polygon prior to subdivision
              //This is a subdivision index mapping.
    pole2 = 0;  //used in special case edge mappings
    pole = false;  //flagged edge that wraps around a pole
  }
  void computeAngles(){
    PVector angleGen1 = PVector.sub(p1,genCenter);
    PVector angleGen2 = PVector.sub(p2,genCenter);
    angle1 = angleGen1.heading();
    angle2 = angleGen2.heading();
  }
  void computeRadius(){
    PVector angleGen1 = PVector.sub(p1,genCenter);
    genR = angleGen1.mag();
  }
  
  PVector getEdgeHalfAnglePt(){
    PVector p1pc = PVector.sub(p1, genCenter);
    float anglediff = angle2 - angle1;
    anglediff *= .5;
    p1pc.rotate(anglediff);
    return PVector.add(p1pc,genCenter);
  }
}

class Circle{
  PVector center;
  float radius;
  Circle(PVector Center, float Radius){
    center = Center;
    radius = Radius;
  }
}

class CircleMap{
  //this class is used in constructing subdiv points arc edge to paired interior arc
  //edge mappings...these maps have a winding around the polygon. 
  //
  HashMap<ArrayList<Integer>,Circle> vertPairToCircle;  //vertex pairs (indexing) to circle
  HashMap<Circle,ArrayList<Circle>> interiorCircles;  // circle to interior circle mapping
  HashMap<Circle, ArrayList<Integer>> circleToVertPair;
  HashMap<Edge,ArrayList<Edge>> interiorEdges;  //Edge to interior edges.
  CircleMap(){
    vertPairToCircle = new HashMap<ArrayList<Integer>,Circle>();
    interiorCircles = new HashMap<Circle,ArrayList<Circle>>();
    circleToVertPair = new HashMap<Circle,ArrayList<Integer>>();
    interiorEdges = new HashMap<Edge,ArrayList<Edge>>();
  }
}

class PolyHolding{
  //pre appending class as an intermediary storage container used when 
  //constructing polygon data...for ease in preappending polygon data when this
  //is finally sent to a given build polygon method call.
  HashMap<ArrayList<Integer>,Integer> ptsToEdge;
  ArrayList<PVector> vertices;
  ArrayList<Edge> parentEdges;
  HashMap<Integer,Integer> vertReIndexing;
  PVector[] verticesArray;
  Edge[] parentEdgesArray;
  PolyHolding(){
    ptsToEdge = new HashMap<ArrayList<Integer>,Integer>(); //reindexed mapped pair to current edge index mapping
    vertices = new ArrayList<PVector>();
    parentEdges = new ArrayList<Edge>();
    vertReIndexing = new HashMap<Integer,Integer>();
  }
  void addVertex(PVector pt, int posMap){
    //posMap is the desired winding order map for the constructed polygon
    //assumed clockwise ordering
    vertices.add(pt);
    vertReIndexing.put(vertices.size()-1,posMap);
  }
  void addEdge(Integer pt1ind, Integer pt2ind, Edge edge){
    //pt1 index is the desired winding order map for the constructed polygon
    //assumed clockwise winding.
    ArrayList<Integer> ptpair = new ArrayList();
    ptpair.add(pt1ind);
    ptpair.add(pt2ind);
    parentEdges.add(edge);
    ptsToEdge.put(ptpair,parentEdges.size()-1);
  }
  
  boolean testPointEqWithError(PVector p1, PVector p2){
    float error = .01;
    PVector p3 = PVector.sub(p1,p2);
    boolean t1 = abs(p3.x) <= error;
    boolean t2 = abs(p3.y) <= error;
    if (t1 && t2){
      return true;
    }
    else{
      return false;
    }
  }
  
  void addPointEdgetoPoly(PointEdgetoPoly pointedgetopoly){
    //when appending a point and edge.  The point on edge is assumed leading on the 
    //clockwise winding.  Filling the polygon thus means edges are never paired with 
    // a non leading point on the clockwise winding.
    
    if (!vertices.contains(pointedgetopoly.PolyHPt)){
      boolean t1 = false;
      for (PVector vertex: vertices){
        if (testPointEqWithError(vertex, pointedgetopoly.PolyHPt)){
          t1 = true;
          break;
        }
      }
      if (!t1){
        addVertex(pointedgetopoly.PolyHPt, pointedgetopoly.PolyHPtIndex);
        addEdge(pointedgetopoly.PolyHPtIndex, pointedgetopoly.PolyHPtIndex2, 
                pointedgetopoly.PolyHParentEdge);
      }
    }
  }
  
  void writeArrayData(){
    //This is called after all data has been populated in non array object type members
    verticesArray = new PVector[vertices.size()];
    parentEdgesArray = new Edge[parentEdges.size()];
    //println("PointsToEdge: ", ptsToEdge);
    //println("ParentEdges size: " ,parentEdges.size());
    for (int i = 0; i < vertices.size(); i++){
      int j = vertReIndexing.get(i);
      int nj = (j+1)%vertices.size();
      verticesArray[j] = vertices.get(i);
      ArrayList<Integer> pts = new ArrayList<Integer>();
      pts.add(j);
      pts.add(nj);
      parentEdgesArray[j] = parentEdges.get(ptsToEdge.get(pts));
    }
  }
}

class PointEdgetoPoly{
  //simple point to edger to polygon index write interface for PolyHolding 
  int PolyHPtIndex;
  int PolyHPtIndex2;  //this is the second PolyH pt index for edge pairing.
  PVector PolyHPt;
  Edge PolyHParentEdge; 
  PointEdgetoPoly(int polyhptindex, int polyhptindex2, PVector polyhpt, 
                  Edge polyhparentedge){
    PolyHPtIndex = polyhptindex;
    PolyHPtIndex2 = polyhptindex2;
    PolyHPt = polyhpt;
    PolyHParentEdge = polyhparentedge; 
  }
  PointEdgetoPoly(){
    PolyHPtIndex = 0;
    PolyHPtIndex2 = 0;
    PolyHPt = new PVector(0.0,0.0,0.0);
    PolyHParentEdge = new Edge();
  }
}

void CircleCircleIntersection(Circle c1, Circle c2, ArrayList<PVector> ipts){
  PVector c1c2 = PVector.sub(c1.center, c2.center);
  float c1c2angle = c1c2.heading();
  float dc1c2 = c1c2.mag();
  float x = (dc1c2*dc1c2 - c1.radius*c1.radius + c2.radius*c2.radius)/(2*dc1c2);
  float y2 = (4*pow(dc1c2,2)*pow(c2.radius,2) -pow(pow(dc1c2,2)-pow(c1.radius,2)+pow(c2.radius,2),2))/(4*pow(dc1c2,2));
  //ipoint.x = x;
  //ipoint.y = y;
  float yn = -1*pow(y2,.5);
  float yp = pow(y2,.5);
  PVector xyn = PVector.add((new PVector(x,yn,0.0)).rotate(c1c2angle),c2.center);
  PVector xyp = PVector.add((new PVector(x,yp,0.0)).rotate(c1c2angle),c2.center);
  ipts.add(xyn);
  ipts.add(xyp);
}

void Polarcoord(float angle, float radius, PVector xy){
  xy.x = radius*cos(angle);
  xy.y = radius*sin(angle);
   //changes to xy extend outside of method/function call since it is wrapped in a container that is call by reference like...
}

void getNGonPoints(int nside, float radius, PVector[] pos){
  float ainc = 2.0*PI/float(nside);  //generates equal angle subidivison for Ngon creation
  float a = -PI/2.0;
  if (nside % 2 == 0){  //checks for polygon side eveness versus oddness
    a += ainc/2.0;
  }
  for (int i = 0; i < nside; i++){
    PVector point = new PVector(0.0,0.0,0.0);
    Polarcoord(a, radius, point);
    pos[i] = point;
    a += ainc;
  }
}

float distPointToLine(PVector p0, PVector p1, PVector p2){
  //p1 and p2 represent points on the line, and p3 is a point elsewhere.
  //see also http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  PVector p0p1 = PVector.sub(p0,p1);
  PVector p0p2 = PVector.sub(p0,p2);
  PVector p2p1 = PVector.sub(p2,p1);
  PVector p0p1Cp0p2 = p0p1.cross(p0p2);
  float dp2p1 = p2p1.mag();
  float dp0p1Cp0p2 = p0p1Cp0p2.mag();
  return dp0p1Cp0p2/dp2p1;
}

float distPointToCircle(Circle ccircle, PVector p, PVector out){
  // assumed that position p is set not the origin which is the Circle center,
  // but that both the circle ccircle and p have a common origin.
  // for derivation of optimization of the distance of a point to a circle see
  //https://docs.google.com/document/d/1ZVrjdMJ2ZRRJZIeBEdDLBcCFGD3YR8hW13OlqV1Lvag/edit?usp=sharing
  //p will be set with the circle ccircle's center as origin.
  PVector porigc = PVector.sub(p, ccircle.center);
  //then point on the circle ccircle with minimal distance to porigc is given by
  float r = pow(porigc.x/porigc.y,2);
  float xpos = r*pow(ccircle.radius,2)/pow((1+r),0.5);
  float xneg = -1*xpos;
  float ypos1 = pow(pow(ccircle.radius,2) - pow(xpos,2),0.5);
  float ypos2 = -1*ypos1;
  PVector v1 = new PVector(xpos,ypos1,0.0);
  PVector v2 = new PVector(xpos,ypos2,0.0);
  PVector v3 = new PVector(xneg,ypos1,0.0);
  PVector v4 = new PVector(xneg,ypos2,0.0);
  PVector v1p = PVector.sub(v1,porigc);
  PVector v2p = PVector.sub(v2,porigc);
  PVector v3p = PVector.sub(v3,porigc);
  PVector v4p = PVector.sub(v4,porigc);
  float dv1p = v1p.mag();
  float dv2p = v2p.mag();
  float dv3p = v3p.mag();
  float dv4p = v4p.mag();
  HashMap<Float,PVector> dtoV = new HashMap<Float,PVector>();
  dtoV.put(dv1p,v1);
  dtoV.put(dv2p,v2);
  dtoV.put(dv3p,v3);
  dtoV.put(dv4p,v4);
  float[] darr = {dv1p,dv2p,dv3p,dv4p};
  //get minimum distance on the set of all possible points for optimization
  float mindist = min(darr);
  PVector minvec = dtoV.get(mindist);
  //set vector original origin
  out = PVector.add(minvec, ccircle.center);
  return mindist;
}

float distPointToCircle(Circle ccircle, PVector p){
  PVector pc = PVector.sub(p,ccircle.center);
  pc.normalize();
  pc = PVector.mult(pc,ccircle.radius);
  PVector k = PVector.add(ccircle.center,pc);
  PVector kp = PVector.sub(k,p);
  return kp.mag();
}

void PolygonCentroid(ArrayList<PVector> verts, PVector PCenter){
  float Cx = 0.0;
  float Cy = 0.0;
  float A = 0.0;
  for(int i=0; i<verts.size(); i++){
    int ni = 0;
    if (i == verts.size()-1){
      ni = 0;
    }
    else{
      ni = i+1;
    }
    Cx += (verts.get(i).x + verts.get(ni).x)*(verts.get(i).x*verts.get(ni).y - verts.get(ni).x*verts.get(i).y);
    Cy += (verts.get(i).y + verts.get(ni).y)*(verts.get(i).x*verts.get(ni).y - verts.get(ni).x*verts.get(i).y);
    A += verts.get(i).x*verts.get(ni).y - verts.get(ni).x*verts.get(i).y;
  }
  A*=.5;
  Cx*=1.0/(6.0*A);
  Cy*=1.0/(6.0*A);
  PCenter.x = Cx;
  PCenter.y = Cy;
}
//given 2 points on a circle find the center of the circle
void getCircleCenter(PVector p1, PVector p2, PVector Centerout){
  //no solution for triangulation without specified radius since
  //substitution of slope of the line given by the chord points
  //leads to cancellation of a necessary solving term.
  //
  PVector p1p2 = PVector.sub(p2,p1);
  float dp1p2 = p1p2.mag();
  p1p2.normalize();
  PVector bisectp1p2 = PVector.add(p1, PVector.mult(p1p2, dp1p2*.5));
  float pslope = -1.0*(p1.x - p2.x)/(p1.y - p2.y);  //slope of line perpendicular to 
  // p1p2 running through the center of the circle from the bisector of p1p2
  //y = bisectp1p2.y + pslope*(x - bisectp1p2.x)
  //distance from center (a,b) to p1 and p2 is the same so
  //we have the following formula reducing the distance formula for center (a,b) and p1 and p2
  //p1.x^2 -p2.x^2 + p1.y^2 -p2.y^2 = 2*a*(p1.x -p2.x)+ 2*b*(p1.y-p2.y)
  // means
  //b = bisectp1p2.y + pslope*(a - bisectp1p2.x)
  //means
  //p1.x^2 -p2.x^2 + p1.y^2 -p2.y^2 = 2*a*(p1.x -p2.x)+ 
  //         2*(bisectp1p2.y + pslope*(a - bisectp1p2.x))*(p1.y-p2.y)
  // means
  //p1.x^2 -p2.x^2 + p1.y^2 -p2.y^2 - 2*(bisectp1p2.y - pslope*bisectp1p2.x)*(p1.y-p2.y)=
  //        2*a*(p1.x -p2.x)+ 2*pslope*a*(p1.y-p2.y)
  //means 
  //a = (p1.x^2 -p2.x^2 + p1.y^2 -p2.y^2 - 2*(bisectp1p2.y - pslope*bisectp1p2.x)*(p1.y-p2.y))/(2*(p1.x -p2.x)+ 2*pslope*a*(p1.y-p2.y))
  float a = (pow(p1.x,2) - pow(p2.x,2) + pow(p1.y,2) - pow(p2.y,2) - 2*(bisectp1p2.y - pslope*bisectp1p2.x)*(p1.y-p2.y))/(2.0*(p1.x -p2.x)+ 2.0*pslope*(p1.y-p2.y));
  float b = bisectp1p2.y + pslope*(a - bisectp1p2.x);
  Centerout.x = a;
  Centerout.y = b;
}

//sdf is subdivision factor
//getNGonSubdivisionPoints works for linear edges, wrap argument is true if 
// to work subdivision over the entirety of the polygon
void getNGonSubdivisionPoints(int sdf, PVector[] pos, PVector[] out,
                              boolean wrap){
  ////println("Pos length: ", pos.length);
  int iterval = pos.length;
  if (!wrap){
    iterval = pos.length-1;
  }
  for (int i=0; i< iterval; i++){
    int ni = i+1;
    if (i == pos.length -1){
      ni = 0;
    }
    PVector pos1 = pos[i];
    PVector pos2 = pos[ni];
    //using vector approach in subdividing points.
    PVector p12 = PVector.sub(pos2,pos1);
    float mp12 = p12.mag();
    p12.normalize();
    mp12 /= float(sdf);
    p12 = PVector.mult(p12, mp12);
    ////println("P12: ", p12);
    PVector newpos = pos1;
    for (int j=0; j < sdf; j++){
      int ind = (sdf)*i + j;
      out[ind] = newpos;
      newpos = PVector.add(newpos,p12);
    }
    
  }
}

void getNGonSubdivisionPoints(int sdf, Edge cedge, ArrayList<PVector> out,
                              HashMap<Integer,Integer> vRemapOut){
  ////println("Pos length: ", pos.length);
  
  //for (int i=0; i< iterval; i++){
  //  int ni = i+1;
  //  if (i == pos.length -1){
  //    ni = 0;
  //  }
  PVector pos1 = cedge.p1;
  PVector pos2 = cedge.p2;
  //using vector approach in subdividing points.
  PVector p12 = PVector.sub(pos2,pos1);
  float mp12 = p12.mag();
  p12.normalize();
  mp12 /= float(sdf);
  p12 = PVector.mult(p12, mp12);
  ////println("P12: ", p12);
  PVector newpos = pos1;
  vRemapOut.put(cedge.p1index, out.size());
  for (int j=0; j < sdf; j++){
    //int ind = (sdf)*i + j;
    out.add(newpos);
    newpos = PVector.add(newpos,p12);
  }
  //}
} //finished getNGonSubdivisionPoints overloaded method

void getNGonSubdivisionPoints(int sdf, Edge cedge, ArrayList<PVector> out){
  ////println("Pos length: ", pos.length);
  
  //for (int i=0; i< iterval; i++){
  //  int ni = i+1;
  //  if (i == pos.length -1){
  //    ni = 0;
  //  }
  PVector pos1 = cedge.p1;
  PVector pos2 = cedge.p2;
  //using vector approach in subdividing points.
  PVector p12 = PVector.sub(pos2,pos1);
  float mp12 = p12.mag();
  p12.normalize();
  mp12 /= float(sdf);
  p12 = PVector.mult(p12, mp12);
  ////println("P12: ", p12);
  PVector newpos = pos1;
  //vRemapOut.put(cedge.p1index, out.size());
  for (int j=0; j < sdf; j++){
    //int ind = (sdf)*i + j;
    out.add(newpos);
    newpos = PVector.add(newpos,p12);
  }
  //}
} //finished getNGonSubdivisionPoints overloaded method

void getCentroidCircles(float frac, Polygon cpoly, ArrayList<Circle> out){
  if (cpoly.vertices.size() == 3){
    for (int i = 0; i < cpoly.vertices.size();i++){
      int ni = (i+1) % cpoly.vertices.size();
      //vector in the direction of the vertex from center
      PVector cv = PVector.sub(cpoly.vertices.get(i),cpoly.center);
      float dcv = cv.mag();
      dcv *= frac;
      cv.normalize();
      //CentroidCircle position is then given by 
      PVector centCircle = PVector.add(cpoly.center,PVector.mult(cv,dcv));
      //find the radius from centroidCircle...
      //There are 2 conditions in producing such radius.  
      // 1) is given by the nearest distance from a point to line
      // 2) the other is given by nearest point to curve (where an edge is an arc)

      // this is done by multiplying frac * distance centroidCircle to nearest line
      //get the edge data
      ArrayList<PVector> edgpts = new ArrayList<PVector>();
      edgpts.add(cpoly.vertices.get(i));
      edgpts.add(cpoly.vertices.get(ni));
      int edgind = cpoly.pointsToEdge.get(edgpts);
      Edge edge = cpoly.edges.get(edgind);
      float d;
      if (edge.linear){
        d = distPointToLine(centCircle, edge.p1,edge.p2);
      }
      else{
        
        d = distPointToCircle(edge.circle, centCircle);
      }
      float radius = d*(1-frac);
      radius *= .8;
      out.add(new Circle(centCircle, radius));
    }
  }
  else if (cpoly.vertices.size() == 4){
    for (int i = 0; i < cpoly.edges.size();i++){
      if (cpoly.edges.get(i).thick){
        ArrayList<PVector> pts = new ArrayList<PVector>();
        if (cpoly.edges.get(i).linear){
          getNGonSubdivisionPoints(2, cpoly.edges.get(i), pts);
        }
        else{
          getCurveSubdivisionPoints(2, cpoly.edges.get(i), pts);
        }
        
        PVector cv = PVector.sub(pts.get(1),cpoly.center);
        float dcv = cv.mag();
        dcv *= 1.4*frac;
        cv.normalize();
        //CentroidCircle position is then given by 
        PVector centCircle = PVector.add(cpoly.center,PVector.mult(cv,dcv));
        float d;
        if (cpoly.edges.get(i).linear){
          d = distPointToLine(centCircle, cpoly.edges.get(i).p1,cpoly.edges.get(i).p2);
        }
        else{
          
          d = distPointToCircle(cpoly.edges.get(i).circle, centCircle);
        }
        float radius = d*.6;
        out.add(new Circle(centCircle, radius));  
      }
    }
  }
  else if (cpoly.vertices.size() >= 5){
    float d;
    if (cpoly.edges.get(0).linear){
      d = distPointToLine(cpoly.center, cpoly.edges.get(0).p1,cpoly.edges.get(0).p2);
    }
    else{
      
      d = distPointToCircle(cpoly.edges.get(0).circle, cpoly.center);
    }
    float radius = d*(1-frac);
    out.add(new Circle(cpoly.center, radius)); 
  }
}

void getCurveSubdivisionPoints(int sdf, Edge cedge, ArrayList<PVector> out,
                               HashMap<Integer,Integer> vRemapOut){
  float p1p2angle = cedge.angle2 - cedge.angle1;
  float p1p2angleinc = p1p2angle/float(sdf);
  float cangle = cedge.angle1;
  vRemapOut.put(cedge.p1index, out.size());
  for (int i =0; i < sdf; i++){
    PVector pos = new PVector(0.0,0.0,0.0);
    Polarcoord(cangle, cedge.genR, pos);
    pos = PVector.add(pos, cedge.genCenter);
    out.add(pos);
    cangle += p1p2angleinc;
  }
}

void getCurveSubdivisionPoints(int sdf, Edge cedge, ArrayList<PVector> out){
  float p1p2angle = cedge.angle2 - cedge.angle1;
  float p1p2angleinc = p1p2angle/float(sdf);
  float cangle = cedge.angle1;
  //vRemapOut.put(cedge.p1index, out.size());
  for (int i =0; i < sdf; i++){
    PVector pos = new PVector(0.0,0.0,0.0);
    Polarcoord(cangle, cedge.genR, pos);
    pos = PVector.add(pos, cedge.genCenter);
    out.add(pos);
    cangle += p1p2angleinc;
  }
}

void setArcincs(Polygon cpoly, float atime){
  //this is run sequentially after polygon arcs have been computed
  float sdf = atime*60.0f;
  ArrayList<Edge> arcs = cpoly.arcs;
  for (int i = 0; i < arcs.size(); i++){
    Edge arc = arcs.get(i);
    float p1p2angle = arc.angle2 - arc.angle1;
    float p1p2angleinc = p1p2angle/sdf;    
    cpoly.arcincs.add(p1p2angleinc);
  }
}

float xfromLine(PVector p1, PVector p2, float y){
  return (p2.x-p1.x)*(y-p1.y)/(p2.y-p1.y)+p1.x;
}

void ptSetinPolygon(Polygon cpoly, ArrayList<PVector> pts, ArrayList<PVector> out){
  //point in convex polygon testing
  for (int i = 0; i < pts.size(); i++){
    //get min max y for x test y interior for pt in pts
    boolean test1 = false;
    ArrayList<ArrayList<PVector>> edges = new ArrayList<ArrayList<PVector>>();
    for (int j = 0; j < cpoly.vertices.size(); j++){
      ArrayList<PVector> edge = new ArrayList<PVector>();
      float yi = cpoly.vertices.get(j).y;
      float nyi = cpoly.vertices.get((j+1)% cpoly.vertices.size()).y;
      float py = pts.get(i).y;
      boolean e1 = yi <= py;
      boolean e2 = py <= nyi;
      boolean e3 = yi >= py;
      boolean e4 = py >= nyi;
      if ((e1 && e2)|| (e3 && e4)){
        test1 = true;
        edge.add(cpoly.vertices.get(j));
        edge.add(cpoly.vertices.get((j+1)% cpoly.vertices.size()));
        edges.add(edge);
      }
    }
    if (edges.size() == 2){
      ArrayList<PVector> edg1 = edges.get(0);
      ArrayList<PVector> edg2 = edges.get(1);
      float x1 = xfromLine(edg1.get(0), edg1.get(1), pts.get(i).y);
      float x2 = xfromLine(edg2.get(0), edg2.get(1), pts.get(i).y);
      boolean e1 = min(x1,x2) <= pts.get(i).x;
      boolean e2 = pts.get(i).x <= max(x1,x2); 
      if ( e1 && e2){
        out.add(pts.get(i));
      }
    }
  }
}

int negToPosMod(int num, int mod){
  if (num < 0){
    if (num%mod == 0){
      return num%mod;
    }
    else{
      int mult = abs(num/mod)+1;
      return num + mult*mod;
    }
  }
  else {
    return num;
  }
}

void writeCircleMapData(HashMap<Integer,ArrayList<Integer>> vertToVertPair,
                        ArrayList<PVector> subdivpts,
                        CircleMap circlemapout, int i, boolean Pole){
    
    int ipi = negToPosMod(i-1,  subdivpts.size()) % subdivpts.size();
    int iipi = negToPosMod(ipi-2,  subdivpts.size()) % subdivpts.size();
    int pole1 = negToPosMod(ipi-1,  subdivpts.size()) % subdivpts.size();
    int pole2 = negToPosMod(ipi-1,  subdivpts.size()) % subdivpts.size();
    boolean pole = true;
    if (!Pole){
      iipi = negToPosMod(ipi-1,  subdivpts.size()) % subdivpts.size();
      pole1 = (i + 1) % subdivpts.size();
      pole2 = negToPosMod(i-4,  subdivpts.size()) % subdivpts.size();
      pole = false;
    }
        //if (vertToVertPair.containsKey(ipi) && vertToVertPair.containsKey(iipi)){
    ArrayList<Integer> pair1 = vertToVertPair.get(i);
    PVector P2;
    Integer P2index;
    if (pair1.get(0) == i){
      P2 = subdivpts.get(pair1.get(1));
      P2index = pair1.get(1);
    }
    else {
      P2 = subdivpts.get(pair1.get(0));
      P2index = pair1.get(0);
    }
    PVector P1 = subdivpts.get(i);
    Integer P1index = i;
    Circle c1 = circlemapout.vertPairToCircle.get(pair1);
    println("index 1: ",P1index);
    println("index 2: ",P2index);
    Edge edg = new Edge(P2, P1, c1.center, c1.radius, P2index, P1index);
    //repeat getting interior edge/arc/circle 1
    edg.interiornp1 = ipi;
    edg.interiornp2 = iipi;
    edg.pole1 = pole1;
    edg.pole2 = pole2;
    edg.pole = pole;
    ArrayList<Edge> iedges = new ArrayList<Edge>();
    ArrayList<Integer> ipair1 = vertToVertPair.get(ipi);
    PVector iP2;
    Integer iP2index;
    if (ipair1.get(0) == ipi){
      iP2 = subdivpts.get(ipair1.get(1));
      iP2index = ipair1.get(1);
    }
    else {
      iP2 = subdivpts.get(ipair1.get(0));
      iP2index = ipair1.get(0);
    }
    PVector iP1 = subdivpts.get(ipi);
    Integer iP1index = ipi;
    Circle ic1 = circlemapout.vertPairToCircle.get(ipair1);

    Edge iedg = new Edge(iP1, iP2, ic1.center, ic1.radius, iP1index, iP2index);
    //
    ArrayList<Integer> iipair1 = vertToVertPair.get(iipi);
    //println(vertToVertPair);
    //println(iipi);
    PVector iiP2;
    Integer iiP2index;
    if (iipair1.get(0) == iipi){
      iiP2 = subdivpts.get(iipair1.get(1));
      iiP2index = iipair1.get(1);
    }
    else {
      iiP2 = subdivpts.get(iipair1.get(0));
      iiP2index = iipair1.get(0);
    }
    PVector iiP1 = subdivpts.get(iipi);
    Integer iiP1index = iipi;
    Circle iic1 = circlemapout.vertPairToCircle.get(iipair1);
    Edge iiedg = new Edge(iiP1, iiP2, iic1.center, iic1.radius, iiP1index, iiP2index);
    iedges.add(iedg);
    iedges.add(iiedg);
    circlemapout.interiorEdges.put(edg,iedges);
        //}
}

void getCircleCenter(Polygon cpoly, ArrayList<PVector> subdivpts, 
                     int i, int ni, PVector Centerout){
  if ((i+3) % subdivpts.size() == ni) {
    Edge pedge = cpoly.getEdgefromSubdivPt(i);
    if (pedge.linear){
      Centerout.x = (subdivpts.get(i).x + subdivpts.get(ni).x)/2.0;
      Centerout.y = (subdivpts.get(i).y + subdivpts.get(ni).y)/2.0;
      PVector cp1 = PVector.sub(Centerout,subdivpts.get(i));
      float radius = cp1.mag();
      cp1.rotate(PI/2.0);
      cp1.normalize();
      PVector nOrth = PVector.mult(cp1, -.4*radius);
      Centerout.x = PVector.add(nOrth, Centerout).x;
      Centerout.y = PVector.add(nOrth, Centerout).y;
      
    }
    else{
      PVector c = pedge.getEdgeHalfAnglePt();
      Centerout.x = c.x;
      Centerout.y = c.y;
    }
  }
  else{
    Centerout.x = subdivpts.get((i+2)%subdivpts.size()).x;
    Centerout.y = subdivpts.get((i+2)%subdivpts.size()).y;
  }
}

void getCircleCenter4S(Polygon cpoly, ArrayList<PVector> subdivpts, 
                       int i, int ni, PVector Centerout){

    Edge pedge = cpoly.getEdgefromSubdivPt(i);
    if (pedge.linear){
      Centerout.x = (subdivpts.get(i).x + subdivpts.get(ni).x)/2.0;
      Centerout.y = (subdivpts.get(i).y + subdivpts.get(ni).y)/2.0;
      PVector cp1 = PVector.sub(Centerout,subdivpts.get(i));
      float radius = cp1.mag();
      cp1.rotate(PI/2.0);
      cp1.normalize();
      PVector nOrth;
      if ((ni - i)==3){
        nOrth = PVector.mult(cp1, .4*radius);
      }
      else{
        nOrth = PVector.mult(cp1, -1.6*radius);
      }
      Centerout.x = PVector.add(nOrth, Centerout).x;
      Centerout.y = PVector.add(nOrth, Centerout).y;
      
    }
    else{
      PVector c = pedge.getEdgeHalfAnglePt();
      Centerout.x = c.x;
      Centerout.y = c.y;
    }
  

}

void getSubdivPolyArcData(Polygon cpoly, ArrayList<PVector> subdivpts, 
                          HashMap<Integer,Integer> vertsRemap,
                          CircleMap circlemapout){
  Collection<Integer> flaggedVerts = vertsRemap.values();
  ArrayList<ArrayList<Integer>> complVertPairs = new ArrayList<ArrayList<Integer>>();
  HashMap<Integer,ArrayList<Integer>> vertToVertPair = new HashMap<Integer,ArrayList<Integer>>();
  //HashMap<Integer,ArrayList<Integer>> vertToPair = new HashMap<Integer,ArrayList<Integer>>();
  //HashMap<Circle,ArrayList<Circle>> interiorCircles = new HashMap<Circle,ArrayList<Circle>>();
  int liter = 2;
  for (int i = 0; i < subdivpts.size(); i++){
    ArrayList<Integer> complpair = new ArrayList<Integer>();
    if (flaggedVerts.contains(i)){
      continue;  //skip original polygon vertex (prior to subidivision)
    }
    //else 
    int ni, pi;
    if (cpoly.vertices.size() >= 5){
      ni = (i+4) % subdivpts.size();
      pi = negToPosMod(i-4,  subdivpts.size()) % subdivpts.size();
      if (flaggedVerts.contains(ni)){
        int ipi = negToPosMod(i-1,  subdivpts.size())% subdivpts.size();
        int iipi = negToPosMod(ipi-2,  subdivpts.size())% subdivpts.size();
        if (vertToVertPair.containsKey(ipi) && vertToVertPair.containsKey(iipi)){
          writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, true);
        }
        continue;  
      }
    }
    else if (cpoly.vertices.size() == 4){
      boolean t1 = i == 7 || i == 10 || i == 12 || i == 15;
      boolean t2 = i == 7 || i == 10 || i == 15;
      if (t1){
        if (t2){
          writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, true);
        }
        else{
          writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, false);
        }
        continue;
      }
      else{
        boolean t3 = i == 1 || i == 9;
        boolean t4 = i == 2 || i == 4;
        if (t3){
          ni = i + 3;
        }
        else{
          ni = (i+4)%subdivpts.size();
        }
        if (t4){
          continue;
        }
      }
    }
    else{
      int ini, ipi;
      ini = (i + 1)% subdivpts.size();
      ipi = negToPosMod(i-1,  subdivpts.size())% subdivpts.size();
      if (flaggedVerts.contains(ini) || flaggedVerts.contains(ipi)){
        ni = (i+3) % subdivpts.size();
        pi = negToPosMod(i-3,  subdivpts.size()) % subdivpts.size();
      }
      else{
        //either 2 points on the most interior of subidivided edge
        ni = (i+4) % subdivpts.size();
        pi = negToPosMod(i-4,  subdivpts.size())% subdivpts.size();
        if (flaggedVerts.contains(negToPosMod(ni-1,  subdivpts.size())% subdivpts.size())){
          ipi = negToPosMod(i-1,  subdivpts.size())% subdivpts.size();
          int iipi = negToPosMod(ipi-2,  subdivpts.size())% subdivpts.size();
          if (vertToVertPair.containsKey(ipi) && vertToVertPair.containsKey(iipi)){
            writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, true);
          }
          continue;  
        }
        //past this point means this is 2nd from pole point 
      }

      //vert pair test to see if Arc already computed on previous pair set for a
      // given vertex iteration  if true for set inclusion, continue again.
      complpair.add(pi);
      complpair.add(i);
      if (complVertPairs.contains(complpair)){
        ipi = negToPosMod(i-1,  subdivpts.size())% subdivpts.size();
        int iipi = negToPosMod(ipi-1,  subdivpts.size())% subdivpts.size();
        if (vertToVertPair.containsKey(ipi) && vertToVertPair.containsKey(iipi)){
          writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, false);
        }
        continue;  //we've already computed an arc for the present iterated vertex index i.
      }
    }

    complpair = new ArrayList<Integer>();
    PVector Centerout = new PVector(0.0,0.0,0.0);
    //getCircleCenter(subdivpts.get(i), subdivpts.get(ni), Centerout);
    if (cpoly.vertices.size()==3){
      getCircleCenter(cpoly, subdivpts, i, ni, Centerout);
    }
    else if (cpoly.vertices.size() == 4){
      getCircleCenter4S(cpoly, subdivpts, i, ni,Centerout);
    }
    else{
      getCircleCenter(cpoly, subdivpts, i, ni, Centerout);
    }
    PVector p1c = PVector.sub(subdivpts.get(i), Centerout);
    float radius = p1c.mag();
    complpair.add(i);
    complpair.add(ni);
    circlemapout.vertPairToCircle.put(complpair,new Circle(Centerout, radius));
    
    complVertPairs.add(complpair);
    vertToVertPair.put(i, complpair);
    vertToVertPair.put(ni, complpair);
  }
  //write the last circlemapout interiors.
  writeCircleMapData(vertToVertPair, subdivpts, circlemapout, liter, true);
  if (cpoly.vertices.size() >= 5){
    writeCircleMapData(vertToVertPair, subdivpts, circlemapout, 5, true);
  }
  else{
    writeCircleMapData(vertToVertPair, subdivpts, circlemapout, 4, false);
  }
}

PVector closestPoint(ArrayList<PVector> pts, PVector pos){
  ArrayList<Float> ds = new ArrayList<Float>();
  float[] dsf = new float[pts.size()];
  int i = 0;
  for (PVector pt : pts){
    PVector postopt = PVector.sub(pt,pos);
    ds.add(postopt.mag());
    dsf[i] = postopt.mag();
    i += 1;
  }
  return pts.get(ds.indexOf(min(dsf)));
}

void getCentCirclesfromPoles(Polygon cpoly, int pole1, int pole2, ArrayList<Circle> centcircles,
                             ArrayList<Circle> out){
  if (cpoly.vertices.size() == 3){
    out.add(centcircles.get(pole1));
    if (pole1 != pole2){
      out.add(centcircles.get(pole2));
    }
  }
  else if (cpoly.vertices.size()== 4){
    if (pole1 == 0 || pole1 == 3){
      out.add(centcircles.get(1));
    }
    else if (pole1 == 1 || pole1 == 2){
      out.add(centcircles.get(0));
    }
    if (pole1 != pole2){
      if (pole2 == 0 || pole2 == 3){
        out.add(centcircles.get(1));
      }
      else if (pole2 == 1 || pole2 == 2){
        out.add(centcircles.get(0));
      }
    }
  }
  else{
    out.add(centcircles.get(0));
  }
}

void getWindingOrder(int v1, int v2, int modop, ArrayList<Integer> out){
  boolean t1 = v2 == (modop-1);
  boolean t2 = v1 == 0;
  boolean t3 = t1 && t2;
  boolean t4 = v1 == (modop-1);
  boolean t5 = v2 == 0;
  boolean t6 = t4 && t5;
  if (v1 % modop < v2 % modop  && !t3){
    out.add(v1);
    out.add(v2);
  }
  else if (v1 % modop > v2 % modop && !t6){
    out.add(v2);
    out.add(v1);
  }
  else if (v1 % modop < v2 % modop  && t3){
    out.add(v2);
    out.add(v1);
  }
  else if (v1 % modop > v2 % modop && t6){
    out.add(v1);
    out.add(v2);
  }
}

void buildEdgefromParent(PVector pt1, PVector pt2, Integer np1index,
                         Integer np2index, Edge parent, Polygon opoly,
                         Integer op1index, Integer op2index,
                          ArrayList<Edge> outs){
  //function call is given from edge data constructed from point subdivision data 
  //on the polygon only...doesn't work for interior subdivision points.
  //Use alternate overloaded method for the other edge construction type for 
  //interior subdivision points.
  //opoly is the original polygon (subdivision point set).
  //np1index is new polygon indexing 1   (On the new polygon point set).
  //op1index is the original polygon indexing 1. (On the set of subivision points set)
  ArrayList<Integer> ptpair = new ArrayList<Integer>();
  ptpair.add(op1index);
  ptpair.add(op2index);
  Boolean ethick = opoly.SubdivPtPairtoThick.get(ptpair);
  if (parent.linear){
    //this constructor call is made..
    //Edge(PVector P1, PVector P2, boolean Thick, int P1index, int P2index)
    Edge out = new Edge(pt1, pt2, ethick, np1index, np2index);
    outs.add(out);
  }
  else{
    //Edge(PVector P1, PVector P2, boolean Thick, PVector GenCenter,
    //   float GenR, float Angle1, float Angle2, int P1index, int P2index)
    PVector CToP1 = PVector.sub(pt1,parent.genCenter);
    PVector CToP2 = PVector.sub(pt2,parent.genCenter);
    float a1 = CToP1.heading();
    float a2 = CToP2.heading();
    Edge out = new Edge(pt1, pt2, ethick, parent.genCenter, 
                   parent.genR, a1, a2, np1index, np2index);
    outs.add(out);
  }
}

void buildEdgefromParent(PVector pt1, PVector pt2, Integer np1index,
                         Integer np2index, Edge parent, Boolean edgeThick,
                          ArrayList<Edge> outs){
  //alternate overloaded method for the other edge construction type for 
  //interior subdivision points.
  //opoly is the original polygon (subdivision point set).
  //np1index is new polygon indexing 1   (On the new polygon point set).
  //op1index is the original polygon indexing 1. (On the set of subivision points set)

  if (parent.linear){
    //this constructor call is made..
    //Edge(PVector P1, PVector P2, boolean Thick, int P1index, int P2index)
    Edge out = new Edge(pt1, pt2, edgeThick, np1index, np2index);
    outs.add(out);
  }
  else{
    //Edge(PVector P1, PVector P2, boolean Thick, PVector GenCenter,
    //   float GenR, float Angle1, float Angle2, int P1index, int P2index)
    PVector CToP1 = PVector.sub(pt1,parent.genCenter);
    //println("parent gencenter: ", parent.genCenter);
    //println("pt2", pt2);
    PVector CToP2 = PVector.sub(pt2,parent.genCenter);
    float a1 = CToP1.heading();
    float a2 = CToP2.heading();
    Edge out = new Edge(pt1, pt2, edgeThick, parent.genCenter, 
                   parent.genR, a1, a2, np1index, np2index);
    outs.add(out);
  }
}

void buildInteriorPolygon(PVector[] verts, Boolean[] edgThcks, Edge[] ParentEdges,
                          ArrayList<Polygon> pouts){
  HashMap<ArrayList<PVector>,Integer> pointsToEdge = new HashMap<ArrayList<PVector>,Integer>();
  ArrayList<Edge> Edges = new ArrayList<Edge>();
  for (int i = 0; i < verts.length; i++){
    int ni = (i+1)%verts.length;
    buildEdgefromParent(verts[i], verts[ni], i, ni, ParentEdges[i], edgThcks[i], Edges);
    ArrayList<PVector> ptpairal = new ArrayList<PVector>();
    PVector[] ptpair = {verts[i],verts[ni]};
    Collections.addAll(ptpairal,ptpair);
    pointsToEdge.put(ptpairal,Edges.size()-1);
    //Edges.add(out);
  }
  ArrayList<PVector> vertsal = new ArrayList<PVector>();
  Collections.addAll(vertsal,verts);
  PVector PCenter = new PVector(0.0,0.0,0.0);
  PolygonCentroid(vertsal, PCenter);
  Polygon pout = new Polygon(vertsal, PCenter);
  pout.edges = Edges;
  pout.pointsToEdge = pointsToEdge;
  pouts.add(pout);
}

PVector getOppositePt(PVector pt, ArrayList<PVector> pts){
  //opposite point in a 2 point set
  int ni = (pts.indexOf(pt)+1 )% pts.size();
  return pts.get(ni);
}

void initializePolyHoldings(Polygon cpoly, ArrayList<PolyHolding> centPolys){
  int iSize = 0;
  if (cpoly.vertices.size() == 3){
    iSize = 7;
  }
  else if (cpoly.vertices.size() == 4){
    iSize = 3;
  }
  else if (cpoly.vertices.size() >= 5){
    iSize = 1;
  }
  for (int i = 0; i < iSize; i++){
    PolyHolding polyholding = new PolyHolding();
    centPolys.add(polyholding);
  }
}

//void getPolyHoldindices(Polygon cpoly, Integer subdivPole){
//  if (cpoly.vertices.size() == 3){
//     int pole = cpoly.subdivPtToPt.get(subdivPole);
//     int prevpole = (pole-1)%3;
//     int nextpole = (pole+1)%3;
//  }
//}

void writeDistantPoints(Polygon cpoly, int pole, ArrayList<PVector> ISPts,
                        ArrayList<PVector> ISPts2, PVector cPt1, PVector cPt2,
                        Edge parentEdge1, Edge parentEdge2,
                        Edge parentEdge3, ArrayList<PolyHolding> centPolys){
  //checks polygon type and indicates whether distant interior subdivision
  //circle intersect points need to be appended to polyholding
  //ISPts is the set of circle circle intersection point where the circle intersection
  //is from the centroid circle relative the interior pole arc.
  //pole is indexed on the original polygon vertices not subdivision vertices.
  //See Dohecahedral Subdivision Rule diagram for modeling help here.
  if (cpoly.vertices.size()==3){
    PVector dPt1 = getOppositePt(cPt1, ISPts);
    PVector dPt2 = getOppositePt(cPt2, ISPts2);
    PointEdgetoPoly pep1 = new PointEdgetoPoly(3, 4, dPt1, parentEdge1);
    PointEdgetoPoly pep2 = new PointEdgetoPoly(1, 2, dPt1, parentEdge3);
  
    PointEdgetoPoly pep3 = new PointEdgetoPoly(1, 2, dPt2, parentEdge3);
    PointEdgetoPoly pep4 = new PointEdgetoPoly(0, 1, dPt2, parentEdge2);
    int ni = (pole+1) % 3;
    int pi = negToPosMod(pole-1, 3) % 3;
    int ni2 = 3+pole;
    int pi2 = 3+pi;
    (centPolys.get(ni)).addPointEdgetoPoly(pep1);
    (centPolys.get(ni2)).addPointEdgetoPoly(pep2);
    (centPolys.get(pi)).addPointEdgetoPoly(pep3);
    (centPolys.get(pi2)).addPointEdgetoPoly(pep4);
  }
  else if (cpoly.vertices.size()==4){
    PVector dPt1 = getOppositePt(cPt1, ISPts);
    PVector dPt2 = getOppositePt(cPt2, ISPts2);
    PointEdgetoPoly pep1 = new PointEdgetoPoly();
    PointEdgetoPoly pep2 = new PointEdgetoPoly();
    //PointEdgetoPoly pep3;
    //PointEdgetoPoly pep4; 
    //pep1 = new PointEdgetoPoly(3, 4, dPt1, parentEdge1);
    //pep2 = new PointEdgetoPoly(0, 1, dPt1, parentEdge3);
    if (pole%2 == 0){
      pep1 = new PointEdgetoPoly(0, 1, dPt1, parentEdge1);
      if (pole < 2){
        pep2 = new PointEdgetoPoly(1, 2, dPt1, parentEdge3);
      }
      else{
        pep2 = new PointEdgetoPoly(3, 0, dPt1, parentEdge3);
      }
    }
    else{
      pep1 = new PointEdgetoPoly(4, 0, dPt2, parentEdge3);
      if (pole < 2){
        pep2 = new PointEdgetoPoly(0, 1, dPt2, parentEdge2);
      }
      else{
        pep2 = new PointEdgetoPoly(2, 3, dPt2, parentEdge2);
      }
      //pep2 = new PointEdgetoPoly(0, 1, dPt2, parentEdge2);
    }
    int ni;
    if (pole == 0){
      ni = 1;
    }
    else if (pole==1){
      ni = 0;
    }
    else if (pole==2){
      ni = 0;
    }
    else{
      ni = 1;
    }
    int ni2 = 2;
    (centPolys.get(ni)).addPointEdgetoPoly(pep1);
    (centPolys.get(ni2)).addPointEdgetoPoly(pep2);
  }
}

void writeDistantPointNonPolePass(Polygon cpoly,int pole, PVector cPt1,
                        Edge parentEdge1, Edge parentEdge2 , 
                        ArrayList<PolyHolding> centPolys){
  //an alterate distant point write method to PolyHolding Class on a non pole pass.
  //Special case write method only...only for 3 gon original polygon type.
  //inputting 'forward' pole on a two pole differenced arc.
  if (cpoly.vertices.size() == 3){
    int npole = (pole+1) % 3;
    int nnpole = (npole+1) % 3;
    int qpole = 3+pole;
    int q2pole = 3+npole;
    PointEdgetoPoly pep1 = new PointEdgetoPoly(npole, nnpole, cPt1, parentEdge1);
    PointEdgetoPoly pep2 = new PointEdgetoPoly(2, 3, cPt1, parentEdge1);
    PointEdgetoPoly pep3 = new PointEdgetoPoly(3, 0, cPt1, parentEdge2);
    PointEdgetoPoly pep4 = new PointEdgetoPoly(2, 3, cPt1, parentEdge2);
    (centPolys.get(6)).addPointEdgetoPoly(pep1);
    (centPolys.get(npole)).addPointEdgetoPoly(pep2);
    (centPolys.get(qpole)).addPointEdgetoPoly(pep4);
    (centPolys.get(q2pole)).addPointEdgetoPoly(pep3);
  }
}

void writeClosePoint(Polygon cpoly, int pole,  PVector Pt1, PVector Pt2,
                        Edge parentEdge1, Edge parentEdge2,  
                        ArrayList<PolyHolding> centPolys){
  if (cpoly.vertices.size() == 3){
    PointEdgetoPoly pep1 = new PointEdgetoPoly(0,1, Pt1, parentEdge1);
    PointEdgetoPoly pep2 = new PointEdgetoPoly(4,0, Pt2, parentEdge2);
    (centPolys.get(pole)).addPointEdgetoPoly(pep1);
    (centPolys.get(pole)).addPointEdgetoPoly(pep2);
  }
  else if (cpoly.vertices.size() == 4){
    int ni;

    if (pole == 0){
      ni = 0;
    }
    else if (pole==1){
      ni = 1;
    }
    else if (pole==2){
      ni = 1;
    }
    else{
      ni = 0;
    }
    if (pole%2 == 0){
      PointEdgetoPoly pep1 = new PointEdgetoPoly(3,4, Pt1, parentEdge1);
      PointEdgetoPoly pep2 = new PointEdgetoPoly(2,3, Pt2, parentEdge2);
      (centPolys.get(ni)).addPointEdgetoPoly(pep1);
      (centPolys.get(ni)).addPointEdgetoPoly(pep2);
    }
    else{
      PointEdgetoPoly pep1 = new PointEdgetoPoly(2,3, Pt1, parentEdge1);
      PointEdgetoPoly pep2 = new PointEdgetoPoly(1,2, Pt2, parentEdge2);
      (centPolys.get(ni)).addPointEdgetoPoly(pep1);
      (centPolys.get(ni)).addPointEdgetoPoly(pep2);
    }
  }
  else{
    PointEdgetoPoly pep1 = new PointEdgetoPoly(pole,(pole+1)%cpoly.vertices.size(), Pt1, parentEdge1);
    (centPolys.get(0)).addPointEdgetoPoly(pep1);
  }

}

void writeClosePoint(Polygon cpoly, int pole,  PVector Pt1, int poleiter,
                        Edge parentEdge1,  
                        ArrayList<PolyHolding> centPolys){
    //poleiter is the iteration on the present pole.  There are 2 forware arc iterations per
    //pole when winding on the edge to edges mappings in the buildInteriorPolyons method
    //for a 3 gon subdivision build, and 3 forward arc iterations per pole on the 4 gon,
    //and 5 forward arc iterations per pole on the 5 gon (given there is only one pole centroid).
    // A forward arc is the 'next' interior edge arc as opposed to 'previous' 
    //arc edge where 'next' is characterized one closest to the next pole (pole+1) mod poleset
    //and 'prev' is characterized as the interior arc closest to the previous pole (pole-1) mod poleset.
    //'forward' arcs are only considered in this write method.
    //The poleiteration is the number of times in winding around the centPoly.
    PointEdgetoPoly pep1 = new PointEdgetoPoly();
    if (cpoly.vertices.size()==3){
      
      if (poleiter == 0){
        pep1 = new PointEdgetoPoly(4, 0, Pt1, parentEdge1);
      }
      else if (poleiter == 1){
        pep1 = new PointEdgetoPoly(0, 1, Pt1, parentEdge1);
      }
      if (pole == 0 && poleiter == 0){
        pep1 = new PointEdgetoPoly(0, 1, Pt1, parentEdge1);
      }
      else if (pole == 0 && poleiter == 1){
        pep1 = new PointEdgetoPoly(4, 0, Pt1, parentEdge1);
      }
      (centPolys.get(pole)).addPointEdgetoPoly(pep1);
    }
    if (cpoly.vertices.size()==4){
      int ni;
      int[] polev1 = {3,1,2};
      int[] polev2 = {1,2,3};
      if (pole == 0){
        ni = 1;
      }
      else if (pole==1){
        ni = 0;
      }
      else if (pole==2){
        ni = 0;
      }
      else{
        ni = 1;
      }
      if (ni == 1){
        pep1 = new PointEdgetoPoly(polev2[poleiter], polev2[poleiter]+1, 
                                   Pt1, parentEdge1);
      }
      else{
        pep1 = new PointEdgetoPoly(polev1[poleiter], polev1[poleiter]+1, 
                                   Pt1, parentEdge1);
      }
      (centPolys.get(ni)).addPointEdgetoPoly(pep1);
    }
    else{
      pep1 = new PointEdgetoPoly(poleiter, poleiter+1, Pt1, parentEdge1);
      (centPolys.get(0)).addPointEdgetoPoly(pep1);
    }
}

void initializePoleiterationMap(Polygon cpoly, HashMap<Integer,Integer> poleiteration){
  if (cpoly.vertices.size() == 3){
    poleiteration.put(0,0);
    poleiteration.put(1,0);
    poleiteration.put(2,0);
  }
  else if (cpoly.vertices.size() == 4){
    poleiteration.put(0,0);
    poleiteration.put(1,0);
  }
  else{
    poleiteration.put(0,0);
  }
}

void iteratePoleiterator(Integer pole, HashMap<Integer,Integer> poleiteration){
  int ival = poleiteration.get(pole)+1;
  poleiteration.put(pole,ival);
}

Boolean[] getThickPolybool(PVector[] verts){
  if (verts.length == 3){
    return new Boolean[] {false,false,false};
  }
  else if (verts.length == 4){
    return new Boolean[] {false,true,false,true};
  }
  else{
    return new Boolean[] {true,true,true,true,true};
  }
}

void buildInteriorPolygons(Polygon cpoly, ArrayList<Circle> centroidcircles, 
                           CircleMap circlemap, ArrayList<Polygon> outPolys){
  
  HashMap<Edge,ArrayList<Edge>> interioredges = circlemap.interiorEdges;
  ArrayList<PolyHolding> centPolys = new ArrayList<PolyHolding>();
  //arraylist polyholding container is addressed in winding order with 
  //interior most polygons around the centroid. For example, on the 3 gon 
  // the poles are used in indexing the 5 gons with reserved addresses 0,1,2
  // then the 4 gons are addressed 3,4,5, and finally the interior most polygon
  //at address 6.
  initializePolyHoldings(cpoly, centPolys);
  HashMap<Integer,Integer> poleiteration = new HashMap<Integer,Integer>();
  initializePoleiterationMap(cpoly, poleiteration);
  //println("starting polygons creation");
  int k = 0;
  for (Map.Entry<Edge,ArrayList<Edge>> me : interioredges.entrySet()) {
    Edge pedge = me.getKey();
    //save edge to parent polygon arc data
    pedge.computeAngles();
    pedge.computeRadius();
    cpoly.arcs.add(pedge);
    Edge iedge = me.getValue().get(0);
    Edge iedge2 = me.getValue().get(1);
    
    if (pedge.pole){
      ArrayList<PVector> ipts = new ArrayList<PVector>();
      CircleCircleIntersection(pedge.circle, iedge.circle, ipts);
      ArrayList<PVector> iptsinpoly = new ArrayList<PVector>();
      //println("pedge center: ", pedge.circle.center);
      //println("iedge center: ", iedge.circle.center);
      //println("pedge radius: ", pedge.circle.radius);
      //println("iedge radius: ", iedge.circle.radius);
      //println("Ipts: ", ipts);
      ptSetinPolygon(cpoly, ipts, iptsinpoly);
      //println("Iptsinpoly: " ,iptsinpoly);
      //PVector ipt1 = closestPoint(iptsinpoly, cpoly.subdivpts.get(pedge.pole1));
      PVector ipt1 = closestPoint(iptsinpoly, pedge.p2);
      ipts = new ArrayList<PVector>();
      CircleCircleIntersection(pedge.circle, iedge2.circle, ipts);
      iptsinpoly = new ArrayList<PVector>();
      ptSetinPolygon(cpoly, ipts, iptsinpoly);
      //PVector ipt2 = closestPoint(iptsinpoly, cpoly.subdivpts.get(pedge.pole1));
      PVector ipt2 = closestPoint(iptsinpoly, pedge.p1);
      ipts = new ArrayList<PVector>();
      ArrayList<Circle> pecCircles = new ArrayList<Circle>();
      int vp1 = cpoly.subdivPtToPt.get(pedge.pole1);
      int vp2 = cpoly.subdivPtToPt.get(pedge.pole2);
      getCentCirclesfromPoles(cpoly, vp1, vp2, centroidcircles,pecCircles);
      Circle peccircle = pecCircles.get(0);
      CircleCircleIntersection(pedge.circle, peccircle, ipts);
      PVector ipt3 = closestPoint(ipts, cpoly.subdivpts.get(pedge.interiornp1));
      PVector ipt4 = getOppositePt(ipt3, ipts);
      ipts = new ArrayList<PVector>();
      CircleCircleIntersection(iedge.circle, peccircle, ipts);
      PVector ipt5 = closestPoint(ipts, cpoly.subdivpts.get(pedge.pole1));
      ArrayList<PVector> ipts2 = new ArrayList<PVector>();
      CircleCircleIntersection(iedge2.circle, peccircle, ipts2);
      PVector ipt6 = closestPoint(ipts2, cpoly.subdivpts.get(pedge.pole1));
      //first polygon is the pole 5 sided polygon
      //get the original edge data...this is inheritance data for parent subdivided
      //edges.
      int opi = cpoly.subdivPtToPt.get(pedge.pole1);
      int nopi = (opi+1)%cpoly.vertices.size();
      int popi = negToPosMod(opi-1, cpoly.vertices.size())%cpoly.vertices.size();
      ArrayList<Integer> polpair = new ArrayList<Integer>();
      getWindingOrder(opi, nopi, cpoly.vertices.size(), polpair);
      PVector[] polpairvec= {cpoly.vertices.get(polpair.get(0)), 
                             cpoly.vertices.get(polpair.get(1))};
      ArrayList<PVector> polpairvecal = new ArrayList<PVector>();
      Collections.addAll(polpairvecal,polpairvec);
      ArrayList<Integer> polpair2 = new ArrayList<Integer>();
      getWindingOrder(opi, popi, cpoly.vertices.size(), polpair2);
      PVector[] polpairvec2= {cpoly.vertices.get(polpair2.get(0)), 
                              cpoly.vertices.get(polpair2.get(1))};
      ArrayList<PVector> polpairvecal2 = new ArrayList<PVector>();
      Collections.addAll(polpairvecal2,polpairvec2);      
      //1rst polygon//
      PVector ept1 = cpoly.subdivpts.get(pedge.pole1);
      PVector ept2 = cpoly.subdivpts.get(pedge.interiornp1);
      PVector ept3 = cpoly.subdivpts.get(pedge.interiornp2);
      PVector ept4 = pedge.p1;
      PVector ept5 = pedge.p2;
      PVector[] verts = {ept1,ept2,ipt5,ipt6,ept3};
      ArrayList<PVector> vertsal = new ArrayList<PVector>();
      Collections.addAll(vertsal,verts);
      Boolean[] edgThcks = {true,true,true,true,true}; 
      
      Edge[] ParentEdges = {cpoly.edges.get(cpoly.pointsToEdge.get(polpairvecal)),iedge,new Edge(peccircle),iedge2,
                            cpoly.edges.get(cpoly.pointsToEdge.get(polpairvecal2))};
      
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //polygon 2
      verts = new PVector[] {ipt5,ipt1,ipt2,ipt6};
      edgThcks = new Boolean[] {false,true,false,true};
      ParentEdges = new Edge[] {iedge,pedge,iedge2, new Edge(peccircle)};
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //check to see that centroid circle intersect distant points need to be
      //written for current polygon type. 
      int vpole = cpoly.subdivPtToPt.get(pedge.pole1);
      writeDistantPoints(cpoly, vpole, ipts, ipts2, ipt5, ipt6, iedge, iedge2,
                        new Edge(peccircle), centPolys);
      //write nearest points 
      //writeClosePoint(cpoly, vp1,  ipt1, poleiteration.get(vp1), iedge,  
      //                centPolys);
      writeClosePoint(cpoly, vp1,  ipt1, ipt2, iedge, pedge, centPolys);
      //iteratePoleiterator(vp1, poleiteration);
      //polygon 3
      verts = new PVector[] {ept2,ept4,ipt3,ipt5};
      edgThcks = new Boolean[] {false,true,false,true};
      ParentEdges = new Edge[] {cpoly.edges.get(cpoly.pointsToEdge.get(polpairvecal)),
                                pedge, new Edge(peccircle), iedge};
      //buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //polygon 4
      verts = new PVector[] {ipt5,ipt3,ipt1};
      edgThcks = new Boolean[] {false,false,false};
      ParentEdges = new Edge[] {new Edge(peccircle),
                                pedge,iedge};
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
    }
    else{
      ArrayList<PVector> ipts = new ArrayList<PVector>();
      CircleCircleIntersection(pedge.circle, iedge.circle, ipts);
      //println("pedge center: ", pedge.circle.center);
      //println("iedge center: ", iedge.circle.center);
      ArrayList<PVector> iptsinpoly = new ArrayList<PVector>();
      ptSetinPolygon(cpoly, ipts, iptsinpoly);
      //println("Ipts: ", ipts);
      //PVector ipt1 = closestPoint(iptsinpoly, cpoly.subdivpts.get(pedge.pole1));
      PVector ipt1 = closestPoint(iptsinpoly, iedge.p2);
      ipts = new ArrayList<PVector>();
      CircleCircleIntersection(pedge.circle, iedge2.circle, ipts);
      iptsinpoly = new ArrayList<PVector>();
      ptSetinPolygon(cpoly, ipts, iptsinpoly);
      //PVector ipt2 = closestPoint(iptsinpoly, cpoly.subdivpts.get(pedge.pole2));
      PVector ipt2 = closestPoint(iptsinpoly, iedge2.p2);
      ipts = new ArrayList<PVector>();
      ArrayList<Circle> pecCircles = new ArrayList<Circle>();
      int vp1 = cpoly.subdivPtToPt.get(pedge.pole1);
      int vp2 = cpoly.subdivPtToPt.get(pedge.pole2);
      getCentCirclesfromPoles(cpoly, vp1, vp2, centroidcircles,pecCircles);
      Circle peccircle = pecCircles.get(0);
      CircleCircleIntersection(pedge.circle, peccircle, ipts);
      PVector ipt3 = closestPoint(ipts, cpoly.subdivpts.get(pedge.interiornp2));
      PVector ipt4 = getOppositePt(ipt3, ipts);
      ipts = new ArrayList<PVector>();
      PVector ipt5;
      PVector ipt6;
      //note we shouldn't have to case structure for 5gons since they shouldn't 
      //have non pole type leading edge maps.
      Circle peccircle2 = pecCircles.get(1);
      CircleCircleIntersection(pedge.circle, peccircle2, ipts);
      ipt5 = closestPoint(ipts, cpoly.subdivpts.get(pedge.interiornp1));
      ipt6 = getOppositePt(ipt5, ipts);
      
      ipts = new ArrayList<PVector>();
      CircleCircleIntersection(iedge.circle, peccircle, ipts);
      PVector ipt7 = closestPoint(ipts, cpoly.subdivpts.get(pedge.interiornp1));
      ArrayList<PVector> ipts2 = new ArrayList<PVector>();
      CircleCircleIntersection(iedge2.circle, peccircle2, ipts2);
      PVector ipt8 = closestPoint(ipts2, cpoly.subdivpts.get(pedge.interiornp2));
      ArrayList<PVector> ipts3 = new ArrayList<PVector>();
      CircleCircleIntersection(peccircle, peccircle2, ipts3);
      PVector ipt9 = closestPoint(ipts3, cpoly.subdivpts.get(pedge.interiornp1));
      PVector ipt10 = getOppositePt(ipt9, ipts3);
      //first polygon is the pole 5 sided polygon
      //get the original edge data...this is inheritance data for parent subdivided
      //edges.
      //println("Pedge pole: ", pedge.pole);
      int opi = cpoly.subdivPtToPt.get(pedge.pole1);
      //int nopi = (opi+1)%cpoly.vertices.size();
      int popi = cpoly.subdivPtToPt.get(pedge.pole2);
      //ArrayList<Integer> polpair = new ArrayList<Integer>();
      //getWindingOrder(opi, nopi, cpoly.vertices.size(), polpair);
      //PVector[] polpairvec= {cpoly.vertices.get(polpair.get(0)), 
      //                       cpoly.vertices.get(polpair.get(1))};
      //ArrayList<PVector> polpairvecal = new ArrayList<PVector>();
      //Collections.addAll(polpairvecal,polpairvec);
      ArrayList<Integer> polpair2 = new ArrayList<Integer>();
      getWindingOrder(opi, popi, cpoly.vertices.size(), polpair2);
      PVector[] polpairvec2= {cpoly.vertices.get(polpair2.get(0)), 
                              cpoly.vertices.get(polpair2.get(1))};
      ArrayList<PVector> polpairvecal2 = new ArrayList<PVector>();
      Collections.addAll(polpairvecal2,polpairvec2);      
      //1rst polygon//
      //PVector ept1 = cpoly.subdivpts.get(pedge.pole1);
      PVector ept2 = cpoly.subdivpts.get(pedge.interiornp1);
      PVector ept3 = cpoly.subdivpts.get(pedge.interiornp2);
      PVector ept4 = pedge.p1;
      PVector ept5 = pedge.p2;
      PVector[] verts = {ept2,ipt7,ipt9,ipt8, ept3};
      ArrayList<PVector> vertsal = new ArrayList<PVector>();
      Collections.addAll(vertsal,verts);
      Boolean[] edgThcks = {true,true,true,true,true}; 
      
      Edge[] ParentEdges = {iedge,new Edge(peccircle),new Edge(peccircle2),iedge2,
                            cpoly.edges.get(cpoly.pointsToEdge.get(polpairvecal2))};
      
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //polygon 2
      verts = new PVector[] {ipt7,ipt1,ipt5,ipt9};
      edgThcks = new Boolean[] {false,true,false,true};
      ParentEdges = new Edge[] {iedge,pedge,new Edge(peccircle2), new Edge(peccircle)};
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //check to see that centroid circle intersect distant points need to be
      //written for current polygon type. 
      int vp3 = (vp2+1)%centroidcircles.size();
      Circle peccircle3 = centroidcircles.get(vp3);
      writeDistantPointNonPolePass(cpoly,vp1, ipt10, new Edge(peccircle),
                                   new Edge(peccircle2), centPolys);
      //int vpole = cpoly.subdivPtToPt.get(pedge.pole1);
      //writeDistantPoints(cpoly, vpole, ipts, ipts2, ipt5, ipt6, iedge, iedge2,
      //                  new Edge(peccircle), centPolys);
      //write nearest points 
      //writeClosePoint(cpoly, vp1,  ipt1, poleiteration.get(vp1), iedge,  
      //                centPolys);
      //iteratePoleiterator(vp1, poleiteration);
      //polygon 3 center triangle 
      verts = new PVector[] {ipt9,ipt5,ipt3};
      edgThcks = new Boolean[] {false,false,false};
      ParentEdges = new Edge[] { new Edge(peccircle2),
                                pedge, new Edge(peccircle)};
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //polygon 4
      verts = new PVector[] {ipt9,ipt3,ipt2,ipt8};
      edgThcks = new Boolean[] {false,true,false,true};
      ParentEdges = new Edge[] {new Edge(peccircle),
                                pedge, iedge2, new Edge(peccircle2)};
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //polygon 5
      verts = new PVector[] {ept2,ept5,ipt4,ipt7};
      edgThcks = new Boolean[] {false,true,false,true};
      ParentEdges = new Edge[] {cpoly.edges.get(cpoly.pointsToEdge.get(polpairvecal2)),
                                pedge, new Edge(peccircle), iedge};
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
      //polygon 6
      verts = new PVector[] {ipt7,ipt4,ipt1};
      edgThcks = new Boolean[] {false,false,false};
      ParentEdges = new Edge[] {new Edge(peccircle),
                                pedge,iedge};
      buildInteriorPolygon(verts, edgThcks, ParentEdges, outPolys);
    }
  }
  //add polygons from centPolys
  println("Interior edges: ", interioredges.size());
  for(PolyHolding centPoly : centPolys){
    println("cent poly vertices: ", centPoly.vertices);
    centPoly.writeArrayData();
    Boolean[] thickArray = getThickPolybool(centPoly.verticesArray);
    //println("Centpolys vertices array: " , centPoly.verticesArray);
    //println("Centpolys vertices length: ", centPoly.verticesArray.length);
    buildInteriorPolygon(centPoly.verticesArray, thickArray, centPoly.parentEdgesArray,
                        outPolys);
  }
  for (Circle centCircle: centroidcircles){
    cpoly.arcs.add(new Edge(centCircle));
  }
}



void buildThickEdgeDat(Edge cedge, HashMap<Integer,Integer> vertsRemap, 
                       int SubdivPtsSize, Boolean lastEdge,
                       HashMap<ArrayList<Integer>,Boolean> SubdivPtPairtoThick){
  //iterated call sequentially done from the last edge subdivision call
  //get starting point.  Note:  This is not to be sequentially called after
  //the construciton of subdivision points on a polygon, but sequentially after the 
  //subdivision points on the edge of a polygon have been formed instead.  
  Integer p1ind = cedge.p1index;
  Integer subdivp1ind = vertsRemap.get(p1ind);
  //test for oddness on modulus 2 indicates thickness or thinness of edge 
  //on starting vertex index on edge point pair.
  int inc = 0;
  for (int i = subdivp1ind; i < SubdivPtsSize; i++){
    int ni = i+1;
    if (lastEdge && i == (SubdivPtsSize-1)){
      ni = 0;
    }
    ArrayList<Integer> ptpair = new ArrayList<Integer>();
    ptpair.add(i);
    ptpair.add(ni);
    if (inc % 2 == 0){
      SubdivPtPairtoThick.put(ptpair, true);
    }
    else{
      SubdivPtPairtoThick.put(ptpair, false);
    }
    inc += 1;
  }
}

void Subdivide(float frac, Polygon cpoly, ArrayList<Polygon> outPolys){
  ArrayList<Edge> cedges = cpoly.edges;
  ArrayList<PVector> subdivpts = new ArrayList<PVector>();
  HashMap<Integer,Integer> vertsRemap = new HashMap<Integer,Integer>();
  HashMap<ArrayList<Integer>,Boolean> SubdivPtPairtoThick = new HashMap<ArrayList<Integer>,Boolean>();
  for (int i = 0; i < cedges.size(); i++){
    Edge cedge = cedges.get(i);
    Boolean lastEdge = false;
    if (i == cedges.size()-1){
      lastEdge = true;
    }
    if (cedge.linear){
      if (cedge.thick){

        getNGonSubdivisionPoints(3,cedge,subdivpts,vertsRemap);
        cpoly.ptToSubdivPt = vertsRemap;
        cpoly.subdivPtToPtbuild();
        cpoly.subdivpts = subdivpts;
      }
      else{
        getNGonSubdivisionPoints(5,cedge,subdivpts,vertsRemap);
        cpoly.ptToSubdivPt = vertsRemap;
        cpoly.subdivPtToPtbuild();
        cpoly.subdivpts = subdivpts;
      }
    }
    else{
      if (cedge.thick){
        getCurveSubdivisionPoints(3,cedge,subdivpts, vertsRemap);
        cpoly.ptToSubdivPt = vertsRemap;
        cpoly.subdivPtToPtbuild();
        cpoly.subdivpts = subdivpts;
      }
      else{
        getCurveSubdivisionPoints(5,cedge,subdivpts, vertsRemap);

      }
    }

    buildThickEdgeDat(cedge, vertsRemap, subdivpts.size(), lastEdge,
               SubdivPtPairtoThick);
  }
  cpoly.ptToSubdivPt = vertsRemap;
  cpoly.subdivPtToPtbuild();
  println("ptToSubdivPt: ", cpoly.ptToSubdivPt);
  cpoly.subdivpts = subdivpts;
  println("Subdivpts: ", cpoly.subdivpts);
  cpoly.SubdivPtPairtoThick = SubdivPtPairtoThick; 
  //build arc/circle data for polygon
  ArrayList<Circle> centcircles = new ArrayList<Circle>();
  getCentroidCircles(frac, cpoly, centcircles);
  CircleMap circlemapout = new CircleMap();
  getSubdivPolyArcData(cpoly, subdivpts, vertsRemap, circlemapout);
  buildInteriorPolygons(cpoly, centcircles, circlemapout, outPolys);
  
}

void buildSubdivisions(ArrayList<ArrayList<Polygon>> polygonfam, int iter, float frac){
  int i = 0;
  while (i < iter){
    ArrayList<Polygon> polys = polygonfam.get(polygonfam.size()-1);
    ArrayList<Polygon> newpolys = new ArrayList<Polygon>();
    for (Polygon poly : polys){
      Subdivide(frac, poly, newpolys);
    }
    polygonfam.add(newpolys);
    i += 1;
  }
}

void initFirstPoly(PVector[] vertices, ArrayList<ArrayList<Polygon>> pfam){
  Boolean[] edgethick = new Boolean[vertices.length];
  Edge[] edges = new Edge[vertices.length];
  ArrayList<Polygon> polys = new ArrayList<Polygon>();
  if (vertices.length == 3){
    edgethick = new Boolean[] {false,false,false};
  }
  else if (vertices.length == 4){
    edgethick = new Boolean[] {false,true,false,true};
  }
  else{
    for(int i = 0; i < vertices.length; i++){
      edgethick[i] = true;
    }
  }
  for (int i = 0; i < vertices.length; i++){
    int ni = (i+1)%vertices.length;
     edges[i] = new Edge(vertices[i], vertices[ni], edgethick[i], i, ni);
  }
  buildInteriorPolygon(vertices, edgethick, edges, polys);
  pfam.add(polys);
}

void createNGon(ArrayList<PVector> pos, PShape s){
  s.beginShape();
  for (int i = 0; i < pos.size(); i++){
    s.vertex(pos.get(i).x, pos.get(i).y);
  }
  s.stroke(255);
  s.strokeWeight(.5);
  s.noFill();
  s.endShape(CLOSE);
}

void setup(){
  size(1080,720);
  getNGonPoints(NGONs1, RNG1, NG1pos);
  initFirstPoly(NG1pos, PFamily);
  //println(PFamily.get(0).get(0).vertices);
  buildSubdivisions(PFamily, SubDivLevel, frac);
  
  for (ArrayList<Polygon> polys : PFamily){
    for (Polygon poly : polys){
      PShape s = createShape();
      createNGon(poly.vertices, s);
      shapes.add(s);
    }
  }
}

void draw(){
  background(0);
  translate(1080.0/2.0, 720.0/2.0);
  for (PShape sh : shapes){
  shape(sh,0.0,0.0);
  }
  shape(shapes.get(0),0.0,0.0);
  int i = 0;
  for (ArrayList<Polygon> polys : PFamily){
    //println("Polygons: ", polys.size());
    for (Polygon poly: polys){
      if (i != PFamily.size()-1){
        //println("Polyarcs size: ", poly.arcs.size());
        for (Edge arci: poly.arcs){
          stroke(255);
          strokeWeight(.5);
          noFill();
          //println("angle1 : ", arci.angle1);
          //println("angle2 : ", arci.angle2);
          //println("radius : ", arci.genR);
          //println("center : ", arci.genCenter);
          //ellipse(arci.genCenter.x, arci.genCenter.y, 2.0*arci.genR, 2.0*arci.genR);
          //arc(arci.genCenter.x,arci.genCenter.y,arci.genR, arci.genR, arci.angle1, arci.angle2);
        }
      }
    }
    i += 1;
  }
  
}