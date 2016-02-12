import java.util.Map;
import java.util.Collections;
import java.util.Collection;
int NGONs1 = 3; //number of polygon sides
float RNG1 = 300;  //maximum radius of a polygon vertex from polygon center for the initializing polygon
float atime = 1.0; //(animation time in seconds)
float frac = .2; //fractional size (recommend that this is < .5, .2 is optimal)
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
}

class Circle{
  PVector center;
  float radius;
  Circle(PVector Center, float Radius){
    center = Center;
    radius = Radius;
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

float distPointToCircle(Circle ccircle, PVector p){
  PVector pc = PVector.sub(p,ccircle.center);
  pc.normalize();
  pc = PVector.mult(pc,ccircle.radius);
  PVector k = PVector.add(ccircle.center,pc);
  PVector kp = PVector.sub(k,p);
  return kp.mag();
}

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
        d = distPointToLine(edge.p1,edge.p2,centCircle);
      }
      else{
        
        d = distPointToCircle(edge.circle, centCircle);
      }
      float radius = d*(1-frac) + d*frac*frac;
      out.add(new Circle(centCircle, radius));
    }
  }
  //else if (cpoly.vertices.size() == 4){
  //  for (int i = 0; i < cpoly.edges.size();i++){
  //    if (cpoly.edges.get(i).thick){
  //      ArrayList<PVector> pts = new ArrayList<PVector>();
  //      if (cpoly.edges.get(i).linear){
  //        getNGonSubdivisionPoints(2, cpoly.edges.get(i), pts);
  //      }
  //      else{
  //        getCurveSubdivisionPoints(2, cpoly.edges.get(i), pts);
  //      }
        
  //      PVector cv = PVector.sub(pts.get(1),cpoly.center);
  //      float dcv = cv.mag();
  //      dcv *= frac;
  //      cv.normalize();
  //      //CentroidCircle position is then given by 
  //      PVector centCircle = PVector.add(cpoly.center,PVector.mult(cv,dcv));
  //      float d;
  //      if (cpoly.edges.get(i).linear){
  //        d = distPointToLine(cpoly.edges.get(i).p1,cpoly.edges.get(i).p2,centCircle);
  //      }
  //      else{
          
  //        d = distPointToCircle(cpoly.edges.get(i).circle, centCircle);
  //      }
  //      float radius = d*(1-frac);
  //      out.add(new Circle(centCircle, radius));  
  //    }
  //  }
  //}
  else if (cpoly.vertices.size() >= 5){
    float d;
    if (cpoly.edges.get(0).linear){
      d = distPointToLine(cpoly.edges.get(0).p1,cpoly.edges.get(0).p2,cpoly.center);
    }
    else{
      
      d = distPointToCircle(cpoly.edges.get(0).circle, cpoly.center);
    }
    float radius = d*(1-frac);
    out.add(new Circle(cpoly.center, radius)); 
  }
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
    println("parent gencenter: ", parent.genCenter);
    println("pt2", pt2);
    PVector CToP2 = PVector.sub(pt2,parent.genCenter);
    float a1 = CToP1.heading();
    float a2 = CToP2.heading();
    Edge out = new Edge(pt1, pt2, edgeThick, parent.genCenter, 
                   parent.genR, a1, a2, np1index, np2index);
    outs.add(out);
  }
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

ArrayList<Circle> circles = new ArrayList<Circle>();

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
  background(0);
  getNGonPoints(NGONs1, RNG1, NG1pos);
  initFirstPoly(NG1pos, PFamily);
  for (ArrayList<Polygon> polys : PFamily){
    for(Polygon poly : polys){
       getCentroidCircles(frac, poly, circles);
       PShape s = createShape();
       createNGon(poly.vertices, s);
       shapes.add(s);
    }
  }
}

void draw(){
  fill(0,2);
  translate(1080.0/2.0, 720.0/2.0);
  for (PShape sh : shapes){
   shape(sh,0.0,0.0);
  }
  for (Circle circ: circles){
    stroke(255);
    strokeWeight(.5);
    noFill();
    ellipse(circ.center.x, circ.center.y, circ.radius, circ.radius);
  }
}