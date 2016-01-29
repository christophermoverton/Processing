import java.util.Map;
float atime = 1.0; //(animation time in seconds)
float frac = .30; //fractional size (recommend that this is < .5)
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
  
  Polygon(ArrayList<PVector> verts, PVector PCent){
    vertices = verts;
    center = PCent;
    arcs = new ArrayList<Edge>();
    arcincs = new ArrayList<Float>();
    spokenorms = new ArrayList<PVector>();
    sincvals = ArrayList<Float>();

    intpolypoints = new ArrayList<PVector>();
    p2s = new ArrayList<PVector>();
    edges = new ArrayList<Edge>();
    pointsToEdge = new HashMap<ArrayList<PVector>,Integer>();
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
  int pole2;  //used in special case edge mappings
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
    //linear = true;
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

void CircleCircleIntersection(Circle c1, Circle c2, PVector ipoint){
  PVector c1c2 = PVector.sub(c1.center, c2.center);
  float dc1c2 = c1c2.mag();
  float x = (dc1c2*dc1c2 - c1.radius*c1.radius + c2.radius*c2.radius)/2*dc1c2;
  float y = (4*pow(dc1c2,2)*pow(c2.radius,2) -pow(pow(dc1c2,2)-pow(c1.radius,2)+pow(c2.radius,2),2))/4*pow(dc1c2,2);
  ipoint.x = x;
  ipoint.y = y;
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

float distPointToLine(PVector p1, PVector p2, PVector p3){
  //p1 and p2 represent points on the line, and p3 is a point elsewhere.
  //see also http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  PVector p0p1 = PVector.sub(p0,p1);
  PVector p0p2 = PVector.sub(p0,p2);
  PVector p2p1 = PVector.sub(p2,p1);
  PVector p0p1Cp0p2 = PVector.cross(p0p1,p0p2);
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
  //println("Pos length: ", pos.length);
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
    //println("P12: ", p12);
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
  //println("Pos length: ", pos.length);
  
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
  //println("P12: ", p12);
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
  //println("Pos length: ", pos.length);
  
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
  //println("P12: ", p12);
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
      float ni = (i+1) % cpoly.vertices.size();
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
      int edgind = cpoly.pointsToEdge.get(edgepts);
      Edge edge = cpoly.edges.get(edgind);
      float d;
      if (edge.linear){
        d = distPointToLine(edge.p1,edge.p2,centCircle);
      }
      else{
        
        d = distPointToCircle(edge.circle, centCircle);
      }
      float radius = d*(1-frac);
      out.add(new Circle(centCircle, radius));
    }
  }
  else if (cpoly.vertices.length == 4){
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
        dcv *= frac;
        cv.normalize();
        //CentroidCircle position is then given by 
        PVector centCircle = PVector.add(cpoly.center,PVector.mult(cv,dcv));
        float d;
        if (cpoly.edges.get(i).linear){
          d = distPointToLine(cpoly.edges.get(i).p1,cpoly.edges.get(i).p2,centCircle);
        }
        else{
          
          d = distPointToCircle(cpoly.edges.get(i).circle, centCircle);
        }
        float radius = d*(1-frac);
        out.add(new Circle(centCircle, radius));  
      }
    }
  }
  else if (cpoly.vertices.length >= 5){
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
  return (p2.x-p1.x)*(y-p1.y)/(p2.y-p1.y)+p1.y;
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
      float nyi = cpoly.vertices.get((j+1)%vertices.size()).y;
      float py = pts.get(i).y;
      if ((yi <= py <= nyi)|| (yi >= py >= nyi)){
        test1 = true;
        edge.add(cpoly.vertices.get(j));
        edge.add(cpoly.vertices.get((j+1)%vertices.size()));
        edgepts.add(edge);
      }
    }
    if (edges.size() == 2){
      ArrayList<PVector> edg1 = edges.get(0);
      ArrayList<PVector> edg2 = edges.get(1);
      float x1 = xfromLine(edg1.get(0), edg1.get(1), pts.get(i).y);
      float x2 = xfromLine(edg2.get(0), edg2.get(1), pts.get(i).y);
      if ( min(x1,x2) <= pts.get(i).x <= max(x1,x2)){
        out.add(pts.get(i));
      }
    }
  }
}

void writeCircleMapData(HashMap<Integer,ArrayList<Integer>> vertToVertPair,
                        ArrayList<PVector> subdivpts,
                        CircleMap circlemapout, int i, boolean Pole){
           
    int ipi = (i - 1) % subdivpts.size();
    int iipi = (ipi - 2) % subdivpts.size();
    int pole1 = (ipi - 1) % subdivpts.size();
    int pole2 = (ipi - 1) % subdivpts.size();
    if (!Pole){
      iipi = (ipi - 1) % subdivpts.size();
      pole1 = (i + 1) % subdivpts.size();
      pole2 = (i - 4) % subdivpts.size();
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
    Edge edg = new Edge(P1, P2, c1.center, c1.radius, P1index, P2index);
    //repeat getting interior edge/arc/circle 1
    edg.interiornp1 = ipi;
    edg.interiornp2 = iipi;
    edg.pole1 = pole1;
    edg.pole2 = pole2;
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
      ni = i+4 % subdivpts.size();
      pi = i-4 % subdivpts.size();
      if (flaggedVerts.contains(ni)){
        int ipi = (i - 1) % subdivpts.size();
        int iipi = (ipi - 2) % subdivpts.size();
        if (vertToVertPair.containsKey(ipi) && vertToVertPair.containsKey(iipi)){
          writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, true);
        }
        continue;  
      }
    }
    else{
      int ini, ipi;
      ini = (i + 1)% subdivpts.size();
      ipi = (i - 1)% subdivpts.size();
      if (flaggedVerts.contains(ini) || flaggedVerts.contains(ipi)){
        ni = i+3 % subdivpts.size();
        pi = i-3 % subdivpts.size();
      }
      else{
        ni = i+4 % subdivpts.size();
        pi = i-4 % subdivpts.size();
        if (flaggedVerts.contains((ni-1)% subdivpts.size())){
          int ipi = (i - 1) % subdivpts.size();
          int iipi = (ipi - 1) % subdivpts.size();
          if (vertToVertPair.containsKey(ipi) && vertToVertPair.containsKey(iipi)){
            writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, true);
          }
          continue;  
        }
      }

    }
    //vert pair test to see if Arc already computed on previous pair set for a
    // given vertex iteration  if true for set inclusion, continue again.
    complpair.add(pi);
    complpair.add(i);
    if (complVertPairs.contains(complpair)){
      int ipi = (i - 1) % subdivpts.size();
      int iipi = (ipi - 1) % subdivpts.size();
      if (vertToVertPair.containsKey(ipi) && vertToVertPair.containsKey(iipi)){
        writeCircleMapData(vertToVertPair, subdivpts, circlemapout, i, false);
      }
      continue;  //we've already computed an arc for the present iterated vertex index i.
    }
    complpair = new ArrayList<Integer>();
    PVector Centerout = new PVector(0.0,0.0,0.0);
    getCircleCenter(subdivpts.get(i), subdivpts.get(ni), Centerout);
    PVector p1c = PVector.sub(Centerout,subdivpts.get(i));
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
  
}

void buildInteriorPolygons(ArrayList<Circle> centroidcircles, CircleMap circlemap){
  HashMap<Edge,ArrayList<Edge>> interioredges = circlemap.interiorEdges;
  for (Map.Entry me : interioredges.entrySet()) {
    Edge pedge = me.getKey();
    Edge iedge = me.getValue().get(0);
    Edge iedge2 = me.getValue().get(1);
    
  }
}

void Subdivide(float frac, Polygon cpoly){
  ArrayList<Edge> cedges = cpoly.edges;
  ArrayList<PVector> subdivpts = new ArrayList<PVector>();
  HashMap<Integer,Integer> vertsRemap = new HashMap<Integer,Integer>();
  for (int i = 0; i < cedges.size(); i++){
    Edge cedge = cedges.get(i);
    if (cedge.linear){
      if (cedge.thick){
        
        getNGonSubdivisionPoints(3,cedge,subdivpts,vertsRemap);
      }
      else{
        getNGonSubdivisionPoints(5,cedge,subdivpts,vertsRemap);
      }
    }
    else{
      if (cedge.thick){
        getCurveSubdivisionPoints(3,cedge,subdivpts, vertsRemap);
      }
      else{
        getCurveSubdivisionPoints(5,cedge,subdivpts, vertsRemap);
      }
    }
  }
  //build arc/circle data for polygon
  ArrayList<Circle> centcircles = new ArrayList<Circle>();
  getCentroidCircles(frac, cpoly, centcircles);
  
}