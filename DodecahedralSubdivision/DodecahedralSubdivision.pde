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
}

class Circle{
  PVector center;
  float radius;
  Circle(PVector Center, float Radius){
    center = Center;
    radius = Radius;
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

void getCentroidCircles(float frac, Polygon cpoly, ArrayList<Circle> out){
  if (cpoly.vertices.size() == 3){
    for (int i = 0; i < cpoly.vertices.size();i++){
      float ni = (i+1) % cpoly.vertices.size();
      //vector in the direction of the vertex from center
      PVector cv = PVector.sub(cpoly.center, cpoly.vertices.get(i));
      float dcv = cv.mag();
      dcv *= frac;
      cv.normalize();
      //CentroidCircle position is then given by 
      PVector centCircle = PVector.add(cpoly.vertices.get(i),PVector.mult(cv,dcv));
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
        PVector ptOnCurve = new PVector(0.0,0.0,0.0);
        d = distPointToCircle(edge.circle, centCircle, ptOnCurve);
      }
      float radius = d*(1-frac);
      out.add(new Circle(centCircle, radius);
    }
  }
  else if (cpoly.vertices.length == 4){
  }
  else if (cpoly.vertices.length >= 5){
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



void Subdivide(Polygon cpoly){
  ArrayList<Edge> cedges = cpoly.edges;
  ArrayList<PVector> out = new ArrayList<PVector>();
  HashMap<Integer,Integer> vertsRemap = new HashMap<Integer,Integer>();
  for (int i = 0; i < cedges.size(); i++){
    Edge cedge = cedges.get(i);
    if (cedge.linear){
      if (cedge.thick){
        
        getNGonSubdivisionPoints(3,cedge,out,vertsRemap);
      }
      else{
        getNGonSubdivisionPoints(5,cedge,out,vertsRemap);
      }
    }
    else{
      if (cedge.thick){
        getCurveSubdivisionPoints(3,cedge,out, vertsRemap);
      }
      else{
        getCurveSubdivisionPoints(5,cedge,out, vertsRemap);
      }
    }
  }
  
}