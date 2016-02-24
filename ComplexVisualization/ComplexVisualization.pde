import java.util.*;

Float Eradius = 40.0;
int NGONs1 = 5; //number of polygon sides
float RNG1 = 300;  //maximum radius of a polygon vertex from polygon center for the initializing polygon
PVector[] NG1pos = new PVector[NGONs1];
int SubDivLevel = 2;
float speed = 1.0; 
float time = 0.0;
float atime = 1.0; //(animation time in seconds)
float frac = .30; //fractional size realtionship of subdivision value greater than 0  and <= 1
PShape s1;
HashMap<Integer, ArrayList<Integer>> internal;
HashMap<PVector,Integer> pointToLabel = new HashMap<PVector,Integer>();
HashMap<Integer,PVector> labelToPoint = new HashMap<Integer,PVector>(); 
ArrayList<Polygon> polygons = new ArrayList<Polygon>();
ArrayList<ArrayList<Polygon>> PFamily = new ArrayList<ArrayList<Polygon>>();
ArrayList<ArrayList<Integer>> complexfamily;
HashMap<Integer, Float> external;
class Polygon{
  PVector[] vertices;
  PVector center;
  PVector[] spokes; //these are PVectors from Polygon center to subdivision point (B-A)
  float[] spokescalars;   // animation scalar increment value
  PVector[] spokenorms;  // normal vectors used in computing animation position... unit vectors in the direction of the spoke
  float[] sincvals;  //  current spoke scalar 'distance' values on animation cycle 
  ArrayList<PVector> intpolypoints; //interior polygon points for Pshape rendering...the same as p2s...no difference... deprecated at the moment...
  ArrayList<PVector> p2s;  //points originating on the smallest interior polygon boundaries
  
  Polygon(PVector[] verts, PVector PCent){
    vertices = verts;
    center = PCent;
    spokes = new PVector[vertices.length];
    spokescalars = new float[vertices.length];
    spokenorms = new PVector[vertices.length];
    sincvals = new float[vertices.length];
    for (int i = 0; i < sincvals.length; i++){
      sincvals[i] = 0.0;
    }
    intpolypoints = new ArrayList<PVector>();
    p2s = new ArrayList<PVector>();
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

//sdf is subdivision factor
void getNGonSubdivisionPoints(int sdf, PVector[] pos, PVector[] out){
  //println("Pos length: ", pos.length);
  for (int i=0; i< pos.length; i++){
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

void PolygonCentroid(PVector[] verts, PVector PCenter){
  float Cx = 0.0;
  float Cy = 0.0;
  float A = 0.0;
  for(int i=0; i<verts.length; i++){
    int ni = 0;
    if (i == verts.length-1){
      ni = 0;
    }
    else{
      ni = i+1;
    }
    Cx += (verts[i].x + verts[ni].x)*(verts[i].x*verts[ni].y - verts[ni].x*verts[i].y);
    Cy += (verts[i].y + verts[ni].y)*(verts[i].x*verts[ni].y - verts[ni].x*verts[i].y);
    A += verts[i].x*verts[ni].y - verts[ni].x*verts[i].y;
  }
  A*=.5;
  Cx*=1.0/(6.0*A);
  Cy*=1.0/(6.0*A);
  PCenter.x = Cx;
  PCenter.y = Cy;
}

ArrayList<ArrayList<Integer>> subdivide(ArrayList<ArrayList<Integer>> complexfamily,
               HashMap<Integer, ArrayList<Integer>> internal,
               HashMap<Integer, Float> external, Integer[] maxlabel){
  ArrayList<Polygon> cpolys = PFamily.get(PFamily.size()-1);
  ArrayList<Polygon> npolys = new ArrayList<Polygon>(); // next polygon subdiv level 
  // 5 point subdivision rule...we track labels only not points or positions.
  ArrayList<ArrayList<Integer>> nout = new ArrayList<ArrayList<Integer>>();
  HashMap<Integer[], Integer> subdivEdges = new HashMap<Integer[],Integer>();
  HashMap<ArrayList<Integer>, Integer> subdivEdges2 = new HashMap<ArrayList<Integer>, Integer>();
  Integer piter = 0;
  for (ArrayList<Integer> complex : complexfamily){
    Polygon cpoly = cpolys.get(piter);
    PVector[] cpolyverts = cpoly.vertices;
    PVector cpolycenter = cpoly.center;
    //pass to compute subdivision vertices
    PVector[] subdverts = new PVector[2*cpolyverts.length];
    //for (int j = 0; j < cpolyverts.length; j++){
    getNGonSubdivisionPoints(2, cpolyverts, subdverts);
    //}
    ArrayList<PVector> intPoly = new ArrayList<PVector>();
    ArrayList<PVector> isubdivPts = new ArrayList<PVector>();
    ArrayList<PVector> subdivPts = new ArrayList<PVector>();
    //next pass to get polygon vertices and center
    for (int j = 0; j < subdverts.length; j+=2){
      PVector nvert;  //next vertex (clockwise from existing vertex)
      PVector pvert;  // previous vertex (counter clockwise from existing vertex)
      PVector cvert = subdverts[j]; //current vertex
      if (j ==0){
        pvert = subdverts[subdverts.length-1];
      }
      else{
        pvert = subdverts[j-1];
      }
      if (j == subdverts.length-1){
        nvert = subdverts[0];
      }
      else{
        nvert = subdverts[j+1];
      }
      PVector cp = PVector.sub(pvert,cpolycenter);
      float cpmag = cp.mag();
      cpmag *= frac;
      cp.normalize();
      cp = PVector.mult(cp,cpmag);
      PVector p2 = PVector.add(cp, cpolycenter);
      PVector cn = PVector.sub(nvert,cpolycenter);
      float cnmag = cn.mag();
      cnmag *= frac;
      cn.normalize();
      cn = PVector.mult(cn,cnmag);
      PVector p1 = PVector.add(cn, cpolycenter);
      PVector[] npolyverts = new PVector[5];
      npolyverts[0] = pvert;
      npolyverts[1] = cvert;
      npolyverts[2] = nvert;
      npolyverts[3] = p1;
      npolyverts[4] = p2;
      intPoly.add(p1);
      isubdivPts.add(p2);
      subdivPts.add(nvert);
      PVector npolycenter = new PVector(0.0,0.0,0.0);
      PolygonCentroid(npolyverts, npolycenter);
      Polygon npoly = new Polygon(npolyverts,npolycenter);
      npolys.add(npoly);
    }
    PVector[] intPoly2 = new PVector[intPoly.size()];
    for (int j = 0; j<intPoly.size(); j++){
      intPoly2[j] = intPoly.get(j);
    }
    PVector npolycenter = new PVector(0.0,0.0,0.0);
    PolygonCentroid(intPoly2, npolycenter);
    Polygon npoly = new Polygon(intPoly2,npolycenter); 
    npolys.add(npoly);    

    //create new 5 complex
    ArrayList<Integer> inner5 = new ArrayList<Integer>();
    Integer[] isubdivlabels = {0,0,0,0,0};
    maxlabel[0] += 1;
    isubdivlabels[0] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[1] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[2] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[3] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[4] = maxlabel[0];
    maxlabel[0] += 1;
    Integer[] subdivlabels = {0,0,0,0,0};
    Boolean[] internallabel = {false,false,false,false,false};
    //get subdivlabel pts set to subdivlabel array
    //rewrite the original internal cycle.  Done once only...
    for (int i = 0; i < complex.size(); i++){
      Integer nlabel = complex.get((i+1)%complex.size());
      Integer ilabel = complex.get(i);
      Integer[] lpair = {nlabel,ilabel};
      Integer[] lpair2 = {ilabel, nlabel};
      ArrayList<Integer> lpair2al = new ArrayList<Integer>();
      Collections.addAll(lpair2al, lpair2);
      ArrayList<Integer> lpairal = new ArrayList<Integer>();
      Collections.addAll(lpairal, lpair);
      if (subdivEdges2.containsKey(lpairal)){
        subdivlabels[i] = subdivEdges2.get(lpairal);

      }
      else if (subdivEdges2.containsKey(lpair2al)){
        subdivlabels[i] = subdivEdges2.get(lpair2al);
      }
      else{
        maxlabel[0] += 1;
        subdivlabels[i] = maxlabel[0];
        subdivEdges.put(lpair2, maxlabel[0]);
        subdivEdges2.put(lpair2al, maxlabel[0]);
      }
      if (internal.containsKey(nlabel)){
        ArrayList<Integer> cycle = internal.get(nlabel);
        Integer ilcind = cycle.indexOf(ilabel);
        if (ilcind != -1){
          cycle.set(ilcind, subdivlabels[i]);
        }
      }
      if (internal.containsKey(ilabel)){
        ArrayList<Integer> cycle = internal.get(ilabel);
        Integer nlcind = cycle.indexOf(nlabel);
        if (nlcind != -1){
          cycle.set(nlcind, subdivlabels[i]);
        }
      }

      if (internal.containsKey(subdivlabels[i])){
        internal.get(subdivlabels[i]).add(isubdivlabels[i]);
      }
      else{
        if (internal.containsKey(ilabel)&&internal.containsKey(nlabel)){
          //set internal
          ArrayList<Integer> cycle = new ArrayList<Integer>();
          cycle.add(nlabel);
          cycle.add(isubdivlabels[i]);
          cycle.add(ilabel);
          internal.put(subdivlabels[i],cycle);
        }
        else{
          external.put(subdivlabels[i], Eradius);
        }
      }
      
    }
    //pass to write internal points
    for (int i = 0; i < isubdivlabels.length; i++){
      ArrayList<Integer> iComplex = new ArrayList<Integer>();
      Integer label1 = isubdivlabels[(i+1)%isubdivlabels.length];
      int pi = i - 1;
      if (i == 0){
        pi = isubdivlabels.length-1;
      }
      Integer label2 = isubdivlabels[pi];
      Integer label3 = subdivlabels[i];
      ArrayList<Integer> cycle = new ArrayList<Integer>();
      Collections.addAll(cycle, new Integer[] {label1,label2,label3});
      internal.put(isubdivlabels[i], cycle);
      label1 = isubdivlabels[i];
      label2 = subdivlabels[i];
      label3 = complex.get(i);
      Integer label4 = subdivlabels[(i+1)%subdivlabels.length];
      Integer label5 = isubdivlabels[(i+1)%isubdivlabels.length];
      ArrayList<Integer> ncomplex = new ArrayList<Integer>();
      Collections.addAll(ncomplex, new Integer[] {label1,label2,label3,label4,label5});
      inner5.add(isubdivlabels[i]);
      nout.add(ncomplex);
      int j = i -1;
      if (j == -1){
        j = 4;
      }
      pointToLabel.put(isubdivPts.get((i)%5), isubdivlabels[i]);
      labelToPoint.put(isubdivlabels[i], isubdivPts.get((i)%5));
      pointToLabel.put(subdivPts.get(j), subdivlabels[i]);
      labelToPoint.put(subdivlabels[i], subdivPts.get(j));
    }
    nout.add(inner5);
    //Integer label6 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label7 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label8 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label9 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label10 = maxlabel[0];
    piter += 1;
  }
  PFamily.add(npolys);
  return nout;
}

void triangulateComplex(ArrayList<ArrayList<Integer>> complexfamily,
                        HashMap<Integer, ArrayList<Integer>> internal, 
                        Integer[] maxlabel){
    Integer tlabelSize = complexfamily.size()+1;
    ArrayList<Integer> newLabels = new ArrayList<Integer>();
    for (int i = 0; i < complexfamily.size(); i++){
      maxlabel[0] += 1;
      newLabels.add(maxlabel[0]);
      ArrayList<Integer> cycle = new ArrayList<Integer>(complexfamily.get(i));
      internal.put(maxlabel[0], cycle);
      for (int j = 0; j < cycle.size(); j++){
        Integer v = cycle.get(j);
        if (internal.containsKey(v)){
          Integer ni1 = j-1;
          if (j == 0){
             ni1 = cycle.size()-1;
          }
          Integer ni2 = (j+1)%cycle.size();
          Integer nlabel1 = cycle.get(ni1);
          Integer nlabel2 = cycle.get(ni2);
          ArrayList<Integer> vpetal = internal.get(v);
          //println(nlabel1);
          //println(nlabel2);
          //println(v);
          //println(vpetal);
          Integer vpind1 = vpetal.indexOf(nlabel1);
          Integer vpind2 = vpetal.indexOf(nlabel2);
          boolean t1 = vpind1 == 0 || vpind1 == (vpetal.size()-1);
          boolean t2 = vpind2 == 0 || vpind2 == (vpetal.size()-1);
          if (t1 && t2){
            vpetal.add(maxlabel[0]);
          }
          else{
            Integer insIndex = max(vpind1,vpind2);
            vpetal.add(insIndex, maxlabel[0]);
          }
        }
      }
    }
}

void setup() {
  size(1080,720);
  getNGonPoints(NGONs1, RNG1, NG1pos);
  PVector PCenter = new PVector(0.0,0.0,0.0);
  PolygonCentroid(NG1pos, PCenter);
  Polygon poly = new Polygon(NG1pos,PCenter);
  polygons.add(poly);
  PFamily.add(polygons);
  internal = new HashMap<Integer, ArrayList<Integer>>();
  ArrayList<Integer> ivgraph = new ArrayList<Integer>();
  Integer[] ivgarr = {1,2,3,4,5};
  Collections.addAll(ivgraph,ivgarr);
  println(ivgraph);
  //internal.put(6, ivgraph);
  external = new HashMap<Integer,Float>();
  external.put(1,Eradius);
  external.put(2,Eradius);
  external.put(3,Eradius);
  external.put(4,Eradius);
  external.put(5,Eradius);
  for (int i = 0; i < 5; i++){
    pointToLabel.put(NG1pos[i], i+1);
    labelToPoint.put(i+1,NG1pos[i]);
  }
  complexfamily = new ArrayList<ArrayList<Integer>>();
  complexfamily.add(ivgraph);
  Integer[] maxlabel = {5};
  complexfamily = subdivide(complexfamily, internal, external, maxlabel);
  //complexfamily = subdivide(complexfamily, internal, external, maxlabel);
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

void draw(){
  background(0);
  //stroke(255);
  //text(complexfamily.toString(),10.0,10.0);
  translate(1080.0/2.0, 720.0/2.0);
  for (Map.Entry<PVector,Integer> ptlabel : pointToLabel.entrySet()){
    PVector pt = ptlabel.getKey();
    noFill();
    stroke(255);
    strokeWeight(10);
    point(pt.x,pt.y);
   //for (
    //ellipse(circv.center.x, circv.center.y, 2*circv.radius, 2*circv.radius);
  }
  strokeWeight(1);
  for (Map.Entry<Integer,ArrayList<Integer>> intern: internal.entrySet()){
    PVector pt1 = labelToPoint.get(intern.getKey());
    for (Integer label2 : intern.getValue()){
      PVector pt2 = labelToPoint.get(label2);
      stroke(50,50,255);
      line(pt1.x,pt1.y,pt2.x,pt2.y);

    }
    fill(180,180,255);
    textSize(12);
    text(intern.getKey().toString(),pt1.x+5.0,pt1.y+10.0);
  }
  if (time > 400){
    for (ArrayList<Integer> complex : complexfamily){
      ArrayList<PVector> pos = new ArrayList<PVector>();
      for (Integer poslabel : complex){
        pos.add(labelToPoint.get(poslabel));
        PVector pt1 = labelToPoint.get(poslabel);
        fill(180,180,255);
        textSize(12);
        text(poslabel.toString(),pt1.x+5.0,pt1.y+10.0);        
      }
      PShape s = createShape();
      createNGon(pos, s);
      shape(s,0.0,0.0);
      
    }
  }
  time += speed;
}