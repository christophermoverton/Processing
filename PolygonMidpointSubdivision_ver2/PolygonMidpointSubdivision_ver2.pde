//Polygon MidPoint Subdivision
//
int NGONs1 = 3;
float RNG1 = 300;
PVector[] NG1pos = new PVector[NGONs1];
int SubDivLevel = 2;
float speed = 1.0; 
float time = 0.0;
float atime = 2.0; //(animation time in seconds)
float frac = .40; //fractional size realtionship of subdivision
PShape s1;

ArrayList<Polygon> polygons = new ArrayList<Polygon>();
ArrayList<ArrayList<Polygon>> PFamily = new ArrayList<ArrayList<Polygon>>();
class Polygon{
  PVector[] vertices;
  PVector center;
  PVector[] spokes; //these are PVectors from Polygon center to subdivision point (B-A)
  float[] spokescalars;   // animation scalar increment value
  PVector[] spokenorms;  // normal vectors used in computing animation position 
  float[] sincvals;  //  current spoke increment scalar values on animation cycle 
  ArrayList<PVector> intpolypoints; //interior polygon points for Pshape rendering
  ArrayList<PVector> p2s;  //points originating on the interior polygon boundaries
  
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
 
}

void getNGonPoints(int nside, float radius, PVector[] pos){
  float ainc = 2.0*PI/float(nside);
  float a = -PI/2.0;
  if (nside % 2 == 0){
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

void Subdivide(ArrayList<ArrayList<Polygon>> pfam, float atime, float frac){
  ArrayList<Polygon> npolys = new ArrayList<Polygon>(); // next polygon subdiv level 
  ArrayList<Polygon> cpolys = pfam.get(pfam.size()-1); //current polygons subdiv level
  for (int i = 0; i < cpolys.size(); i++){
    Polygon cpoly = cpolys.get(i);
    PVector[] cpolyverts = cpoly.vertices;
    PVector cpolycenter = cpoly.center;
    //pass to compute subdivision vertices
    PVector[] subdverts = new PVector[2*cpolyverts.length];
    //for (int j = 0; j < cpolyverts.length; j++){
    getNGonSubdivisionPoints(2, cpolyverts, subdverts);
    //}
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
      PVector npolycenter = new PVector(0.0,0.0,0.0);
      PolygonCentroid(npolyverts, npolycenter);
      Polygon npoly = new Polygon(npolyverts,npolycenter);
      npolys.add(npoly);
    }
    //pass to add existing polygon spoke data
    for (int j = 1; j < subdverts.length; j+=2){
      PVector pos2 = subdverts[j];
      PVector pos1 = cpolycenter;
      PVector p12 = PVector.sub(pos2,pos1);
      int k = (j-1)/2;
      cpoly.spokes[k] = p12;
      float mp12 = p12.mag();
      mp12 /= (atime*60); //seconds * 60 frames per seconds yields inc scalar per frame
      cpoly.spokescalars[k]=mp12;
      p12.normalize();
      cpoly.spokenorms[k] = p12;
    }
  }
  println("Number of polys: ", npolys.size());
  pfam.add(npolys);
}

void createNGon(PVector[] pos, PShape s){
  s.beginShape();
  for (int i = 0; i < pos.length; i++){
    s.vertex(pos[i].x, pos[i].y);
  }
  s.stroke(255);
  s.strokeWeight(2);
  s.noFill();
  s.endShape(CLOSE);
}

void createNGon(ArrayList<PVector> pos, PShape s){
  s.beginShape();
  for (int i = 0; i < pos.size(); i++){
    s.vertex(pos.get(i).x, pos.get(i).y);
  }
  s.stroke(255);
  s.strokeWeight(2);
  s.noFill();
  s.endShape(CLOSE);
}

void setup() {
  size(1080,720);
  getNGonPoints(NGONs1, RNG1, NG1pos);
  s1 = createShape();
  createNGon(NG1pos, s1);
  PVector PCenter = new PVector(0.0,0.0,0.0);
  PolygonCentroid(NG1pos, PCenter);
  Polygon poly = new Polygon(NG1pos,PCenter);
  polygons.add(poly);
  PFamily.add(polygons);
  int i = 0;
  while (i < SubDivLevel){
    Subdivide(PFamily, atime, frac);
    i++;
  }
}


boolean nextLevel = false;
boolean continueAnim = true;
boolean recordShapes = false;
int clevel = 0;
float maxscalar = 0.0; //max scalar
float pscalar = 0.0; //present scalar
ArrayList<Polygon> polys;
ArrayList<PVector[]> linelist = new ArrayList<PVector[]>();
ArrayList<PShape> shapeslist = new ArrayList<PShape>();
void draw(){
  background(0);
  translate(1080.0/2.0, 720.0/2.0);
  shape(s1,0.0,0.0);
  for (int i = 0; i < linelist.size(); i++){
    PVector[] ptpair = linelist.get(i);
    PVector p1 = ptpair[0];
    PVector p2 = ptpair[1];
    line(p1.x,p1.y,p2.x,p2.y);
    stroke(255);
    strokeWeight(2);      
  }
  if (continueAnim){

    ArrayList<Polygon> polys = PFamily.get(clevel);
    for (int i = 0; i < polys.size(); i++){
      Polygon poly = polys.get(i);
      float[] spscalars = poly.spokescalars;
      for (int j = 0; j < spscalars.length; j++){
        if (poly.sincvals[j] == 0.0){
          poly.sincvals[j] = spscalars[j]*(atime*60)*frac;
          float spscalar2 = poly.sincvals[j];
          PVector spvector2 = PVector.mult(poly.spokenorms[j], spscalar2);
          PVector p2 = PVector.add(spvector2,poly.center);
          poly.p2s.add(p2);
          recordShapes = true;
        }
        float spscalar = poly.sincvals[j];
        PVector spvector = PVector.mult(poly.spokenorms[j], spscalar);
        PVector p1 = PVector.add(spvector, poly.center);
        
        line(p1.x,p1.y,poly.p2s.get(j).x,poly.p2s.get(j).y);
        stroke(255);
        strokeWeight(2);
        poly.sincvals[j] += spscalars[j];
        if (poly.sincvals[j] >= spscalars[j]*(atime*60)){
          nextLevel = true;
        }
        //if (nextLevel){
        //  PVector[] ptpair = new PVector[2];
        //  ptpair[0] = p1;
        //  ptpair[1] = poly.center;
        //  linelist.add(ptpair);
        //}
        
      }
      if (recordShapes){
      //add interior polygon shape
        PShape s2 = createShape();
        createNGon(poly.p2s, s2);
        recordShapes = false;
        shapeslist.add(s2);
      }
    }

    if(nextLevel){
      for (int i = 0; i < polys.size(); i++){
        Polygon poly = polys.get(i);
        float[] spscalars = poly.spokescalars;
        for (int j = 0; j < spscalars.length; j++){
          //float spscalar = poly.sincvals[j];
          PVector spvector = PVector.mult(poly.spokenorms[j],spscalars[j]*(atime*60));
          PVector p1 = PVector.add(spvector, poly.center);
          PVector[] ptpair = new PVector[2];
          ptpair[0] = p1;
          ptpair[1] = poly.p2s.get(j);
          linelist.add(ptpair);
          //line(p1.x,p1.y,poly.center.x,poly.center.y);
          //stroke(255);
          //strokeWeight(2);
          //poly.sincvals[j] += spscalars[j];
        }
      }
      clevel += 1;
      if (clevel > PFamily.size()-2){
        continueAnim = false;
      }
      nextLevel = false;
    }
  }
  for (int i = 0; i < shapeslist.size(); i++){
    shape(shapeslist.get(i),0.0,0.0);
  }
}