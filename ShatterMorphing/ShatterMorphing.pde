//I want to build a morphing function that works between two animation states
//For this problem I have visualized at the moment some important characteristics when
// working between different geometries:
//1) The geometries are approximated by root lines in this case...neglecting curves
// although curves could be used.
//2) The initial and final configurations between two independent state geometries
// are such to have the same number approximating set of lines, for example, or given
// by whatever approximating geometry that is used.  
//3) The algorithm applies randomly applies one to one mapping between final and initial 
// geometry configurations given to either final and initial state represented 
// approximation.
//4) Describing an animation then between final initial states, in the interim
// is a process that should resemble destruction of the final state of one approximation
// state and reintegration of an initial state in so far as phase animation.
// How this then works:
// To do this we simply need to measure the global or the model() coordinate positions
// of the final state relative to the model() coordinate positions of the initial 
// configuration of the new state geometry.  In this case we form a mapping of  the basis 
// in describing all approximating elements in the final state, we'll call this
// S1 relative the initial state S2.  With the mapping then we measure the distances 
// on x and y axis between such elements and form a step division measure that 
// forms the basis of mapping such elements when redrawing an animation until having 
// reached configuration S2.  

//  Example:  A square S1 is given by a mapping of 30 line segments which is mapped into
// a triangle S2 which also has an approximating representation of 30 line segments.
// Randomly all elements in S1 are mapped to the elements of S2 and a metric is done
// given for the X and Y distances between S1 and S2 elements for each given mapping.
//  It is decided that there should be a total of 100 frames in step subdividing 
// such animation that completes the animation morphology from S1 to S2 so all such 
// state S1S2 metrics (for each axis respectively) are divided by 100, and a given 
// increment by such is added or subtracted in moving such element from position S1
// to S2.

//  Visually this might appear as a square shattering into line pieces and then 
// reassembling into a triangle.  

// Part of this problem in working by pure polygons is forming the basis of some
// subdivision of the polygon into parts.  The more obvious direct choice, for instance,
// is lines composed by the length of each side of the polygon direct, but there 
// is a given discrepancy, as in the example above,
//in the number of elements of S1 in this case given by 4 lines relative S2 which is
// 3 lines.  Which means that lower order polygons where the n the number of sides
// of such polygon are smaller relative to higher order value of n for a higher order
//polygon type.  This implies that immediately some subdivision process exists in
//finding the lcm (least common multiple) between the lower s1n type and a given
// higher order s2n type.

int NGONs1 = 3;
int NGONs2 = 7;
float RNG1 = 120;
float RNG2 = 100;
PVector[] NG1pos = new PVector[NGONs1];
PVector[] NG2pos = new PVector[NGONs2];
PShape s1;

// recursive implementation of gcd
int gcd(int p, int q) {
  if (q == 0) {return p;}
  else {return gcd(q, p % q);}
}

int lcm(int p, int q) {
  return int(abs(p*q)/gcd(p,q));
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

float linterp(float x, float p1x, float p1y, float p2x, float p2y){
  return p1y+(p2y-p1y)*(x-p1x)/(p2x-p1x);
}
//sdf is subdivision factor
void getNGonSubdivisionPoints(int sdf, PVector[] pos, PVector[] out){
  for (int i=0; i< pos.length; i++){
    int ni = i+1;
    if (i == pos.length -1){
      ni = 0;
    }
    PVector pos1 = pos[i];
    PVector pos2 = pos[ni];
    ////pick smallest x and smallest y
    //float sx = min(pos1.x,pos2.x);
    //float sy = min(pos1.y,pos2.y);
    ////get distance between positions
    //float dx = abs(pos1.x-pos2.x);
    //float dy = abs(pos1.y-pos2.y);
    ////get subdivision increment
    //float sbdx = dx/float(sdf);
    //float sbdy = dy/float(sdf);
    PVector p12 = pos1.sub(pos2);
    float mp12 = p12.mag();
    p12.normalize();
    mp12 /= float(sdf);
    p12 = PVector.mult(p12, mp12);
    for (int j=0; j < sdf; j++){
    }
    
  }
}

void setup(){
  size(1280, 720);
  PVector xy = new PVector(0.0,0.0,0.0);
  Polarcoord(2.0, 1.0, xy);
  println(xy);
  getNGonPoints(NGONs1, RNG1, NG1pos);
  int ng12lcm = lcm(NGONs1, NGONs2);
  PVector[] ng1posf = new PVector[ng12lcm];
  PVector[] ng2posf = new PVector[ng12lcm];
  println(ng12lcm);
  int s1f = ng12lcm/NGONs1;
  int s2f = ng12lcm/NGONs2;
  s1 = createShape();
  createNGon(NG1pos, s1);
  s1.setStroke(255);
  s1.setFill(0,0);
}

void createNGon(PVector[] pos, PShape s){
  s.beginShape();
  for (int i = 0; i < pos.length; i++){
    s.vertex(pos[i].x, pos[i].y);
  }
  s.strokeWeight(10);
  s.noFill();
  s.endShape(CLOSE);
}

void draw(){
  background(0);
  translate(1280.0/2.0, 720.0/2.0);
  //rotate(PI/4.0);
  shape(s1,0,0);
}