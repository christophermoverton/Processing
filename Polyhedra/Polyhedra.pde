
///***************************************************************************
// * Quaternion class written by BlackAxe / Kolor aka Laurent Schmalen in 1997
// * Translated to Java(with Processing) by RangerMauve in 2012
// * this class is freeware. you are fully allowed to use this class in non-
// * commercial products. Use in commercial environment is strictly prohibited
// */

//class Quaternion {
//  float W, X, Y, Z;      // components of a quaternion

//  // default constructor
//   Quaternion() {
//    W = 1.0;
//    X = 0.0;
//    Y = 0.0;
//    Z = 0.0;
//  }

//  // initialized constructor

//  Quaternion(float w, float x, float y, float z) {
//    W = w;
//    X = x;
//    Y = y;
//    Z = z;
//  }

//  // quaternion multiplication
//  Quaternion qmult (Quaternion q) {
//    println("Checking multiplication on w:");
//    println(W,q.W, X,q.X, Y,q.Y, Z,q.Z);
//    println(W*q.W, X*q.X, Y*q.Y, Z*q.Z);
//    float w = W*q.W - (X*q.X + Y*q.Y + Z*q.Z);

//    float x = W*q.X + q.W*X + Y*q.Z - Z*q.Y;
//    float y = W*q.Y + q.W*Y + Z*q.X - X*q.Z;
//    float z = W*q.Z + q.W*Z + X*q.Y - Y*q.X;
//    println("Inside multiplication: ");
//    println(w,x,y,z);
//    W = w;
//    X = x;
//    Y = y;
//    Z = z;
//    return this;
//  }

//  // conjugates the quaternion
//  Quaternion conjugate () {
//    X = -X;
//    Y = -Y;
//    Z = -Z;
//    return this;
//  }

//  // inverts the quaternion
//  Quaternion reciprical () {
//    float norme = sqrt(W*W + X*X + Y*Y + Z*Z);
//    if (norme == 0.0)
//      norme = 1.0;

//    float recip = 1.0 / norme;

//    W =  W * recip;
//    X = -X * recip;
//    Y = -Y * recip;
//    Z = -Z * recip;

//    return this;
//  }

//  // sets to unit quaternion
//  Quaternion qnormalize() {
//    float norme = sqrt(W*W + X*X + Y*Y + Z*Z);
//    if (norme == 0.0)
//    {
//      W = 1.0; 
//      X = Y = Z = 0.0;
//    }
//    else
//    {
//      float recip = 1.0/norme;

//      W *= recip;
//      X *= recip;
//      Y *= recip;
//      Z *= recip;
//    }
//    return this;
//  }

//  // Makes quaternion from axis
//  Quaternion fromAxis(float Angle, float x, float y, float z) { 
//    float omega, s, c;
//    int i;

//    s = sqrt(x*x + y*y + z*z);

//    if (abs(s) > Float.MIN_VALUE)
//    {
//      c = 1.0/s;

//      x *= c;
//      y *= c;
//      z *= c;

//      omega = -0.5f * Angle;
//      s = (float)sin(omega);

//      X = s*x;
//      Y = s*y;
//      Z = s*z;
//      W = (float)cos(omega);
//    }
//    else
//    {
//      X = Y = 0.0f;
//      Z = 0.0f;
//      W = 1.0f;
//    }
//    qnormalize();
//    return this;
//  }

//  Quaternion fromAxis(float Angle, PVector axis) {
//    return this.fromAxis(Angle, axis.x, axis.y, axis.z);
//  }

//  // Rotates towards other quaternion
//  void slerp(Quaternion a, Quaternion b, float t)
//  {
//    float omega, cosom, sinom, sclp, sclq;
//    int i;


//    cosom = a.X*b.X + a.Y*b.Y + a.Z*b.Z + a.W*b.W;


//    if ((1.0f+cosom) > Float.MIN_VALUE)
//    {
//      if ((1.0f-cosom) > Float.MIN_VALUE)
//      {
//        omega = acos(cosom);
//        sinom = sin(omega);
//        sclp = sin((1.0f-t)*omega) / sinom;
//        sclq = sin(t*omega) / sinom;
//      }
//      else
//      {
//        sclp = 1.0f - t;
//        sclq = t;
//      }

//      X = sclp*a.X + sclq*b.X;
//      Y = sclp*a.Y + sclq*b.Y;
//      Z = sclp*a.Z + sclq*b.Z;
//      W = sclp*a.W + sclq*b.W;
//    }
//    else
//    {
//      X =-a.Y;
//      Y = a.X;
//      Z =-a.W;
//      W = a.Z;

//      sclp = sin((1.0f-t) * PI * 0.5);
//      sclq = sin(t * PI * 0.5);

//      X = sclp*a.X + sclq*b.X;
//      Y = sclp*a.Y + sclq*b.Y;
//      Z = sclp*a.Z + sclq*b.Z;
//    }
//  }

//  Quaternion exp()
//  {                               
//    float Mul;
//    float Length = sqrt(X*X + Y*Y + Z*Z);

//    if (Length > 1.0e-4)
//      Mul = sin(Length)/Length;
//    else
//      Mul = 1.0;

//    W = cos(Length);

//    X *= Mul;
//    Y *= Mul;
//    Z *= Mul; 

//    return this;
//  }

//  Quaternion log()
//  {
//    float Length;

//    Length = sqrt(X*X + Y*Y + Z*Z);
//    Length = atan(Length/W);

//    W = 0.0;

//    X *= Length;
//    Y *= Length;
//    Z *= Length;

//    return this;
//  }
//}
class Quaternion {
  float W,X,Y,Z;
  Quaternion(float w, float x, float y, float z){
    W = w;
    X = x;
    Y = y;
    Z = z;
  }
  Quaternion conjugate(){
    float x = -X;
    float y = -Y;
    float z = -Z;
    return new Quaternion(W,x,y,z);
  }
  void qnorm(){
    float qmag = pow(W*W+X*X+Y*Y+Z*Z,.5);
    W *= 1/qmag;
    X *= 1/qmag;
    Y *= 1/qmag;
    Z *= 1/qmag;
   
  }
}


Quaternion qmult2 (Quaternion p, Quaternion q) {
  //println("Checking multiplication on w:");
  //println(p.W,q.W, p.X,q.X, p.Y,q.Y, p.Z,q.Z);
  //println(p.W*q.W, p.X*q.X, p.Y*q.Y, p.Z*q.Z);
  float w = p.W*q.W - (p.X*q.X + p.Y*q.Y + p.Z*q.Z);
  float x = p.W*q.X + q.W*p.X + p.Y*q.Z - p.Z*q.Y;
  float y = p.W*q.Y + q.W*p.Y + p.Z*q.X - p.X*q.Z;
  float z = p.W*q.Z + q.W*p.Z + p.X*q.Y - p.Y*q.X;
  
  //float w = p.W*q.W - (p.X*q.X + p.Y*q.Y + p.Z*q.Z);

  //float x = p.W*q.X + q.W*p.X + p.Y*q.Z - p.Z*q.Y;
  //float y = p.W*q.Y + q.W*p.Y + p.Z*q.X - p.X*q.Z;
  //float z = p.W*q.Z + q.W*p.Z + p.X*q.Y - p.Y*q.X;
  //println("Inside multiplication: ");
  //println(w,x,y,z);
  return new Quaternion(w,x,y,z);
  
}
//Example of rotating PVector about a directional PVector
PVector Qrotate(PVector v, PVector r, float a) {
  Quaternion Q1 = new Quaternion(0.0, v.x, v.y, v.z);
  Quaternion Q2 = new Quaternion(cos(a / 2.0), r.x * sin(a / 2.0), r.y * sin(a / 2.0), r.z * sin(a / 2.0));
  //println(Q2.W,Q2.X, Q2.Y, Q2.Z);
  Quaternion Q2C = Q2.conjugate();
  //println(Q2C.W, Q2C.X, Q2C.Y, Q2C.Z);
  Quaternion q2q1 = qmult2 (Q2, Q1);
  //Quaternion q2q1 = Q2.qmult(Q1);
  //println(q2q1.W, q2q1.X, q2q1.Y, q2q1.Z);
  Quaternion q2q1q2inv = qmult2 (q2q1, Q2C);
  //Quaternion q2q1q2inv = q2q1.qmult(Q2C);
  //println(q2q1q2inv.W, q2q1q2inv.X, q2q1q2inv.Y, q2q1q2inv.Z);
  //Quaternion Q3 = (Q2.qmult(Q1)).qmult(Q2.conjugate());
 // println(Q3.W,Q3.X,Q3.Y,Q3.Z);
  return new PVector(q2q1q2inv.X, q2q1q2inv.Y, q2q1q2inv.Z);
}

PVector sphericalToC(PVector p){
  //p in cartesian
  float r = p.x;
  float theta = p.y;
  float phi = p.z;
  return new PVector(r*cos(theta)*sin(phi), r*sin(theta)*sin(phi), r*cos(phi));
}

PVector cToSpherical(PVector p){
  //p in spherical coordinates system
  float r = pow(p.x*p.x + p.y*p.y+ p.z*p.z,.5);
  float theta = atan2(p.y,p.x);
  float phi = acos(p.z/r);
  return new PVector(r,theta,phi);
}

PVector polarCoord(float r, float theta){
  return new PVector(r*cos(theta), r*sin(theta));
}

class Disk{
  PVector origin; //spherical coordinates origin  
  float r;
  PVector dpoint;  //on the disk in spherical coordinates
  Disk(){
    origin = new PVector(0.0,0.0,0.0);
    dpoint = new PVector(0.0,0.0,0.0);
    r = 0.0;
  }
  Disk(float R, PVector Origin, PVector Dpoint){
    r = R;
    origin = Origin;
    dpoint = Dpoint;
  }
}

Disk getDisk(float deltaphi, PVector sp){
  //sp is a spherical coordinate point
  PVector pt1sp = new PVector(sp.x, sp.y, sp.z+deltaphi);
  PVector pt2sp = new PVector(sp.x, sp.y, sp.z-deltaphi);
  PVector pt1 = sphericalToC(pt1sp);//in cartesian
  PVector pt2 = sphericalToC(pt2sp);
  PVector pt1pt2 = pt1.sub(pt2);
  float pt1pt2mag = pt1pt2.mag();
  pt1pt2mag *= .5;
  pt1pt2.normalize();
  pt1pt2.mult(pt1pt2mag);
  PVector dcenter = pt1pt2.add(pt2);
  println("Dcenter: ",dcenter);
  println("Disk Radius: ", pt1pt2mag);
  PVector dorigin = cToSpherical(dcenter); //in spherical coordinates
  return new Disk(pt1pt2mag, dorigin, pt1sp);
}
ArrayList<Disk> disks = new ArrayList<Disk>();
ArrayList<ArrayList<PVector>> diskspoints = new ArrayList<ArrayList<PVector>>();

float phidiv = 5;
float thetadiv = 5;
float thetadiv2 = 30.0;
float deltaphi = PI/6.0;
float sphereRad = 200.0;
void setup(){
  PVector zvec = new PVector(0.0,0.0,1.0);
  PVector xvec = new PVector(1.0,0.0,0.0);
 
  //println(Qrotate(xvec,zvec,theta));
  //get disks
  float theta = 0.0;
  float phi = -PI;
  float thetainc = 2.0*PI/thetadiv;
  float theta2inc = 2.0*PI/thetadiv2;
  float phiinc = PI/phidiv;
  for (int i = 0; i < phidiv; i++ ){
    theta = 0.0;
    for (int j = 0; j < thetadiv; j++){
      //PVector sp = sphericalToC(new PVector(sphereRad,theta, phi));
      disks.add( getDisk(deltaphi, new PVector(sphereRad,theta, phi)));
      theta += thetainc;
    }
    phi += phiinc;
  }
  
  //generate diskpoints
  for( Disk disk : disks){
    //theta = disk.origin.y;
    //phi = disk.origin.z;
    ArrayList<PVector> diskpoints = new ArrayList<PVector>();
    PVector corigin = sphericalToC(disk.origin);
    PVector pt1 = sphericalToC(disk.dpoint);
    PVector p1toCorigin = pt1.sub(corigin);
    PVector ncorigin = new PVector(corigin.x,corigin.y,corigin.z);
    ncorigin.normalize();
    p1toCorigin.normalize();
    theta = 0.0;
    for (int i = 0; i < thetadiv2; i++){
      PVector point = Qrotate(p1toCorigin, ncorigin, theta);
      point = point.mult(disk.r);
      point = point.add(corigin);
      theta += theta2inc;
      diskpoints.add(point);
      //println(point);
    }
    diskspoints.add(diskpoints);
  }
  size(1920,1080, P3D);
  background(0);
  lights();
}

float t = 0.0;
boolean DrawSphere = true;
void draw(){
  fill(0,80);
  noStroke();
  beginShape();
  vertex(0,0,-200);
  vertex(1920,0, -200);
  vertex(1920,1080,-200);
  vertex(0,1080,-200);
  endShape();
  //rect(0,0, 1920,1080);
  translate(1920/2, 1080/2, 0);
  rotateY(t);
  noFill();
  stroke(200);
  strokeWeight(.3);
  sphereDetail(20);
  if (DrawSphere){
    sphere(200);
  }
  int fillc = 0;
  for (ArrayList<PVector> diskpoints: diskspoints){
    fill(fillc%255,60);
    beginShape();
    
    for (PVector point: diskpoints){
      vertex(point.x,point.y,point.z);
      //println(point);
    }
    endShape();
    fillc += 10;
  }
  t += .001;
}