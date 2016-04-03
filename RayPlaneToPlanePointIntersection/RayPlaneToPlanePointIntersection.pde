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

class rayTest{
  ArrayList<Float> rayscalars;
  ArrayList<Boolean> intersection;
  ArrayList<Float> rayscalarincs;
  ArrayList<Float> lastDistances;
}

ArrayList<PVector> generateRays(Disk disk){
    PVector corigin = sphericalToC(disk.origin);
    PVector pt1 = sphericalToC(disk.dpoint);
    PVector p1toCorigin = pt1.sub(corigin);
    PVector ncorigin = new PVector(corigin.x,corigin.y,corigin.z);
    ncorigin.normalize();
    p1toCorigin.normalize();
    ArrayList<PVector> rays = new ArrayList<PVector>();
    rays.add(p1toCorigin);
    float theta = PI/2.0;
    for (int i = 0; i < 3; i++){
      PVector point = Qrotate(p1toCorigin, ncorigin, theta);
      rays.add(point);
      theta += PI/2.0;
    }
    return rays;
}

Float PointToPlaneDist(Disk disk, PVector ray){
    PVector corigin = sphericalToC(disk.origin);
    PVector pt1 = sphericalToC(disk.dpoint);
    //PVector p1toCorigin = pt1.sub(corigin);
    PVector ncorigin = new PVector(corigin.x,corigin.y,corigin.z);
    ncorigin.normalize();
    PVector raytopt1 = ray.sub(pt1);
    return raytopt1.dot(ncorigin);
}

ArrayList<PVector> rayTester(Disk disk,  Disk disk2, 
                             ArrayList<PVector> rays, rayTest raytest){
  int i = 0;
  for (PVector ray: rays){
    ArrayList<Float> rayscalars = raytest.rayscalars;
    ArrayList<Float> rayscalarincs = raytest.rayscalarincs;
    Float rayscalar = rayscalars.get(i);
    Float rayscalarinc = rayscalarincs.get(i);
    rayscalar += rayscalarinc;
    PVector cray = ray.mult(rayscalar);
    Float dist = PointToPlaneDist(disk2, cray);
    if 
    i += 1;
  }
}

ArrayList<Disk> disks = new ArrayList<Disk>();
ArrayList<ArrayList<PVector>> diskspoints = new ArrayList<ArrayList<PVector>>();

float phidiv = 3;
float thetadiv = 3;
float thetadiv2 = 30.0;
float deltaphi = PI/3.0;
float sphereRad = 200.0;
float tolerance = 1.0e-8;
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
boolean DrawSphere = false;
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
  t += .007;
}