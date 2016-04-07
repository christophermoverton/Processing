import java.util.Map;
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
  ArrayList<PVector> intersectionpoints;
  ArrayList<Boolean> stopTest;
  ArrayList<PVector> currentPoints;
  PVector startPoint;
  rayTest(int num){
    rayscalars = new ArrayList<Float>();
    intersection = new ArrayList<Boolean>();
    rayscalarincs = new ArrayList<Float>();
    lastDistances = new ArrayList<Float>();
    intersectionpoints = new ArrayList<PVector>();
    stopTest = new ArrayList<Boolean>();
    currentPoints = new ArrayList<PVector>();
    startPoint = new PVector(0.0,0.0,0.0);
    for (int i = 0; i < num; i++){
      rayscalars.add(0.0);
      intersection.add(false);
      rayscalarincs.add(.5);
      lastDistances.add(9.0e10);
      intersectionpoints.add(new PVector(0.0,0.0,0.0));
      currentPoints.add(new PVector(0.0,0.0,0.0));
      stopTest.add(false);
      
    }
  }
}

class irayTest{
  Float rayscalar;
  Boolean intersection;
  Float rayscalarinc;
  Float lastDistance;
  PVector intersectionpoint;
  Boolean stopTest;
  PVector currentPoint;
  PVector startPoint;
  irayTest(){
    rayscalar = 0.0;
    intersection = false;
    rayscalarinc = 2.0;
    lastDistance = 1.0e10;
    intersectionpoint = new PVector(0.0,0.0,0.0);
    stopTest = false;
    currentPoint = new PVector(0.0,0.0,0.0);
    startPoint = new PVector(0.0,0.0,0.0);
  }
  irayTest(Float Rayscalar){
    rayscalar = Rayscalar;
    intersection = false;
    rayscalarinc = 1.1;
    lastDistance = 1.0e10;
    intersectionpoint = new PVector(0.0,0.0,0.0);
    stopTest = false;
    currentPoint = new PVector(0.0,0.0,0.0);
    startPoint = new PVector(0.0,0.0,0.0);
  }
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
    float theta = PI/3.0;
    for (int i = 0; i < 6; i++){
      PVector point = Qrotate(p1toCorigin, ncorigin, theta);
      rays.add(point);
      theta += PI/3.0;
    }
    return rays;
}

Float PointToPlaneDist(Disk disk, PVector ray){
    PVector rayc = new PVector(ray.x, ray.y, ray.z);
    PVector corigin = sphericalToC(disk.origin);
    PVector pt1 = sphericalToC(disk.dpoint);
    //PVector p1toCorigin = pt1.sub(corigin);
    PVector ncorigin = new PVector(corigin.x,corigin.y,corigin.z);
    ncorigin.normalize();
    PVector raytopt1 = rayc.sub(pt1);
    Float d = raytopt1.dot(ncorigin);
    //ncorigin.mult(d);
    //ray.add(ncorigin);
    return raytopt1.dot(ncorigin);
}

Float PlaneAndLineIntersection(Disk disk, PVector ray, PVector rayOrigin){
  //returns u parameter (scalar) distance 
  PVector rayc = new PVector(ray.x, ray.y, ray.z);
  PVector corigin = sphericalToC(disk.origin);
  PVector pt1 = sphericalToC(disk.dpoint);
  pt1 = new PVector(pt1.x,pt1.y,pt1.z);
    //PVector p1toCorigin = pt1.sub(corigin);
  PVector ncorigin = new PVector(corigin.x,corigin.y,corigin.z);
  ncorigin.normalize();
  return (ncorigin.dot(pt1.sub(rayOrigin)))/(ncorigin.dot(rayc));
}

void findIntersection(Disk disk, Disk disk2, PVector ray, irayTest iraytest){
  int bcond = 0;
  Boolean bcondc = bcond < 100000; 
  Float d = 1.0e14;
  Float lastDistance = d;
  PVector corigin = sphericalToC(disk.origin);
  int rcount = 0;
  while ((abs(d) > tolerance ) && (bcondc)){
    Float rayscalar = iraytest.rayscalar;
    Float rayscalarinc = iraytest.rayscalarinc;
    //Float rayscalar = rayscalars.get(i);
    //Float rayscalarinc = rayscalarincs.get(i);
    rayscalar += rayscalarinc;
    PVector rayc = new PVector(ray.x, ray.y, ray.z);
    rayc.mult(rayscalar);
    rayc.add(corigin);
    d = PointToPlaneDist(disk2, rayc);
    if (bcond == 0){
      lastDistance = abs(d);
    }
    Boolean t1 = lastDistance < 0;
    Boolean t2 = d < 0;
    //if ((t1 && t2)||(!t1 && !t2)){
    if (abs(d) <= lastDistance){
      lastDistance = abs(d);
      iraytest.rayscalar = rayscalar;
    }
    else{
      if (rayscalarinc > 1.0){
        //println("hit");
        iraytest.rayscalarinc = 1.0/rayscalarinc;
      }
      iraytest.rayscalarinc = pow(rayscalarinc, 2);
      if (rcount > 10){
        iraytest.stopTest = true;
        println("hit");
        break;
      }
      else{
        rcount += 1;
      }
    }
    println(rcount);
    iraytest.currentPoint.x = rayc.x;
    iraytest.currentPoint.y = rayc.y;
    iraytest.currentPoint.z = rayc.z;
    bcond += 1;
    bcondc = bcond < 10;
    if (d < tolerance){
      println("hit d less than tolerance");
    }
  }
  iraytest.intersectionpoint.x = iraytest.currentPoint.x;
  iraytest.intersectionpoint.y = iraytest.currentPoint.y;
  iraytest.intersectionpoint.z = iraytest.currentPoint.z;
}

void findIntersection(Disk disk2, PVector ray, irayTest iraytest){
   float u = PlaneAndLineIntersection(disk2, ray, iraytest.startPoint);
   PVector rayc = new PVector(ray.x,ray.y,ray.z);
   PVector ip = rayc.mult(u);
   ip.add(iraytest.startPoint);
   iraytest.intersectionpoint = new PVector(ip.x,ip.y,ip.z);
}

void stopRayInitialization(Disk disk, Disk disk2, ArrayList<PVector> rays, 
                           rayTest raytest){
  int i = 0;
  PVector corigin = sphericalToC(disk.origin);
  raytest.startPoint.x = corigin.x;
  raytest.startPoint.y = corigin.y;
  raytest.startPoint.z = corigin.z;
  for (PVector ray: rays){
    
      ArrayList<Float> rayscalars = raytest.rayscalars;
      ArrayList<Float> rayscalarincs = raytest.rayscalarincs;
      Float rayscalar = rayscalars.get(i);
      Float rayscalarinc = rayscalarincs.get(i);
      rayscalar += rayscalarinc;
      PVector rayc = new PVector(ray.x, ray.y, ray.z);
      rayc.mult(rayscalar);
      rayc.add(corigin);
      Float dist = PointToPlaneDist(disk2, rayc);
      irayTest iraytest = new irayTest(rayscalar);
      iraytest.startPoint.x = raytest.startPoint.x;
      iraytest.startPoint.y = raytest.startPoint.y;
      iraytest.startPoint.z = raytest.startPoint.z;
      //if (dist < tolerance){
      //  raytest.intersection.set(i, true);
      //  raytest.intersectionpoints.set(i, cray);
      //  break;
      //}
      //else{
        //rayscalar += dist;
        //cray = ray.mult(rayscalar);
        //cray = cray.add(corigin);
        //raytest.intersectionpoints.set(i,rayc);
        if (abs(dist) < 500.0){
          raytest.lastDistances.set(i,abs(dist));
          //findIntersection(disk, disk2, ray, iraytest);
          findIntersection(disk2, ray, iraytest);
          if (iraytest.stopTest){
            raytest.lastDistances.set(i,1.0e7);
            raytest.intersectionpoints.set(i, corigin);
            raytest.stopTest.set(i,true);           
          }
          else{
            raytest.intersectionpoints.set(i, iraytest.intersectionpoint);
          }
        }
        else{
          raytest.lastDistances.set(i,1.0e7);
          raytest.intersectionpoints.set(i, corigin);
          raytest.stopTest.set(i,true);
        }
        int j = 0;
        //ArrayList<Boolean> dtest = new ArrayList<Boolean>();
        float dist2 = dist;
        float rayscalar2 = rayscalar;
        //while (j < 10){
        //  rayc = new PVector(ray.x, ray.y, ray.z);
        //  rayscalarinc = pow(rayscalarinc, .5);
        //  rayscalar = rayscalars.get(i);
        //  rayscalar += rayscalarinc;
        //  rayc.mult(rayscalar);
        //  rayc.add(corigin);
        //  dist = PointToPlaneDist(disk2, rayc);
        //  if (dist > raytest.lastDistances.get(i)){
        //    raytest.stopTest.set(i,true);
        //    raytest.intersectionpoints.set(i, corigin);
        //    break;
        //  }
        //  j+=1;
        //}
        raytest.rayscalarincs.set(i,1.0);
      
      i += 1;
    
  } 
}

void rayTester(Disk disk, Disk disk2, 
               ArrayList<PVector> rays, rayTest raytest){
  int i = 0;
  PVector corigin = sphericalToC(disk.origin);
  for (PVector ray: rays){
    if(!raytest.stopTest.get(i)){
      PVector rayc = new PVector(ray.x, ray.y, ray.z);
      ArrayList<Float> rayscalars = raytest.rayscalars;
      ArrayList<Float> rayscalarincs = raytest.rayscalarincs;
      Float rayscalar = rayscalars.get(i);
      Float rayscalarinc = rayscalarincs.get(i);
      rayscalar += rayscalarinc;
      PVector cray = rayc.mult(rayscalar);
      cray = cray.add(corigin);
      Float dist = PointToPlaneDist(disk2, cray);
      if (abs(dist) > raytest.lastDistances.get(i)){
        raytest.stopTest.set(i,true);
      }
      else{
        raytest.lastDistances.set(i, abs(dist));
        raytest.rayscalars.set(i,rayscalar);
        raytest.currentPoints.set(i, cray);
        
      }
      i += 1;
    }
  }
}

void planePlaneIntersection(){
}

ArrayList<Disk> disks = new ArrayList<Disk>();
ArrayList<ArrayList<PVector>> diskspoints = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> diskRays = new ArrayList<ArrayList<PVector>>();
ArrayList<HashMap<Integer, rayTest>> diskRayTestMaps = new ArrayList<HashMap<Integer, rayTest>>();
//ArrayList<rayTest> diskRayTest = new ArrayList<rayTest>();

float phidiv = 3;
float thetadiv = 3;
float thetadiv2 = 30.0;
float deltaphi = PI/3.0;
float sphereRad = 200.0;
float tolerance = 1.0e-8;
float sInvert = 1.8;
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
    ArrayList<PVector> rays = generateRays(disk);
    diskRays.add(rays);
    int j = 0;
    HashMap<Integer, rayTest> diskraytestmap = new HashMap<Integer, rayTest>();
    for (Disk disk2: disks){
      if (disk == disk2){
        
        continue;
      }
      rayTest raytests = new rayTest(rays.size());
      stopRayInitialization(disk, disk2, rays, raytests);
      println(raytests.stopTest);
      diskraytestmap.put(j,raytests);
      j += 1;
    }
    diskRayTestMaps.add(diskraytestmap);
  }
  
  size(1920,1080, P3D);
  background(0);
  lights();
}

float t = 0.02;
float tinc = .007;
boolean DrawSphere = false;
boolean DrawDisks = true;
boolean Diskfill = false;
void draw(){
  fill(0,80);
  noStroke();
  beginShape();
  vertex(-5000,-3000,-200);
  vertex(5000,-3000, -200);
  vertex(5000,3000,-200);
  vertex(-5000,3000,-200);
  endShape();
  //rect(0,0, 1920,1080);
  translate(1920/2, 1080/2, 0);
  rotateY(t);
  //scale(t);
  noFill();
  stroke(200);
  strokeWeight(.3);
  sphereDetail(20);
  if (DrawSphere){
    sphere(200);
  }
  int fillc = 0;
  if (DrawDisks){
    for (ArrayList<PVector> diskpoints: diskspoints){
      if (Diskfill){
        fill(fillc%255,2);
      }
      else{
        noFill();
      }
      beginShape();
      
      for (PVector point: diskpoints){
        vertex(point.x,point.y,point.z);
        //println(point);
      }
      endShape();
      fillc += 10;
    }
  }
  int i = 0;
  stroke(255);
  strokeWeight(1.0);
  for(HashMap<Integer, rayTest> diskRayTestMap: diskRayTestMaps){
    for(Map.Entry<Integer, rayTest> me: diskRayTestMap.entrySet()){
      Integer j = me.getKey();
      rayTest raytest = me.getValue();
      rayTester(disks.get(i), disks.get(j), diskRays.get(i), raytest);
      int k = 0;
      for (Boolean stopt : raytest.stopTest){
        PVector sp = raytest.startPoint;
        if (stopt){
          
          PVector ep = raytest.intersectionpoints.get(k);
          if (ep.x != sp.x && ep.y != sp.y){
            line(sp.x,sp.y,sp.z, ep.x,ep.y,ep.z);
              textSize(5);
              text("Ray "+ String.valueOf(i)+ " "+ j.toString()+ " "+ String.valueOf(k),
                  ep.x,ep.y,ep.z);
                  fill(255);
          }
        }
        else{
          PVector ep = raytest.currentPoints.get(k);
          line(sp.x,sp.y,sp.z, ep.x,ep.y,ep.z);
        }
        k += 1;
      }
    }
    i+= 1;
  }
  t += tinc;
  //println(t);
  if (t > 3){
    tinc = -tinc;
  }
  if (t < .01){
    tinc = -tinc;
  }
}