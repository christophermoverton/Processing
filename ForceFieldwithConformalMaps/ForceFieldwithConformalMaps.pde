/*
Conformal mapping of the complex plane z to 1/z 
*/
PVector cmult(PVector a, PVector b){
  return new PVector(a.x*b.x - (a.y*b.y), a.x*b.y + b.x*a.y, 0.0); 
}

PVector cdiv(PVector NumA, PVector DenomA){
  //get conjugate of DenomA
  PVector cNumA = new PVector(NumA.x, NumA.y, 0.0);
  PVector cDenomA = new PVector(DenomA.x,DenomA.y, 0.0);
  PVector ConjDenomA = new PVector(DenomA.x, -DenomA.y, 0.0);
  PVector pdenom = cmult(cDenomA,ConjDenomA);
  PVector pnum = cmult(cNumA,ConjDenomA);
  if (pdenom.x != 0.0){
    return PVector.mult(pnum, 1.0/pdenom.x);
  }
  else{
    return new PVector(0.0,0.0,0.0);
  }
}

ArrayList<PVector> nthroot(PVector c, Float nroot){
  //complex nth root
  ArrayList<PVector> roots = new ArrayList<PVector>();
  Float r = c.mag();
  Float theta = c.heading();
  Float rroot = pow(r, 1.0/nroot);
  for (int i = 0; i < nroot; i++){
    PVector cn = new PVector(cos((theta + i*2*PI)/nroot), sin((theta+i*2*PI)/nroot), 0.0);
    roots.add(PVector.mult(cn, rroot));
  }
  return roots;
}

PVector creq(PVector c, Float n){
  //Cuachy Riemannn equations
  //ArrayList<PVector> roots = new ArrayList<PVector>();
  Float r = c.mag();
  Float theta = c.heading();
  Float rroot = pow(r, n);
  //for (int i = 0; i < nroot; i++){
  //  PVector cn = new PVector(cos((theta + i*2*PI)/nroot), sin((theta+i*2*PI)/nroot), 0.0);
  //  roots.add(PVector.mult(cn, rroot));
  //}
  PVector cn = new PVector(cos(theta*n), sin(theta*n),0.0);
  cn = PVector.mult(cn, rroot);
  return cn;
  
}



PVector ccos(PVector z){
  //complex cos function
  // cos(z) = cos(a)cosh(b) - isin(a) sinh(b)  making use 
  // of cosh x = (e^x + e^-x)/2 for instance.
  //Float a = cos(z.x)*(exp(z.y)+exp(-z.y))/2.0;
  float a = cos(z.x)*(float)Math.cosh(z.y);
  //Float b = 1.0*(sin(z.x)*(exp(-z.y)-exp(z.y))/2.0);
  float b = -1.0*(sin(z.x)*(float)Math.sinh(z.y));
  return new PVector(a,b,0.0);
}

PVector csin(PVector z){
  //complex cos function
  // cos(z) = cos(a)cosh(b) - isin(a) sinh(b)  making use 
  // of cosh x = (e^x + e^-x)/2 for instance.
  Float a = sin(z.x)*(exp(z.y)+exp(-z.y))/2.0;
  Float b = (cos(z.x)*(exp(-z.y)-exp(z.y))/-2.0);
  return new PVector(a,b, 0.0);
}

color EarthColor = color(200,200,255);
Float cplaneinc = 9.0;
Float cplanelineinc = 30.0;
Integer dimX = 1080;
Integer dimY = 720;
Float time = 1.0;
Float timeinc = .5;
boolean popm = true;
Float ForceC = .1;  //force field constant varies the strength of the field
Float minM = .01; // min mass (for random distribution minimum)
Float maxM = 2.0; // max mass
Float minvel = .01;
Float maxvel = 2.0;
Body[] bodies = new Body[20];
PVector[] BForce = new PVector[20];
double total_fx;
double total_fy;


class Body{
  PVector position;
  Float mass;
  PVector velocity;
  PVector force;
  color pcolor;
  
  Body(PVector pos, PVector vel, Float m, color pc){
    position = pos;
    mass = m;
    velocity = vel;
    pcolor = pc;
    force = new PVector(0.0,0.0,0.0);
  }
}

void genRandVelocity(Body b2){
  //This stores the computed random velocity in b2.
  //Randomization boundaries are set as follows orbits are assigned minimized 
  // so as neither to exceed escape velocity as initialization and/or equal to
  //orbit set at eccentricity zero (circular orbit).
  
  //b1 is assumed to be a star mass in such system.
  PVector cvelocity = new PVector(0.0,0.0);

  float rvelocitymag = random(minvel,.8*maxvel);
  b2.velocity.x = rvelocitymag*cos(cvelocity.y);
  b2.velocity.y = rvelocitymag*sin(cvelocity.y);
}

void genRandMass(Body b2){
  float rfloat = random(minM,maxM);
  b2.mass = rfloat;
}

void genRandPosition(Body b2){
  //generates initially in pixel dimension then converts to m using scale conversion
  float distance = random(100 , .7 * dimX/2.0);
  float angle = random(0, 2*PI);
  b2.position.x = (distance*cos(angle));
  b2.position.y = (distance*sin(angle));
}

void attraction(Body b1, PVector force, Integer Ftype){
  //float d = PVector.dist(b1.position, b2.position);
  float f = ForceC*b1.mass;
  PVector fvec = new PVector(0.0,0.0,0.0);
  if (Ftype == 0){
    // ftype 0 is 1/z
    fvec = cdiv(new PVector(1.0, 0.0, 0.0), b1.position);
  }
  else if (Ftype == 1){
    fvec = creq(b1.position, .5);
  }
  else {
     //-(2 i sin(2 theta))/r^3-(2 cos(2 theta))/r^3-(2 sin(2 theta) ( dtheta)/( dr))/r^2
     //+(2 i cos(2 theta) ( dtheta)/( dr))/r^2
     
    PVector fvec1 = creq(PVector.mult(b1.position,1.0e-6), -2.0);
    PVector fvec2 = new PVector(b1.position.x*1.0e-6 + .01*1.0e-6, b1.position.y*1.0e-6, 0.0);
    fvec2 = creq(fvec2,-2.0);
    PVector fvec3 = new PVector(b1.position.x*1.0e-6, b1.position.y*1.0e-6+.01*1.0e-6, 0.0);
    fvec3 = creq(fvec3,-2.0);
    
    PVector fv1v2 = PVector.sub(fvec1,fvec2);

    PVector fv1v3 = PVector.sub(fvec1,fvec3);
    Float r = b1.position.heading();
    Float th = b1.position.mag();
    Float a = -(2*cos(2*th))/pow(r,3.0)-(2.0*sin(2.0*th))/pow(r,2.0);
    Float b = -2*sin(2*th)/pow(r,3.0)+(2*cos(2.0*th))/pow(r,2.0);
    //if (fvec1.heading() > 0 && fvec1.heading() < PI/2.0){
    //  fv1v2 = PVector.mult(fv1v2, -1.0);
    //  fv1v3 = PVector.mult(fv1v3, 1.0);
    //}
    //else if (fvec1.heading() > 0 && fvec1.heading() > PI/2.0){
    //}
    //else if (fvec1.heading() < 0 && fvec1.heading() > -PI/2.0){
      
    //}
    //else{
    //  fv1v2 = PVector.mult(fv1v2, -1.0);
    //}
    fvec = PVector.add(fv1v2,fv1v3);
    fvec.normalize();
    fvec.x = a;
    fvec.y = b;
    fvec.normalize();
  }
  //if (increment < 5){
  //  println("Distance:");
  //  println(d);
  //  println("Force: ");
  //  println(f);
  //}
  //float ra = PVector.angleBetween(b1.position, b2.position);
  //float dx = (b2.position.x - b1.position.x);
  //float dy = (b2.position.y - b1.position.y);
  //float rdydx = (float)((double)dy/ (double) dx);
  //ra = atan2(dy,dx);
  ///double fx = cos(ra)*(float)f;
  //double fy = sin(ra)*(float)f;
  force.x = fvec.x*f;
  force.y = fvec.y*f;
  force.z = f;
}

ArrayList<ArrayList<PVector>> ConformalMap = new ArrayList<ArrayList<PVector>>(); 
ArrayList<ArrayList<PVector>> ConformalMap2 = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> ConformalMap3 = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> ConformalMap4 = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> ConformalMap5 = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> ConformalMap6 = new ArrayList<ArrayList<PVector>>();
void setup(){
  size(1080,720);
  
  PVector a = new PVector(1.0, 2.5, 0.0);
  PVector b = new PVector(1.0, -2.0, 0.0);
  println(cdiv(a,b));
  PVector z = new PVector(-dimX,-dimY,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    ArrayList<PVector> nlist3 = new ArrayList<PVector>();
    
    while(z.y < dimY){
      z.y += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist2.add(creq(z,-2.0));
      nlist3.add(creq(z,2.0));
      
    }
    z.x += cplanelineinc;
    z.y = -dimY;
    ConformalMap.add(nlist);
    ConformalMap5.add(nlist2);
    ConformalMap6.add(nlist3);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    ArrayList<PVector> nlist3 = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist2.add(creq(z,-2.0));
      nlist3.add(creq(z,2.0));
    }
    z.y += cplanelineinc;
    z.x = -dimX;
    ConformalMap.add(nlist);
    ConformalMap5.add(nlist2);
    ConformalMap6.add(nlist3);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.y < dimY){
      z.y += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      ArrayList<PVector> roots = nthroot(z, 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));
      //nlist.add(nroot(z, 2.0));
    }
    z.x += cplanelineinc;
    z.y = -dimY;
    ConformalMap2.add(nlist);
    ConformalMap2.add(nlist2);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      ArrayList<PVector> roots = nthroot(z, 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));     
      //nlist.add(nroot(z, 2.0));
    }
    z.y += cplanelineinc;
    z.x = -dimX;
    ConformalMap2.add(nlist);
    ConformalMap2.add(nlist2);
  }
  //println(ConformalMap);
  z = new PVector(0.0,0.0,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    while(z.y < dimY){
      z.y += cplaneinc;
      PVector sz = PVector.mult(z,1.0/50.0);
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist.add(ccos(sz));
    }
    z.x += cplanelineinc;
    z.y = 0.0;
    ConformalMap3.add(nlist);
  }
  z = new PVector(0.0,0.0,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      PVector sz = PVector.mult(z,1.0/50.0);
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist.add(ccos(sz));
    }
    z.y += cplanelineinc;
    z.x = 0.0;
    ConformalMap3.add(nlist);
    //println(ConformalMap3);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.y < dimY){
      z.y += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      ArrayList<PVector> roots = nthroot(cdiv(new PVector(1.0, 0.0, 0.0), z), 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));
    }
    z.x += cplanelineinc;
    z.y = -dimY;
    ConformalMap4.add(nlist);
    ConformalMap4.add(nlist2);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      ArrayList<PVector> roots = nthroot(cdiv(new PVector(1.0, 0.0, 0.0), z), 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));
    }
    z.y += cplanelineinc;
    z.x = -dimX;
    ConformalMap4.add(nlist);
    ConformalMap4.add(nlist2);
  }
  for (int i = 0; i < bodies.length; i++){
    Body binit = new Body(new PVector(0.0,0.0), new PVector(0.0,0.0),0.0, EarthColor);
    genRandMass(binit);
    genRandPosition(binit);
    //genRandVelocity(binit);
    bodies[i] = binit;
    println(binit.mass);
    println(binit.position);
    println(binit.velocity);
  }
}

void draw(){
  background(0);
  fill(0,2);
  noStroke();
  rect(0,0, 1920,1080);
  translate(dimX/2, dimY/2);
  
  noFill();
  stroke(255,20);
  strokeWeight(.5*1/time);
  //println(time);
  //strokeWeight(1.5*1/time);
  //if (time < 2){
  // scale(time);
  // for (ArrayList<PVector> curvepts : ConformalMap){
  //   beginShape();
  //   for(PVector curvept : curvepts){
  //     curveVertex(10000*curvept.x, 10000*curvept.y);
  //   }
  //   endShape();
  // }
  // fill(255);
  // textSize(12);
  // text("1/z",-100.0,-100.0);  
  //}
  
  //else if (time < 3.0){
  // scale(time);
  // for (ArrayList<PVector> curvepts : ConformalMap2){
  //   beginShape();
  //   for(PVector curvept : curvepts){
  //     curveVertex(5.0*curvept.x, 5.0*curvept.y);
  //  }
  //  endShape();
  //  }
  // fill(255);
  // textSize(6);
  // text("square root z",-100.0,-100.0);  
  //}
  //else if (time < 6.0){
  // strokeWeight(.1*1/time);
  // scale(6.0*time);
  // for (ArrayList<PVector> curvepts : ConformalMap3){
  //   beginShape();
  //   for(PVector curvept : curvepts){
  //     curveVertex(3.0*curvept.x, 3.0*curvept.y);
  //   }
  //   endShape();
  // }
  // fill(255);
  // textSize(1);
  // text("Cos(z)",-5.0,-5.0);  
  //}
  //else if (time < 9.0){
  
  // //if (popm){
  // //  pushMatrix();
  // //  popMatrix();
  // //  translate(dimX/2, dimY/2);
  // //  popm = false;
  // //}
  // scale(time);
  // for (ArrayList<PVector> curvepts : ConformalMap4){
  //   beginShape();
  //   for(PVector curvept : curvepts){
  //     curveVertex(400.0*curvept.x, 400.0*curvept.y);
  //   }
  //   endShape();
  // }
  // fill(255);
  // textSize(5);
  // text("1/sqroot(z)",-30.0,-30.0); 
  //}
  //else if (time < 12.0) {
  // scale(time);
  // for (ArrayList<PVector> curvepts : ConformalMap5){
  //   beginShape();
  //   for(PVector curvept : curvepts){
  //     curveVertex(1000000.0*curvept.x, 1000000.0*curvept.y);
  //   }
  //   endShape();
  // }
  // fill(255);
  // textSize(5);
  // text("1/z^2",-30.0,-30.0); 
  //}
  //else if (time < 15.0){
  // scale(time);
  // for (ArrayList<PVector> curvepts : ConformalMap6){
  //   beginShape();
  //   for(PVector curvept : curvepts){
  //     curveVertex(.0001*curvept.x, .0001*curvept.y);
  //   }
  //   endShape();
  // }
  // fill(255);
  // textSize(5);
  // text("z^2",-10.0,-10.0); 
  //}
  //else{
  //  time = 0.0;
  //}
  //scale(6.0*time);
  //for (ArrayList<PVector> curvepts : ConformalMap3){
  //  beginShape();
  //  for(PVector curvept : curvepts){
  //    curveVertex(curvept.x, curvept.y);
  //  } 
  //  endShape();
  //}
  for (ArrayList<PVector> curvepts : ConformalMap5){
    beginShape();
    for(PVector curvept : curvepts){
      curveVertex(10000000.0*curvept.x, 10000000.0*curvept.y);
    }
    endShape();
  }
   
  fill(255);
  textSize(5);
  text("1/z^2",-30.0,-30.0);
  for (int i = 0; i < bodies.length; i++){
    //total_fx = 0.0;
    //total_fy = 0.0;
    Body b1 = bodies[i];
    //PVector force = new PVector(0.0,0.0,0.0);
    attraction(b1, b1.force, 3);
    //for (int j = 0; j < bodies.length; j++){
    //  Body b2 = bodies[j];
    //  if (i == j){
    //    continue;
    //  }
    //  PVector force = new PVector(0.0,0.0);
    //  attraction(b1,b2,force);
    //  total_fx += force.x;
    //  total_fy += force.y;
    //}
    //if (increment < 5){
    //  println("Total force x: ");
    //  println(total_fx);
    //}
    //BForce[i] = new PVector((float)total_fx, (float)total_fy);
  }
  for (int i = 0; i < bodies.length; i++) {
    //double fx = BForce[i].x;
    //double fy = BForce[i].y;
    bodies[i].velocity.x += bodies[i].force.x/ bodies[i].mass * timeinc;
    bodies[i].velocity.y += bodies[i].force.y/ bodies[i].mass * timeinc;
    
    bodies[i].position.x += bodies[i].velocity.x * timeinc;
    bodies[i].position.y += bodies[i].velocity.y * timeinc;
    float rd;
    if (i == 0){
      rd = 4;
    }
    else{
      rd = 4;
    }
    //fill(bodies[i].pcolor);
    //println(bodies[i].position.x);
    //println(bodies[i].position.y);
    float bpx = bodies[i].position.x ;
    float bpy = bodies[i].position.y ;
    //println(bpx);
    //println(bpy);
    ellipse(bpx, bpy,rd,rd);
    
  }
  time += timeinc;
  //saveFrame();
}