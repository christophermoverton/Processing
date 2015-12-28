float G =  6.67428e-11;
float SunMass = 1.98892e30;
float EarthMass = 5.9742e24;
float VenusMass = 4.8675e24;
float AU = 1.496e11;
float SCALE = 100 / AU;
float Earthpx = -1*AU;
float Earthvy = 29.783*1000;
float Venuspx = .723 * AU;
float Venusvy = -35.02*1000;
float timestep = 24*3600;
color EarthColor = color(200,200,255);
color VenusColor = color(255, 200,200);
color SunColor = color(255,255,200);
int increment = 0;
float rd;
Body[] bodies = new Body[20];
PVector[] BForce = new PVector[20];
double total_fx;
double total_fy;
float minPixeldim = 1080;

class Body{
  PVector position;
  double mass;
  PVector velocity;
  color pcolor;
  
  Body(PVector pos, PVector vel, double m, color pc){
    position = pos;
    mass = m;
    velocity = vel;
    pcolor = pc;
  }
}

void attraction(Body b1, Body b2, PVector force){
  float d = PVector.dist(b1.position, b2.position);
  double f = (double)G*(double)b1.mass*(double)b2.mass/((double)d*(double)d);
  //if (increment < 5){
  //  println("Distance:");
  //  println(d);
  //  println("Force: ");
  //  println(f);
  //}
  float ra = PVector.angleBetween(b1.position, b2.position);
  float dx = (b2.position.x - b1.position.x);
  float dy = (b2.position.y - b1.position.y);
  //float rdydx = (float)((double)dy/ (double) dx);
  ra = atan2(dy,dx);
  double fx = cos(ra)*(float)f;
  double fy = sin(ra)*(float)f;
  force.x = (float)fx;
  force.y = (float)fy;
  force.z = (float)f;
}

void getCircularOrbit(Body b1, Body b2, PVector velocity){
 //assuming that the mass of b1 is greater relative to all other bodies in 
 // a given system, which are considered for intial starting point,
 //negligible in so far as gravitational interactions.
 //This computes the necessary tangential velocity for in defining a circular 
 //orbit.  Or in another words force due to centripetal acceleration
 PVector force = new PVector(0.0,0.0);
 attraction(b1,b2,force);
 //float fmag = force.mag();
 float fmag = force.z;
 println("fmag: ");
 println(fmag);
 double d = PVector.dist(b1.position, b2.position);
 println("Distance: ");
 println(d);
 float v1 = pow((float)((double)fmag*d/ (double)b2.mass),.5);
 float dx = (b2.position.x - b1.position.x);
 float dy = (b2.position.y - b1.position.y);
  //float rdydx = (float)((double)dy/ (double) dx);
 float ra = atan2(dy,dx);
 //float ra = PVector.angleBetween(b1.position, b2.position);
 ra += PI/2; // represents angle of velocity vector (relative to force vector)
 velocity.x = v1;
 velocity.y = ra;
}

void getEscapeVelocity(Body b1, Body b2, PVector velocity){
  float d = PVector.dist(b1.position, b2.position);
  float ve = pow((float) ( (double)2.0 * (double) G* (double) b1.mass/ (double) d ), .5);
 float dx = (b2.position.x - b1.position.x);
 float dy = (b2.position.y - b1.position.y);
  //float rdydx = (float)((double)dy/ (double) dx);
 float ra = atan2(dy,dx); 
 ra += PI/2; // represents angle of velocity vector (relative to force vector)
 velocity.x = ve;
 velocity.y = ra;
}

void genRandPosition(Body b2){
  //generates initially in pixel dimension then converts to m using scale conversion
  float distance = random(100 , .7 * minPixeldim/2.0);
  distance = (float)((double) distance / (double) SCALE); 
  float angle = random(0, 2*PI);
  b2.position.x = (float)((double)distance* (double)cos(angle));
  b2.position.y = (float)((double)distance* (double)sin(angle));
}

void genRandVelocity(Body b1, Body b2){
  //This stores the computed random velocity in b2.
  //Randomization boundaries are set as follows orbits are assigned minimized 
  // so as neither to exceed escape velocity as initialization and/or equal to
  //orbit set at eccentricity zero (circular orbit).
  
  //b1 is assumed to be a star mass in such system.
  PVector cvelocity = new PVector(0.0,0.0);
  PVector evelocity = new PVector(0.0,0.0);
  getCircularOrbit(b1,b2, cvelocity);
  getEscapeVelocity(b1,b2,evelocity);
  float minvel = cvelocity.x;
  float maxvel = evelocity.x;
  println("Minimum velocity: ");
  println(minvel);
  println(maxvel);
  float rvelocitymag = random(minvel,.8*maxvel);
  b2.velocity.x = rvelocitymag*cos(cvelocity.y);
  b2.velocity.y = rvelocitymag*sin(cvelocity.y);
}

void genRandMass(Body b2){
  float rfloat = random((float)1e9,(float)1.9e27);
  b2.mass = rfloat;
}

void setup() {
  size(1920,1080);
  background(0);
  Body Earth = new Body(new PVector(Earthpx,0.0),new PVector(0.0,Earthvy),
  EarthMass, EarthColor);
  Body Sun = new Body(new PVector(0.0,0.0), new PVector(0.0,0.0), 
  SunMass, SunColor);
  Body Venus = new Body(new PVector(Venuspx,0.0), new PVector(0.0, Venusvy), 
  VenusMass, VenusColor);
  bodies[0] = Sun;
  bodies[1] = Venus;
  bodies[2] = Earth;
  println(bodies[1].position.x);
  println(bodies[2].position.x);
  println(PVector.dist(bodies[1].position, bodies[2].position));
  for (int i = 3; i < 20; i++){
    Body binit = new Body(new PVector(0.0,0.0), new PVector(0.0,0.0),0.0, EarthColor);
    genRandMass(binit);
    genRandPosition(binit);
    genRandVelocity(bodies[0],binit);
    println(binit.velocity);
  }
  
}

void draw() {
  fill(0,2);
  noStroke();
  rect(0,0, 1920,1080);
  fill(255);
  translate(1920.0/2.0, 1080.0/2.0);
  for (int i = 0; i < bodies.length; i++){
    total_fx = 0.0;
    total_fy = 0.0;
    Body b1 = bodies[i];
    for (int j = 0; j < bodies.length; j++){
      Body b2 = bodies[j];
      if (i == j){
        continue;
      }
      PVector force = new PVector(0.0,0.0);
      attraction(b1,b2,force);
      total_fx += force.x;
      total_fy += force.y;
    }
    //if (increment < 5){
    //  println("Total force x: ");
    //  println(total_fx);
    //}
    BForce[i] = new PVector((float)total_fx, (float)total_fy);
  }
  for (int i = 0; i < bodies.length; i++) {
    double fx = BForce[i].x;
    double fy = BForce[i].y;
    bodies[i].velocity.x += (float)((double)fx/ (double)bodies[i].mass * (double)timestep);
    bodies[i].velocity.y += (float)((double)fy/ (double)bodies[i].mass * (double)timestep);
    
    bodies[i].position.x += (float)((double)bodies[i].velocity.x * (double) timestep);
    bodies[i].position.y += (float)((double)bodies[i].velocity.y * (double)timestep);
    if (i == 0){
      rd = 50;
    }
    else{
      rd = 1;
    }
    //fill(bodies[i].pcolor);
    //println(bodies[i].position.x);
    //println(bodies[i].position.y);
    float bpx = (float)((double)bodies[i].position.x * (double)SCALE);
    float bpy = (float)((double)bodies[i].position.y * (double)SCALE);
    //println(bpx);
    //println(bpy);
    ellipse(bpx, bpy,rd,rd);
    
  }
  increment += 1;
}