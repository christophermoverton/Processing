float G =  6.67428e-11;
float SunMass = 1.98892e30;
float EarthMass = 5.9742e24;
float VenusMass = 4.8675e24;
float AU = 1.496e11;
float SCALE = 250 / AU;
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
Body[] bodies = new Body[3];
PVector[] BForce = new PVector[3];
double total_fx;
double total_fy;

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
  if (increment < 5){
    println("Distance:");
    println(d);
    println("Force: ");
    println(f);
  }
  float ra = PVector.angleBetween(b1.position, b2.position);
  double fx = cos(ra)*(float)f;
  double fy = sin(ra)*(float)f;
  force.x = (float)fx;
  force.y = (float)fy;
}

//void circularorbit(Body b1, Body b2, PVector velocity){
//  //assuming that the mass of b1 is greater relative to all other bodies in 
//  // a given system, which are considered for intial starting point,
//  //negligible in so far as gravitational interactions.
//  //This computes the necessary tangential velocity for in defining a circular 
//  //orbit.  Or in another words force due to centripetal acceleration
//  PVector force = new PVector(0.0,0.0);
//  attraction(b1,b2,force);
//  float fmag = force.mag();
//  float d = PVector.dist(b1.position, b2.position);
//  println("Distance: ");
//  println(d);
//  float v1 = pow(fmag*d/b2.mass,.5);
//  float ra = PVector.angleBetween(b1.position, b2.position);
//  ra += PI/2; // represents angle of velocity vector (relative to force vector)
//  velocity.x = cos(ra)*v1;
//  velocity.y = sin(ra)*v1;
//}

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
    if (increment < 5){
      println("Total force x: ");
      println(total_fx);
    }
    BForce[i] = new PVector((float)total_fx, (float)total_fy);
  }
  for (int i = 0; i < bodies.length; i++) {
    double fx = BForce[i].x;
    double fy = BForce[i].y;
    bodies[i].velocity.x += (float)((double)fx/(double)bodies[i].mass * (double)timestep);
    bodies[i].velocity.y += (float)((double)fy/(double)bodies[i].mass * (double)timestep);
    
    bodies[i].position.x += bodies[i].velocity.x * timestep;
    bodies[i].position.y += bodies[i].velocity.y * timestep;
    if (i == 0){
      rd = 100;
    }
    else{
      rd = 30;
    }
    //fill(bodies[i].pcolor);
    //println(bodies[i].position.x);
    //println(bodies[i].position.y);
    ellipse(bodies[i].position.x * SCALE/AU, bodies[i].position.y * SCALE/AU,rd,rd);
    
  }
  increment += 1;
}