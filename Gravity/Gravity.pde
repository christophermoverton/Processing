float G =  6.67428e-11;
float SunMass = 1.98892e30;
float EarthMass = 5.9742e24;
float VenusMass = 4.7885e24;
float AU = (149.6e6 * 1000);
float SCALE = 250 / AU;
float Earthpx = -1*AU;
float Earthvy = 29.783*1000;
float Venuspx = 4.8685e24;
float Venusvy = -35.02*1000;
Body[] bodies = new Body[3];
class Body{
  PVector position;
  float mass;
  PVector velocity;
  
  Body(PVector pos, PVector vel, float m){
    position = pos;
    mass = m;
    velocity = vel;
  }
}

void attraction(Body b1, Body b2, PVector force){
  float d = b1.position.dist(b2.position);
  float f = G*b1.mass*b2.mass/pow(d,2);
  float ra = PVector.angleBetween(b1.position, b2.position);
  float fx = cos(ra)*f;
  float fy = sin(ra)*f;
  force.x = fx;
  force.y = fy;
}

void circularorbit(Body b1, Body b2, PVector velocity){
  //assuming that the mass of b1 is greater relative to all other bodies in 
  // a given system, which are considered for intial starting point,
  //negligible in so far as gravitational interactions.
  //This computes the necessary tangential velocity for in defining a circular 
  //orbit.  Or in another words force due to centripetal acceleration
  PVector force = new PVector(0.0,0.0);
  attraction(b1,b2,force);
  float fmag = force.mag();
  float d = b1.position.dist(b2.position);
  float v1 = pow(fmag*d/b2.mass,.5);
  float ra = PVector.angleBetween(b1.position, b2.position);
  ra += PI/2; // represents angle of velocity vector (relative to force vector)
  velocity.x = cos(ra)*v1;
  velocity.y = sin(ra)*v1;
}

void setup() {
  size(1920,1080);
  background(0);
  
}