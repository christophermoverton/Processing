ParticleSystem ps;
float MaximumParticles = 100000;
float diffParticleVelocity = .0001;
float minParticleVelocity = .001f;
float minParticleSize = 0.01f;
float maxParticleSize = 8.0f;
float turbulence = 1.0f;
float maxparticleLifetime = 2.0f;
float minparticleLifetime = 100.0f;
float particleStrength = .1f;
float minstrength =.1f;
float tstart = 0.0f;
float tstop = 1.0f;

//control point data P0 (4,0)
//control point data P1 (5,4)
//P2 (3.47,7.5)
//P3 (0,10)
void setup() {
  size(640,360);
  ps = new ParticleSystem(new PVector(width/2,150));

}

void draw() {
  background(0);
  //fill(0);
  //rect(0,0,640,360);
  for (int i =0; i < 4000; i++){
    ps.addParticle();
    if (ps.particles.size()>MaximumParticles){
      break;
    }
  }
  //ps.addParticle();
  ps.run();
  //saveFrame();
}

PVector bezierCurve(PVector P0, PVector P1, PVector P2, PVector P3, float t){
  float ot = (1-t);
  PVector c1 = (P0.copy()).mult(ot*ot*ot);
  PVector c2 = (P1.copy()).mult(3.0f*(ot*ot)*t);
  PVector c3 = (P2.copy()).mult(3.0f*ot*t*t);
  PVector c4 = (P3.copy()).mult(t*t*t);
  return ((c1.add(c2)).add(c3)).add(c4);
}


// A class to describe a group of Particles
// An ArrayList is used to manage the list of Particles 

class ParticleSystem {
  ArrayList<Particle> particles;
  PVector origin;

  ParticleSystem(PVector location) {
    origin = location.copy();
    particles = new ArrayList<Particle>();
  }

  void addParticle() {
    particles.add(new Particle(origin));
  }

  void run() {
    for (int i = particles.size()-1; i >= 0; i--) {
      Particle p = particles.get(i);
      p.run();
      if (p.isDead()) {
        particles.remove(i);
      }
    }
  }
}



// A simple Particle class

class Particle {
  PVector location;
  PVector velocity;
  PVector acceleration;
  float lifespan;
  float redval;
  float greenval;
  float blueval;
  float size;
  float t;
  float tstep;
  float scale;
  PVector P0;
  PVector P1;
  PVector P2;
  PVector P3;
  PVector origin;
  float strength;
  
  void Stept(){
    t+=tstep;
    t=t%1;
  }

  Particle(PVector l) {
    acceleration = new PVector(0,0.05);
    velocity = new PVector(random(-1,1),random(1,2));
    location = l.copy();
    origin = l.copy();
    tstep = random(minParticleVelocity,  minParticleVelocity+diffParticleVelocity);
    //float rval = random(1,20);
    //location.y += rval;
    //generate random parameter t 0 to 1
    t = random(tstart,tstop);
    //control point data P0 (4,0)
    //control point data P1 (5,4)
    //P2 (3.47,7.5)
    //P3 (0,10)
    scale = random(.01, 10.0f);
    float rval1 = random(0,1);
    if (rval1 > .5){
    P0 = new PVector(4.0,0.0,0.0);
    P1 = new PVector(5.0,4.0,0.0);
    P2 = new PVector(3.47,7.5,0.0);
    P3 = new PVector(0.0,10.0,0.0);
    }
    else{
    P0 = new PVector(-4.0,0.0,0.0);
    P1 = new PVector(-5.0,4.0,0.0);
    P2 = new PVector(-3.47,7.5,0.0);
    P3 = new PVector(0.0,10.0,0.0);
    }
    P0.mult(scale); P1.mult(scale); P2.mult(scale);
    P3.mult(scale);
    //P0a.mult(scale); P1a.mult(scale); P2a.mult(scale);
    //P3a.mult(scale);
    PVector c1 = bezierCurve(P0, P1, P2, P3, t);
    //PVector c2 = bezierCurve(P0a, P1a, P2a, P3a, t);
    location.y = c1.y + origin.y;
    location.x = c1.x + origin.x;
    /*
    if (random(0,1) > .5){
      float rval = c1.x;
      location.y += c1.y;
      location.x += random(0,rval);
    }
    else{
      float rval = c2.x;
      location.y += c2.y;
      location.x += random(rval,0);
    }
    */
    lifespan = random(minparticleLifetime,maxparticleLifetime)/velocity.mag();
    PVector loc = location.copy();
    /*
    if ((loc.sub(l.copy())).mag() > 5){
      
      redval = random(0,255);
    
      greenval = random(0,255);
      blueval = random(200,255);
    }
    else{
      redval = random(200,255);
    
      greenval = random(0,10);
      blueval = random(0,30);
    }
    */
    PVector colval = new PVector(scale/10.0f*255.0f,200.0f, (1.0f-scale/10.0f)*255.0f);
    redval = colval.x;
    greenval = colval.y;
    blueval = colval.z;
    size = random(minParticleSize,maxParticleSize);
    strength = random(minstrength,particleStrength);
  }

  void run() {
    updateFlowPosition();
    display();
  }

  // Method to update location
  void update() {
    velocity.add(acceleration);
    location.add(velocity);
    lifespan -= 1.0;
  }
  
  void updateFlowPosition(){
    Stept();
    PVector c1 = bezierCurve(P0, P1, P2, P3, t);
    location.x = origin.x + c1.x;
    location.y = origin.y + c1.y;
    lifespan -= 1.0;
    P3.y += random(-scale*turbulence,scale*turbulence);
    P3.x += random(-scale*.4*turbulence,scale*.4*turbulence);
  }

  // Method to display
  void display() {
    //stroke(255,lifespan);
    noStroke();
    fill(redval, greenval, blueval,lifespan*strength);
    ellipse(location.x,location.y,size*random(.1,1), size*random(.1,1));
  }
  
  // Is the particle still useful?
  boolean isDead() {
    if (lifespan < 0.0) {
      return true;
    } else {
      return false;
    }
  }
}