ParticleSystem ps;
float MaximumParticles = 320;  //120-170 excellent
float diffParticleVelocity = .0001;
float minParticleVelocity = .001f;
float minParticleSize = 0.01f;
float maxParticleSize = 1.0f;
float turbulence = .5f;
float maxparticleLifetime = 100.0f;  //1.0f excellent  10.0f max 30.0f min
float minparticleLifetime = 30.0f;  // 4.0f excellent  //30.0f min 100.0f max slower
float particleStrength = 2.0f;   //2.0f excellent
float minstrength =.1f;   //.1 excellent
float tstart = 0.0f;
float tstop = 1.0f;
int glowRadius = 6;
boolean glowing = true;
boolean randomP3 = false;
float minRandomt3 = .95f;
float glowdampening = 5f;  //1.4 - 2.5 excellent
float glowRaddiv = 375.0f;  //300 - 400.0f dampening on glow radius alpha
//control point data P0 (4,0)
//control point data P1 (5,4)
//P2 (3.47,7.5)
//P3 (0,10)
void setup() {
  size(640,360);
  background(0);
  ps = new ParticleSystem(new PVector(width/2,150));

}

void draw() {
  //background(0,0,0,50);
  fill(0,100);
  rect(0,0,640,360);
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
    if (randomP3){
      //get three random control points greater than P0
      float t3 = random(minRandomt3,1);
      float t2 = random(0,t3);
      float t1 = random(0,t2);
      P1 = bezierCurve(P0,P1,P2,P3,t1);
      P2 = bezierCurve(P0,P1,P2,P3,t2);
      P3 = bezierCurve(P0,P1,P2,P3,t3);
    }
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
    if(glowing){
      //float glowRadius = 100.0;
      //float glowRadius = 100.0 + 15 * sin(frameCount/(3*frameRate)*TWO_PI); 
      //strokeWeight(2);
      //noFill();
      //fill(255,0);

      for(int i = 0; i <= glowRadius; i++){
        noFill();
        strokeWeight(i*glowdampening);
        stroke(200*(1-.1*float(i)/glowRadius),200*(1-.1*float(i)/glowRadius),255,255.0*(1-float(i)/glowRadius)*(float(i)/glowRaddiv));


        bezier(origin.x+P0.x,origin.y+P0.y,origin.x+P1.x,origin.y+P1.y,
           origin.x+P2.x,origin.y+P2.y,origin.x+P3.x,origin.y+P3.y);
   
      }
    strokeWeight(size*.05f);
    noFill();
    //fill(redval, greenval, blueval,lifespan*strength);
    
    stroke(100,100,255,255.0*(1-size/maxParticleSize)*strength);
    bezier(origin.x+P0.x,origin.y+P0.y,origin.x+P1.x,origin.y+P1.y,
           origin.x+P2.x,origin.y+P2.y,origin.x+P3.x,origin.y+P3.y);
    //ellipse(location.x,location.y,size*random(.1,1), size*random(.1,1));

    }
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

class Thing
{
  PVector loc;
  boolean glowing;
  float size;
  Thing(int x, int y, boolean glow){
    loc = new PVector(x,y);
    glowing = glow;
  }
  void display(){
    if(glowing){
      //float glowRadius = 100.0;
      float glowRadius = 100.0 + 15 * sin(frameCount/(3*frameRate)*TWO_PI); 
      strokeWeight(2);
      fill(255,0);
      for(int i = 0; i < glowRadius; i++){
        stroke(255,255.0*(1-size/maxParticleSize));
        ellipse(loc.x,loc.y,i,i);
      }
    }
    //irrelevant
    stroke(0);
    fill(0);
    ellipseMode(CENTER);
    ellipse(loc.x,loc.y,40,40);
    stroke(0,255,0);
    line(loc.x,loc.y+30,loc.x,loc.y-30);
    line(loc.x+30,loc.y,loc.x-30,loc.y);
  }
}