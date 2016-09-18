ParticleSystem ps;
float MaximumParticles = 300;  //120-170 excellent
float diffParticleVelocity = .01;
float minParticleVelocity = .07f;
float minParticleSize = 0.1f;
float maxParticleSize = 20.0f;
float minScale = 0.01f;
float maxScale = 15.0f;
float turbulence = 0.0f;
float maxparticleLifetime = 100.0f;  //1.0f excellent  10.0f max 30.0f min
float minparticleLifetime = 1.0f;  // 4.0f excellent  //30.0f min 100.0f max slower
float particleStrength = 2.0f;   //2.0f excellent
float minstrength =.1f;   //.1 excellent
float tstart = .8f;
float tstop = 1.0f;
int glowRadius = 8;
boolean glowing = true;
boolean randomP3 = false;
float minRandomt3 = .95f;
float glowdampening = 2.4f;  //1.4 - 2.5 excellent
float glowRaddiv = 400.0f;  //300 - 400.0f dampening on glow radius alpha
boolean oturb = true;
PVector oturbvec = new PVector(0.0,0.0,0.0);
PVector oturbvec2 = new PVector(0.0,0.0,0.0);
PVector oturbvec3 = new PVector(0.0,0.0,0.0);
float oturbfac_x = 0.0f;
float oturbfac_y = 0.0f;
float oturbfac2_x = 0.0f;
float oturbfac2_y = 1.0f;
float oturbfac3_x = 0.2;
float oturbfac3_y = 2.0f;
boolean oturb_random = false;
boolean cyclic_dampening = true;
float oturb_frequency = .2f;
boolean flowField = true;
boolean flowFieldrandomStart = false;

//control point data P0 (4,0)
//control point data P1 (5,4)
//P2 (3.47,7.5)
//P3 (0,10)
void setup() {
  size(640,360);
  background(0);
  ps = new ParticleSystem(new PVector(width/2,150));

}
int t = 0;
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
  float cdamp = 1.0f;
  if (cyclic_dampening){
    cdamp = exp(2.0/(t%100+1));
  }
  if (oturb_random){
    oturbvec.x = random(-cdamp*oturbfac_x, cdamp*oturbfac_x);
    oturbvec.y = random(-cdamp*oturbfac_y, cdamp*oturbfac_y);
    oturbvec2.x = random(-cdamp*oturbfac2_x, cdamp*oturbfac2_x);
    oturbvec2.y = random(-cdamp*oturbfac2_y, cdamp*oturbfac2_y);
    oturbvec3.x = random(-cdamp*oturbfac3_x, cdamp*oturbfac3_x);
    oturbvec3.y = random(-cdamp*oturbfac3_y, cdamp*oturbfac3_y);
  }
  else{
    oturbvec.x = cdamp*oturbfac_x*sin(oturb_frequency*t);
    oturbvec.y = cdamp*oturbfac_y*cos(oturb_frequency*t);
    oturbvec2.x = cdamp*oturbfac2_x*sin(oturb_frequency*t);
    oturbvec2.y = cdamp*oturbfac2_y*cos(oturb_frequency*t);
    oturbvec3.x = cdamp*oturbfac3_x*sin(oturb_frequency*t);
    oturbvec3.y = cdamp*oturbfac3_y*cos(oturb_frequency*t);
  }
  ps.run();
  filter(BLUR, .5);
  t+=1;
  //filter(DILATE);
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
  PVector P4;
  PVector P5;
  PVector P6;
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
    if (flowFieldrandomStart){
      
      t = random(tstart,tstop);
    }
    else{
      t = tstart;
    }
    //control point data P0 (4,0)
    //control point data P1 (5,4)
    //P2 (3.47,7.5)
    //P3 (0,10)
    scale = random(minScale, maxScale);
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
    if (oturb){
      if (!((P3.y + oturbvec3.y) < P2.y)){
         P3.y += oturbvec3.y;
      }
      //if (!((P3.x + oturbvec3.x) < P2.x)){
         P3.x += oturbvec3.x;
      //}
      if (!((P2.y + oturbvec2.y) < P1.y)){
         P2.y += oturbvec2.y;
      }
      //if (!((P2.x + oturbvec2.x) < P1.x)){
         P2.x += oturbvec2.x;
      //}
      P1.y += oturbvec.y;
      P1.x += oturbvec.x;
    }
    if (flowField){
      P6 = bezierCurve(P0,P1,P2,P3,t).copy();
      P5 = bezierCurve(P0,P1,P2,P3,2.0f/3.0f*t);
      P4 = bezierCurve(P0,P1,P2,P3,1.0f/3.0f*t);
    }
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

        if (flowField){
          bezier(origin.x+P0.x,origin.y+P0.y,origin.x+P4.x,origin.y+P4.y,
                 origin.x+P5.x,origin.y+P5.y,origin.x+P6.x,origin.y+P6.y);            
        }
        else{
          bezier(origin.x+P0.x,origin.y+P0.y,origin.x+P1.x,origin.y+P1.y,
             origin.x+P2.x,origin.y+P2.y,origin.x+P3.x,origin.y+P3.y);
        }
      }
    }
    strokeWeight(size*.05f);
    noFill();
    //fill(redval, greenval, blueval,lifespan*strength);
    
    stroke(100,100,255,255.0*(1-size/maxParticleSize)*strength);
    if (flowField){
        bezier(origin.x+P0.x,origin.y+P0.y,origin.x+P4.x,origin.y+P4.y,
               origin.x+P5.x,origin.y+P5.y,origin.x+P6.x,origin.y+P6.y);      
    }
    else{
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