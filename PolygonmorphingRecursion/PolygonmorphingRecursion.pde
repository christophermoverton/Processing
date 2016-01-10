float bfillA = 1.0;
float nbfillA = 1.0;
float angle = -.5;
int NGONs1 = 3;
int NGONs2 = 7;
float RNG1 = 50;
float RNG2 = 150;
PVector[] NG1pos = new PVector[NGONs1];
PVector[] NG2pos = new PVector[NGONs2];
PVector[] ng1posf;
PVector[] ng2posf;
int ng12lcm;
PShape s1;
PShape s2;
float speed = 1.0; 
float time = 0.0;
float atime = 2.0; //(animation time in seconds)
float rtime = 2.0; // (rest time in seconds)
boolean t3 = false; //controls animation cycle
boolean shape = true;  //controls destination animation shape
// recursive implementation of gcd
float maxAngle = 120.0;
int TopLevel = 20;
int gcd(int p, int q) {
  if (q == 0) {return p;}
  else {return gcd(q, p % q);}
}

int lcm(int p, int q) {
  return int(abs(p*q)/gcd(p,q));
}

void Polarcoord(float angle, float radius, PVector xy){
  xy.x = radius*cos(angle);
  xy.y = radius*sin(angle);
 
}

void getNGonPoints(int nside, float radius, PVector[] pos){
  float ainc = 2.0*PI/float(nside);
  float a = -PI/2.0;
  if (nside % 2 == 0){
    a += ainc/2.0;
  }
  for (int i = 0; i < nside; i++){
    PVector point = new PVector(0.0,0.0,0.0);
    Polarcoord(a, radius, point);
    pos[i] = point;
    a += ainc;
  }
}

float linterp(float x, float p1x, float p1y, float p2x, float p2y){
  return p1y+(p2y-p1y)*(x-p1x)/(p2x-p1x);
}
//sdf is subdivision factor
void getNGonSubdivisionPoints(int sdf, PVector[] pos, PVector[] out){
  //println("Pos length: ", pos.length);
  for (int i=0; i< pos.length; i++){
    int ni = i+1;
    if (i == pos.length -1){
      ni = 0;
    }
    PVector pos1 = pos[i];
    PVector pos2 = pos[ni];
    ////pick smallest x and smallest y
    //float sx = min(pos1.x,pos2.x);
    //float sy = min(pos1.y,pos2.y);
    ////get distance between positions
    //float dx = abs(pos1.x-pos2.x);
    //float dy = abs(pos1.y-pos2.y);
    ////get subdivision increment
    //float sbdx = dx/float(sdf);
    //float sbdy = dy/float(sdf);
    //using vector approach in subdividing points.
    PVector p12 = PVector.sub(pos2,pos1);
    float mp12 = p12.mag();
    p12.normalize();
    mp12 /= float(sdf);
    p12 = PVector.mult(p12, mp12);
    //println("P12: ", p12);
    PVector newpos = pos1;
    for (int j=0; j < sdf; j++){
      int ind = (sdf)*i + j;
      out[ind] = newpos;
      newpos = PVector.add(newpos,p12);
    }
    
  }
}

void setup(){
  size(1280, 720);
  PVector xy = new PVector(0.0,0.0,0.0);
  Polarcoord(2.0, 1.0, xy);
  println(xy);
  getNGonPoints(NGONs1, RNG1, NG1pos);
  getNGonPoints(NGONs2, RNG2, NG2pos);
  ng12lcm = lcm(NGONs1, NGONs2);
  ng1posf = new PVector[ng12lcm];
  ng2posf = new PVector[ng12lcm];
  println(ng12lcm);
  int s1f = ng12lcm/NGONs1;
  int s2f = ng12lcm/NGONs2;
  println(s1f);
  getNGonSubdivisionPoints(s1f, NG1pos, ng1posf);
  getNGonSubdivisionPoints(s2f, NG2pos, ng2posf);
  s1 = createShape();
  s2 = createShape();
  println(ng1posf);
  createNGon(ng1posf, s1);
  createNGon(ng2posf, s2);
  s1.setStroke(255);
  s1.setFill(0,0);
}

void createNGon(PVector[] pos, PShape s){
  s.beginShape();
  for (int i = 0; i < pos.length; i++){
    s.vertex(pos[i].x, pos[i].y);
  }
  s.stroke(255);
  s.strokeWeight(5);
  s.noFill();
  s.endShape(CLOSE);
}

void getITOP(PVector[] pos1, PVector[] pos2, PVector[] itop, float atime){
  //This procedure creates an iterated measure vector between two polygon
  //states.
  for (int i=0; i < pos1.length; i++){
    PVector p12 = PVector.sub(pos2[i],pos1[i]);
    float p12m = p12.mag();
    float frames = atime*24.0;
    float p12inc = p12m/frames;
    p12.normalize();
    p12 = PVector.mult(p12, p12inc);
    itop[i] = p12;
  }
}

void updatePosPoints(PVector[] pos, PVector[] itop){
  for (int i=0; i < pos.length; i++){
    pos[i] = PVector.add(pos[i],itop[i]);
  }
}

void recurDraw(int Level, PShape s){
  //fill(int((70.0-Level)/70.0*255.0), int(nbfillA*(70.0-Level)/70.0*255.0));
  //ellipse(0,0,50.0,15.0);
  shape(s, 0,0);
  //stroke(255);
  //strokeWeight(1);
  //noStroke();
  
  //noFill();
  if (Level > 0){
    Level -= 1;
    //float x = 0;
    //float y = 0;
    //Polarcoord(angle, 100.0, x, y);
    translate(0.0,50.0);
    scale(.85);
    rotate(angle);
    recurDraw(Level, s);
  }
}

void bigRecurDraw(int Level, float iangle, PShape s){
  pushMatrix();
  rotate(iangle);
  float x = 0.0;
  float y = 0.0;
  PVector xy = new PVector(x,y,0.0);
  Polarcoord(iangle, 300.0, xy);
  translate(xy.x,xy.y);
  recurDraw(int(TopLevel), s);
  popMatrix();
  if (Level > 0){
    Level -= 1;
    iangle += 1.0*PI/float(int(TopLevel/2.0));
    nbfillA -= 1/70.0;
    if (nbfillA < 0){
      nbfillA = 0.0;
    }
    bigRecurDraw(Level, iangle, s);
    
  }
}

boolean stopacycle = true;
PVector[] pos;
PVector[] itop;
PShape s;
void draw(){
  background(0, 8);
  translate(1280.0/2.0, 720.0/2.0);
  rotate(2*angle*.01);
  //scale(angle);
  //rotate(PI/4.0);
  //shape(s1,0,0);
  //shape(s2,0,0);
  boolean t1 = (time - rtime*24) % (rtime*24.0 + atime*24.0) == 0; //start of animation cycle
  boolean t2 = (time) % (rtime*24.0 + atime*24.0) == 0; //start of rest cycle
  if (t2){
    //set destination shape
    if (shape){
      NGONs2 = int(random(3,10));
      RNG2 = random(5,100);
      NG2pos = new PVector[NGONs2];
      getNGonPoints(NGONs2, RNG2, NG2pos);
      ng12lcm = lcm(NGONs1, NGONs2);
      ng1posf = new PVector[ng12lcm];
      ng2posf = new PVector[ng12lcm];
      int s1f = ng12lcm/NGONs1;
      int s2f = ng12lcm/NGONs2;
      getNGonSubdivisionPoints(s1f, NG1pos, ng1posf);
      getNGonSubdivisionPoints(s2f, NG2pos, ng2posf);
      itop = new PVector[ng12lcm];
      //shape 2 is destination
      getITOP(ng1posf, ng2posf, itop, atime);
      pos = ng1posf.clone();
      //print("Hit shape!");
      s = createShape();
      createNGon(pos, s);
      //print(itop);
    }
    else{
      //println("Hit shape2!");
      NGONs1 = int(random(3,10));
      NG1pos = new PVector[NGONs1];
      RNG1 = random(5,100);
      getNGonPoints(NGONs1, RNG1, NG1pos);
      ng12lcm = lcm(NGONs1, NGONs2);
      ng1posf = new PVector[ng12lcm];
      ng2posf = new PVector[ng12lcm];
      int s1f = ng12lcm/NGONs1;
      int s2f = ng12lcm/NGONs2;
      getNGonSubdivisionPoints(s1f, NG1pos, ng1posf);
      getNGonSubdivisionPoints(s2f, NG2pos, ng2posf);
      itop = new PVector[ng12lcm];
      getITOP(ng2posf, ng1posf, itop, atime);
      pos = ng2posf.clone();
      s = createShape();
      createNGon(pos, s);
    }
    stopacycle = true;
    
  }
  if (t1){
    stopacycle = false;
    if (shape){
      shape = false;
    }
    else{
      shape = true;
    }
  }
  if (stopacycle != true) {
    //println("pos:", pos);
    //println("Itop: ", itop);
    updatePosPoints(pos, itop);
    //println("hitting stop acycle!");
    s = createShape();
    createNGon(pos, s);
  }
  //rotate(time*.01);
  int k = 1;
  while (k < 5){
      bigRecurDraw(int(TopLevel/2.0), 0.0, s);
  }
  if (angle > maxAngle){
    speed = -1*speed;
  }
  else if (angle <-1*maxAngle){
    speed = -1*speed;
  }
  angle += speed*.05;
  //shape(s,0.0,0.0);
  //int j = 1;
  //while (j < 30){
  //  scale(.85);
  //  shape(s,0.0,0.0);
  //  j+= 1;
  //}
  time += speed;
  //println(time);
  //saveFrame();
}