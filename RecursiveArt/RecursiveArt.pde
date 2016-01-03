float angle = -.5;
float speed = .003;
float bigX = 0.0;
float bigY = 0.0;
float maxAngle = 2.0;
float bfillA = 1.0;
float nbfillA = 1.0;

void setup() {
  size(1280, 720);
}

void Polarcoord(float angle, float radius, float x, float y){
  x = radius*cos(angle);
  y = radius*sin(angle);
}

void recurDraw(int Level){
  fill(int((70.0-Level)/70.0*255.0), int(nbfillA*(70.0-Level)/70.0*255.0));
  ellipse(0,0,50.0,15.0);
  //stroke(255);
  //strokeWeight(1);
  noStroke();
  
  //noFill();
  if (Level > 0){
    Level -= 1;
    //float x = 0;
    //float y = 0;
    //Polarcoord(angle, 100.0, x, y);
    translate(0.0,50.0);
    scale(.95);
    rotate(angle);
    recurDraw(Level);
  }
}

void bigRecurDraw(int Level, float iangle){
  pushMatrix();
  rotate(iangle);
  float x = 0.0;
  float y = 0.0;
  Polarcoord(angle, 50.0, x, y);
  translate(x,y);
  recurDraw(70);
  popMatrix();
  if (Level > 0){
    Level -= 1;
    iangle += PI/15.0;
    nbfillA -= 1/70.0;
    if (nbfillA < 0){
      nbfillA = 0.0;
    }
    bigRecurDraw(Level, iangle);
    
  }
}

void draw() {
  background(0);
  translate(1280.0/2.0, 720.0/2);
  rotate(2*angle);
  scale(2.0*angle);
  bigRecurDraw(70, 0.0);
  if (angle > maxAngle){
    speed = -1*speed;
  }
  else if (angle <-1*maxAngle){
    speed = -1*speed;
  }
  angle += speed;
  nbfillA = bfillA;
  saveFrame();
}