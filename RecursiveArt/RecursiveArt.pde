float angle = -.5;
float speed = .003;
float bigX = 0.0;
float bigY = 0.0;

void setup() {
  size(1280, 720);
}

void Polarcoord(float angle, float radius, float x, float y){
  x = radius*cos(angle);
  y = radius*sin(angle);
}

void recurDraw(int Level){
  ellipse(0,0,100,70);
  stroke(255);
  strokeWeight(1);
  noFill();
  if (Level > 0){
    Level -= 1;
    //float x = 0;
    //float y = 0;
    //Polarcoord(angle, 100.0, x, y);
    translate(0.0,100);
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
  Polarcoord(angle, 1.0, x, y);
  translate(x,y);
  recurDraw(30);
  popMatrix();
  if (Level > 0){
    Level -= 1;
    iangle += PI/15;
    bigRecurDraw(Level, iangle);
  }
}

void draw() {
  background(0);
  translate(1280.0/2.0, 720.0/2);
  bigRecurDraw(30, 0.0);
  if (angle > .5){
    speed = -1*speed;
  }
  else if (angle <-.5){
    speed = -1*speed;
  }
  angle += speed;
  //saveFrame();
}