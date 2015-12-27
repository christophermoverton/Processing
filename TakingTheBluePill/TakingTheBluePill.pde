float amp = 1.0;
float time = 1.0;
float time2 = 1.0;
float speed = .05;
float damp = .5;
float res = 10;
float max_time = 5;
float speed2 = speed;
float level = 10;
void setup() {
  size(1280,720);
  background(0);
}

void drawEllipse(float level){
  ellipse(0,0, 500/level, 500/level);
  if (level > 0){
    level -= 1;
    drawEllipse(level);
  }
}

void draw(){
  
  fill(0,2);
  noStroke();
  rect(0,0, 1280,720);
  pushMatrix();
  translate(640, 400);
  scale(.02*time+.06*time*time);
  rotate(.2*time + .5*time*time);
  noFill();
  stroke(255-20*time, 255-20*time, 255);
  strokeWeight(.5+.2*time);
  //noStroke();
  //smooth(20);
  beginShape();
  // Exterior part of shape, clockwise winding
  for (int i = 0; i < res*4; i++){
    curveVertex(200*(damp*abs(cos(i*PI/res+i*time))+damp)*cos(i*PI/res),200*(damp*abs(sin(i*PI/res+i*time))+damp)*sin(i*PI/res));
  }
  endShape();
  popMatrix();
  translate(640, 400);
  scale(.009*time2+.01*time*time + .001*time2*time2);
  //rotate(.2*time + .5*time*time);
  //ellipse(0,0,500,500);
  drawEllipse(level);
  //// Interior part of shape, counter-clockwise winding
  //beginContour();
  //vertex(-20, -20);
  //vertex(-20, 20);
  //vertex(20, 20);
  //vertex(20, -20);
  //endContour();
  //endShape(CLOSE);
  boolean t1 = time > max_time;
  boolean t2 = time < 0;
  if (t2){
    speed = -1*speed; 
  }
  if (t1){
    speed = -1*speed;
  }
  time += speed;
  time2 += speed2;
  //if (time2 % max_time == 0){
  //  level = 2;
  //}
  
}