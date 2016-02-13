//Fixed Chord varies the center of the circle along the bisector line of such chord.
PVector p1 = new PVector(50.0, 150.0, 0.0);
PVector p2 = new PVector(-100.0, -200.0, 0.0);

PVector center;
PVector cp1orth;
float dinc = 1.0;
float maxradius = 0.0;
float orad = 0.0;
float radius = 0.0;
void setup(){
  size(1080,720);
  background(0);
  center = new PVector((p1.x+p2.x)/2.0, (p1.y+p2.y)/2.0, 0.0);
  PVector cp1 = PVector.sub(center,p1);
  radius = cp1.mag();
  maxradius = 2.0*radius;
  cp1.rotate(PI/2.0);
  cp1.normalize();
  cp1orth = cp1;
}

void draw(){
  fill(0,40);
  noStroke();
  rect(0,0, 1920,1080);
  translate(1080.0/2.0, 720.0/2.0);
  stroke(255);
  strokeWeight(.5);
  noFill();
  line(p1.x,p1.y,p2.x,p2.y);
  PVector nOrth = PVector.mult(cp1orth, orad);
  PVector nCenter = PVector.add(nOrth, center);
  PVector nCenterp1 = PVector.sub(nCenter, p1);
  float nr = nCenterp1.mag();
  ellipse(nCenter.x, nCenter.y, 2*nr, 2*nr);
  line(nCenter.x, nCenter.y, p1.x, p1.y);
  line(nCenter.x, nCenter.y, p2.x, p2.y);
  line(nCenter.x, nCenter.y, center.x, center.y);
  orad += dinc;
  if (abs(orad) > maxradius){
    dinc *= -1;
  }
  //saveFrame();
}