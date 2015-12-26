FloatList c1 = new FloatList();
boolean EqualSubdivision = true;
int Levels = 7;
float p1x = 0;
float p2x = 1280;
float p1y = 0;
float p2y = 720;
int arrelement = 0;
boolean iterate = false;
boolean initialize = true;
float boundary = .35;  //set this for any value greater than 0 and less than .5.
float x1 = 0;
float x2 = 0;
float x3 = 0;
float x4 = 0;
float y1 = 0;
float y2 = 0;
float y3 = 0;
float y4 = 0;
float x5 = 0;
float y5 = 0;
float m1 = 0;
float m2 = 0;
float m3 = 0;
float m4 = 0;
float speed = 15;
FloatList rasterLines = new FloatList();

void setup() {
  size(1280, 720);
  //translate(100,50);
  drawLines(Levels,p1x,p1y,p2x,p2y,c1, boundary, EqualSubdivision);
}

void draw() {
  background(0);
  if (iterate){
    arrelement += 6;
    initialize = true;
    iterate = false;
  }
  if (initialize){
    x5 = c1.get(arrelement);
    y5 = c1.get(arrelement+1);
    x1 = x5;
    y1 = y5;
    x2 = x5;
    y2 = y5;
    x3 = x5;
    y3 = y5;
    x4 = x5;
    y4 = y5;
    m1 = c1.get(arrelement+2);
    m2 = c1.get(arrelement+3);
    m3 = c1.get(arrelement+4);
    m4 = c1.get(arrelement+5);
    initialize = false;
  }
  boolean t1 = abs(x1 - x5) < m1;
  boolean t2 = abs(x2 - x5) < m2;
  boolean t3 = abs(y3 - y5) < m3;
  boolean t4 = abs(y4 - y5) < m4;
  if (t1){
    x1 -= speed;
  }
  else{
    if (abs(x1-x5) > m1){
      x1 = x5 - m1;
    }
  }
  if (t2){
    x2 += speed;
  }
  else{
    if (abs(x2-x5) > m2){
      x2 = x5 + m2;
    }
  }
  if (t3){
    y3 -= speed;
  }
  else{
    if (abs(y3 - y5) > m3){
      y3 = y5 - m3;
    }
  }
  if (t4){
    y4 += speed;
  }
  else{
    if (abs(y4 - y5) > m4){
      y4 = y5 + m4;
    }
  }
  if (rasterLines.size() > 0) {
      for (int i = 0; i < rasterLines.size()-1; i += 10){
        float rx1 = rasterLines.get(i);
        float ry1 = rasterLines.get(i+1);
        float rx2 = rasterLines.get(i+2);
        float ry2 = rasterLines.get(i+3);
        float rx3 = rasterLines.get(i+4);
        float ry3 = rasterLines.get(i+5);
        float rx4 = rasterLines.get(i+6);
        float ry4 = rasterLines.get(i+7);
        float rx5 = rasterLines.get(i+8);
        float ry5 = rasterLines.get(i+9);
        line(rx1,ry1,rx5,ry5);
        stroke(255);
        strokeWeight(2);
        line(rx2,ry2,rx5,ry5);
        stroke(255);
        strokeWeight(2);  
        line(rx3,ry3,rx5,ry5);
        stroke(255);
        strokeWeight(2);
        line(rx4,ry4,rx5,ry5);
        stroke(255);
        strokeWeight(2);        
      }
  }
  line(x1,y1,x5,y5);
  stroke(255);
  strokeWeight(2);
  line(x2,y2,x5,y5);
  stroke(255);
  strokeWeight(2);  
  line(x3,y3,x5,y5);
  stroke(255);
  strokeWeight(2);
  line(x4,y4,x5,y5);
  stroke(255);
  strokeWeight(2);
  if (!t1 && !t2 && !t3 && !t4){
    if (arrelement+6 < c1.size()-1){
      iterate = true;
      rasterLines.append(x1);
      rasterLines.append(y1);
      rasterLines.append(x2);
      rasterLines.append(y2);
      rasterLines.append(x3);
      rasterLines.append(y3);
      rasterLines.append(x4);
      rasterLines.append(y4);
      rasterLines.append(x5);
      rasterLines.append(y5);
    }
  }
  saveFrame();
}

void drawLines(int level, float x1, float y1, float x2, float y2, FloatList c2, float boundary, boolean eqsub){
  float dist = abs(x1-x2)*boundary;
  float disty = abs(y1-y2)*boundary;
  float x3 = random(x1 + dist,x2 - dist);
  float y3 = random(y1 + disty,y2 - disty);
  float dx1 = abs(x3 - x1);
  float dx2 = abs(x3 - x2);
  float dy1 = abs(y3 - y1);
  float dy2 = abs(y3 - y2);
  c2.append(x3);
  c2.append(y3);
  c2.append(dx1);
  c2.append(dx2);
  c2.append(dy1);
  c2.append(dy2);
  int i = int(random(1,4));
  if (eqsub){
    i = 4;
  }
  int[] k1 = {1,2,3,4};
  IntList k = new IntList(k1);
  if (level > 0){
    level -= 1;
  }
  for (int j = 0; j < i; j++){
      int l = int(random(0,k.size()-1));
      int m = k.get(l);
      k.remove(l);
      float nx1 = 0;
      float nx2 = 0;
      float ny1 = 0;
      float ny2 = 0;
      if (m == 1){
        nx1 = min(x1,x3);
        nx2 = max(x1,x3);
        ny1 = min(y2,y3);
        ny2 = max(y2,y3);
      }
      else if (m == 2){
        nx1 = min(x2,x3);
        nx2 = max(x2,x3);
        ny1 = min(y2,y3);
        ny2 = max(y2,y3);
      }
      else if (m == 3){
        nx1 = min(x1,x3);
        nx2 = max(x1,x3);
        ny1 = min(y1,y3);
        ny2 = max(y1,y3);
      }
      else{
        nx1 = min(x2,x3);
        nx2 = max(x2,x3);
        ny1 = min(y1,y3);
        ny2 = max(y1,y3);
      }
      if (level > 0){
          
          drawLines(level, nx1, ny1, nx2, ny2, c2, boundary, eqsub);
      }
  }
}