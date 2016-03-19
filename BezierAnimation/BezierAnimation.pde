PVector[] points = {new PVector(-1.0, 0.0, 0.0), 
new PVector(1.0818086862564087, 1.4234715700149536, 0.0), 
new PVector(2.1098151206970215, 0.39557957649230957, 0.0), 
new PVector(1.8972643613815308, -1.706305742263794, 0.0), 
new PVector(-0.1573883295059204, -2.828098773956299, 0.0), 
new PVector(-2.058530807495117, -2.6982064247131348, 0.0), 
new PVector(-3.5936172008514404, -1.848006248474121, 0.0), 
new PVector(-2.9205403327941895, -0.501854658126831, 0.0)};

PVector[] lpoints = {new PVector(-1.7051784992218018, -0.32142385840415955, 0.0), new PVector(0.10144931077957153, 1.5147082805633545, 0.0), new PVector(1.919292688369751, 0.9302026629447937, 0.0), new PVector(2.3533387184143066, -1.0190889835357666, 0.0), new PVector(0.7353513240814209, -2.632430076599121, 0.0), new PVector(-1.3450032472610474, -2.908806085586548, 0.0), new PVector(-3.390713691711426, -2.5023648738861084, 0.0), new PVector(-3.376924514770508, -0.8719509243965149, 0.0)};

PVector[] rpoints = {new PVector(-0.10408270359039307, 0.40836358070373535, 0.0), new PVector(1.6469234228134155, 1.3708794116973877, 0.0), new PVector(2.386686325073242, -0.3813459277153015, 0.0), new PVector(1.3918957710266113, -2.467799663543701, 0.0), new PVector(-0.8840961456298828, -2.9873769283294678, 0.0), new PVector(-2.7156028747558594, -2.50426983833313, 0.0), new PVector(-3.76764178276062, -1.2867815494537354, 0.0), new PVector(-2.3186075687408447, -0.013728529214859009, 0.0)};

void CubicBezier(float t, PVector P1, PVector P2, PVector P3, PVector P4, PVector out){
  float c1 = pow((1-t),3);
  float c2 = 3.0*pow((1-t),2)*t;
  float c3 = 3.0*(1-t)*t*t;
  float c4 = t*t*t;
  PVector nP1 = PVector.mult(P1, c1);
  PVector nP2 = PVector.mult(P2, c2);
  PVector nP3 = PVector.mult(P3, c3);
  PVector nP4 = PVector.mult(P4, c4);
  nP1 = PVector.add(nP1,nP2);
  nP1 = PVector.add(nP1, nP3);
  nP1 = PVector.add(nP1, nP4);
  out.x = nP1.x;
  out.y = nP1.y;
}

float tstep = .0001;  //should be 0<tstep< 1
int curvelen = 500;//in point steps;
boolean showTrack = true;

ArrayList<PVector> curvepoints = new ArrayList<PVector>();
void setup(){
  size(1080,720);
  int steps = (int)(1/tstep);
  for(int i = 1; i < points.length; i++){
    int j = (i-1)%points.length;
    PVector P1 = points[j];
    PVector P4 = points[i];
    PVector P3 = lpoints[i];
    PVector P2 = rpoints[j]; 
   
    
    float pstep = 0.0;
    for (int k = 0; k < steps; k++){
      PVector out = new PVector(0.0,0.0,0.0);
      CubicBezier(pstep, P1, P2, P3, P4, out);
      curvepoints.add(out);
      pstep += tstep;
    }
  }
  PVector P1 = points[points.length-1];
  PVector P4 = points[0];
  PVector P3 = lpoints[0];
  PVector P2 = rpoints[points.length-1];
  float pstep = 0.0;
  for (int k = 0; k < steps; k++){
    PVector out = new PVector(0.0,0.0,0.0);
    CubicBezier(pstep, P1, P2, P3, P4, out);
    curvepoints.add(out);
    pstep += tstep;
  }
  println(curvepoints.size());
  
}
int astep = 0;
void draw(){
  fill(255,20);
  rect(0,0, 1920,1080);
  if (astep > (curvepoints.size()-curvelen-1)){
    astep = 0;
    println("hit");
  }
  translate(1080.0/2.0, 720.0/2.0);
  scale(100.0);
  strokeWeight(.01);
  stroke(120,40);
  noFill();
  if (showTrack){
    beginShape();
    vertex(points[0].x,points[0].y);
    for(int i = 1; i < points.length; i++){
      int j = (i-1)%points.length;
      i = i % points.length;
      float x1 = points[i].x;
      float y1 = points[i].y;
      float x2 = lpoints[i].x;
      float y2 = lpoints[i].y;
      float x3 = rpoints[j].x;
      float y3 = rpoints[j].y;
      float x4 = points[j].x;
      float y4 = points[j].y;
      bezierVertex(x3, y3, x2, y2, x1, y1);
      //bezier(x1, y1, x2, y2, x3, y3, x4, y4);
    }
    endShape();
    PVector lp1 = points[points.length-1];
    PVector lp2 = points[0];
    PVector lpc1 = lpoints[0];
    PVector lpc2 = rpoints[points.length-1];
    bezier(lp2.x,lp2.y, lpc1.x, lpc1.y, lpc2.x, lpc2.y, lp1.x, lp1.y);
  }
  int astepi = astep;
  stroke(0);
  strokeWeight(.02);
  beginShape();
  int k = 0;
  while(k < (curvelen)){
    curveVertex(curvepoints.get(astepi).x, curvepoints.get(astepi).y);
    astepi += 1;
    k += 1;
  }
  endShape();

  astep += 100;

}