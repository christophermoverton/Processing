PVector[] points = {new PVector(2.798079252243042, -0.39110803604125977, 0.0), 
new PVector(1.921189308166504, -2.4318735599517822, 0.0), 
new PVector(-0.7254283428192139, -2.8942344188690186, 0.0), 
new PVector(-1.9371321201324463, -1.2520554065704346, 0.0), 
new PVector(-2.8618533611297607, 0.6771044731140137, 0.0), 
new PVector(-2.052269458770752, 2.136425256729126, 0.0), 
new PVector(1.0797173976898193, 1.4189691543579102, 0.0)};

PVector[] lpoints = {new PVector(2.6267125606536865, 0.5681036114692688, 0.0), 
new PVector(2.601461410522461, -1.894078016281128, 0.0), 
new PVector(0.2482459545135498, -3.284349203109741, 0.0), 
new PVector(-1.5265755653381348, -1.9348978996276855, 0.0), 
new PVector(-2.8867032527923584, -0.15774792432785034, 0.0), 
new PVector(-2.6472978591918945, 1.871018886566162, 0.0), 
new PVector(-0.009972095489501953, 2.040407657623291, 0.0)};

PVector[] rpoints = {new PVector(2.9505887031555176, -1.244767665863037, 0.0), 
new PVector(1.0983455181121826, -3.082380533218384, 0.0), 
new PVector(-1.46503484249115, -2.5979018211364746, 0.0), 
new PVector(-2.3675060272216797, -0.5362522006034851, 0.0), 
new PVector(-2.842468500137329, 1.3283523321151733, 0.0), 
new PVector(-0.9066309928894043, 2.647425651550293, 0.0), 
new PVector(1.9261479377746582, 0.9362587928771973, 0.0)};

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
boolean showTrack = false;

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
    println(k);
    curveVertex(curvepoints.get(astepi).x, curvepoints.get(astepi).y);
    astepi += 1;
    k += 1;
  }
  endShape();

  astep += 100;

}