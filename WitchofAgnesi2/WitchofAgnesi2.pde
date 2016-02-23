void ltolintersect(PVector p1, PVector p2, PVector p3, PVector p4, PVector out){
  //p1 and p2 belong to L1 and p3 and p4 belong to L2
  float detp1p2 = p1.x*p2.y - p1.y*p2.x;
  float detp3p4 = p3.x*p4.y - p3.y*p4.x;
  float outxnum = detp1p2*(p3.x-p4.x) - (detp3p4*(p1.x-p2.x));
  float outdenom = (p1.x-p2.x)*(p3.y-p4.y) - (p1.y-p2.y)*(p3.x-p4.x);
  float outynum = detp1p2*(p3.y-p4.y) - (detp3p4*(p1.y-p2.y));
  out.x = outxnum/outdenom;
  out.y = outynum/outdenom;
}

PVector refVec;
PVector refVec2;
float CircleR = 300.0;
float CircleR2 = CircleR;
PVector pO; 
PVector pM;
PVector lpO2; 
PVector lpM2;
float angle = .2;
float ainc = .009;
float rinc = .001;
PVector prevP;
PVector prevP2;
PVector prevP3;
boolean start = true;
int count = 0;
ArrayList<ArrayList<PVector>> pointset = new ArrayList<ArrayList<PVector>>();
ArrayList<PVector> pts = new ArrayList<PVector>();
//ArrayList<PVector> pointset = new ArrayList<PVector>();

void setup(){
  size(1080,720);
  smooth(8);

}


void draw(){
  background(0);
  refVec = new PVector(1.0,0.0,0.0);
  refVec2 = new PVector(0.0,-1.0,0.0);
  pO = new PVector(0.0,-1.0*CircleR,0.0);
  pM = new PVector(0.0, CircleR, 0.0);
  lpO2 = PVector.add(pO, refVec);
  lpM2 = PVector.add(pM, refVec);
  int nc = int(angle/(PI));
  if (nc != count){
    //pointset = new ArrayList<ArrayList<PVector>>();
    pointset.add(pts);
    pts = new ArrayList<PVector>();
    count = nc;
    //rinc *= -1.0;
  }
  translate(1080.0/2.0, 720.0/2.0);
  scale(1.0+1.0/(CircleR*.01));
  stroke(255);
  noFill();
  ellipse(0.0,0.0, CircleR*2.0, CircleR*2.0);
  float OrigA = 2*angle;
  PVector pA = new PVector(pO.x, pO.y, 0.0);
  float sw = 10;
  //println(CircleR);
  if (CircleR < 1){
    sw = CircleR*.1;
  }
  else{
    sw -= 9.0*((300 - CircleR)/300);
  }
  strokeWeight(sw);
  
  pA.rotate(OrigA);
  point(pA.x, pA.y);
  PVector lpA2 = PVector.add(pA, refVec);
  PVector pN = new PVector(0.0,0.0,0.0);
  ltolintersect(pM, lpM2, pO, pA, pN);
  point(pN.x, pN.y);
  point(pO.x, pO.y);
  point(pM.x, pM.y);
  PVector lpN2 = PVector.add(pN,refVec2);
  PVector pP = new PVector(0.0,0.0,0.0);
  ltolintersect(pA, lpA2, pN, lpN2, pP);
  point(pP.x, pP.y);
  if (start){
    prevP = pP;
    prevP2 = pP;
    prevP3 = pP;
    start = false;
  }
  stroke(255);
  if (CircleR < 1){
    strokeWeight(sw);
  }
  else{
    strokeWeight(sw*.3);
  }
  line(pP.x, pP.y, pA.x, pA.y);
  line(pP.x, pP.y, pN.x, pN.y);
  line(pA.x, pA.y, pN.x, pN.y);
  line(pO.x, pO.y, pA.x, pA.y);
  line(pM.x, pM.y, pN.x, pN.y);
  curve(prevP3.x, prevP3.y, prevP2.x, prevP2.y, prevP.x,prevP.y, pP.x,pP.y);
  //for (ArrayList<PVector> pts : pointset){
  // curve(pts.get(0).x, pts.get(0).y, pts.get(1).x, pts.get(1).y,
  //       pts.get(2).x, pts.get(2).y, pts.get(3).x, pts.get(3).y);
  //}
  for (ArrayList<PVector> opts : pointset){
    beginShape();
    for (PVector pt : opts){
      curveVertex(pt.x, pt.y);
    }
    endShape();
  }
  beginShape();
  for(PVector pt : pts){
    curveVertex(pt.x, pt.y);
  }
  endShape();
  //ArrayList<PVector> pts = new ArrayList<PVector>();
  //PVector pts1 = new PVector(prevP.x, prevP.y,0.0);
  PVector pts2 = new PVector(pP.x, pP.y, 0.0);
  //PVector pts3 = new PVector(prevP2.x, prevP2.y, 0.0);
  //PVector pts4 = new PVector(prevP3.x, prevP3.y, 0.0);
  //pts.add(pts4);
  //pts.add(pts3);
  //pts.add(pts1);
  //pts.add(pts2);
  //pointset.add(pts);
  pts.add(pts2);
  prevP3.x = prevP2.x;
  prevP3.y = prevP2.y;
  prevP2.x = prevP.x;
  prevP2.y = prevP.y;
  prevP.x = pP.x;
  prevP.y = pP.y;
  angle += ainc;
  CircleR += -1.0*CircleR*rinc;//(cos(2.0*angle))*CircleR*rinc;
  //saveFrame();
  
}