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
float CircleR = 200.0;
PVector pO; 
PVector pM;
PVector lpO2; 
PVector lpM2;
float angle = .2;
float ainc = .01;
PVector prevP;
PVector prevP2;
PVector prevP3;
boolean start = true;
int count = 0;
//ArrayList<ArrayList<PVector>> pointset = new ArrayList<ArrayList<PVector>>();
ArrayList<PVector> pointset = new ArrayList<PVector>();

void setup(){
  size(1080,720);
  smooth(8);
  refVec = new PVector(1.0,0.0,0.0);
  refVec2 = new PVector(0.0,-1.0,0.0);
  pO = new PVector(0.0,-1.0*CircleR,0.0);
  pM = new PVector(0.0, CircleR, 0.0);
  lpO2 = PVector.add(pO, refVec);
  lpM2 = PVector.add(pM, refVec);
}


void draw(){
  background(0);
  int nc = int(angle/(2.0*PI));
  if (nc != count){
    //pointset = new ArrayList<ArrayList<PVector>>();
    pointset = new ArrayList<PVector>();
    count = nc;
  }
  translate(1080.0/2.0, 720.0/2.0);
  stroke(255);
  noFill();
  ellipse(0.0,0.0, CircleR*2.0, CircleR*2.0);
  float OrigA = 2*angle;
  PVector pA = new PVector(pO.x, pO.y, 0.0);
  strokeWeight(10.0);
  
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
  strokeWeight(1.0);
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
  beginShape();
  for (PVector pt : pointset){
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
  pointset.add(pts2);
  prevP3.x = prevP2.x;
  prevP3.y = prevP2.y;
  prevP2.x = prevP.x;
  prevP2.y = prevP.y;
  prevP.x = pP.x;
  prevP.y = pP.y;
  
  angle += ainc;

}