/*
Conformal mapping of the complex plane z to 1/z 
*/
PVector cmult(PVector a, PVector b){
  return new PVector(a.x*b.x - (a.y*b.y), a.x*b.y + b.x*a.y, 0.0); 
}

PVector cdiv(PVector NumA, PVector DenomA){
  //get conjugate of DenomA
  PVector cNumA = new PVector(NumA.x, NumA.y, 0.0);
  PVector cDenomA = new PVector(DenomA.x,DenomA.y, 0.0);
  PVector ConjDenomA = new PVector(DenomA.x, -DenomA.y, 0.0);
  PVector pdenom = cmult(cDenomA,ConjDenomA);
  PVector pnum = cmult(cNumA,ConjDenomA);
  if (pdenom.x != 0.0){
    return PVector.mult(pnum, 1.0/pdenom.x);
  }
  else{
    return new PVector(0.0,0.0,0.0);
  }
}

Float cplaneinc = 9.0;
Float cplanelineinc = 50.0;
Integer dimX = 1080;
Integer dimY = 720;
Float time = 1.0;
Float timeinc = .01;

ArrayList<ArrayList<PVector>> ConformalMap = new ArrayList<ArrayList<PVector>>(); 
void setup(){
  size(1080,720);
  
  PVector a = new PVector(1.0, 2.5, 0.0);
  PVector b = new PVector(1.0, -2.0, 0.0);
  println(cdiv(a,b));
  PVector z = new PVector(-dimX,-dimY,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    while(z.y < dimY){
      z.y += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
    }
    z.x += cplanelineinc;
    z.y = -dimY;
    ConformalMap.add(nlist);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
    }
    z.y += cplanelineinc;
    z.x = -dimX;
    ConformalMap.add(nlist);
  }
  //println(ConformalMap);
}

void draw(){
  background(0);
  translate(dimX/2, dimY/2);
  scale(time);
  noFill();
  stroke(255);
  strokeWeight(1.5*1/time);
  for (ArrayList<PVector> curvepts : ConformalMap){
    beginShape();
    for(PVector curvept : curvepts){
      curveVertex(10000*curvept.x, 10000*curvept.y);
    }
    endShape();
  }
  time += timeinc;
  saveFrame();
}