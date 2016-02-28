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

ArrayList<PVector> nthroot(PVector c, Float nroot){
  //complex nth root
  ArrayList<PVector> roots = new ArrayList<PVector>();
  Float r = c.mag();
  Float theta = c.heading();
  Float rroot = pow(r, 1.0/nroot);
  for (int i = 0; i < nroot; i++){
    PVector cn = new PVector(cos((theta + i*2*PI)/nroot), sin((theta+i*2*PI)/nroot), 0.0);
    roots.add(PVector.mult(cn, rroot));
  }
  return roots;
}

PVector creq(PVector c, Float n){
  //Cuachy Riemannn equations
  //ArrayList<PVector> roots = new ArrayList<PVector>();
  Float r = c.mag();
  Float theta = c.heading();
  Float rroot = pow(r, n);
  //for (int i = 0; i < nroot; i++){
  //  PVector cn = new PVector(cos((theta + i*2*PI)/nroot), sin((theta+i*2*PI)/nroot), 0.0);
  //  roots.add(PVector.mult(cn, rroot));
  //}
  PVector cn = new PVector(cos(theta*n), sin(theta*n),0.0);
  cn = PVector.mult(cn, rroot);
  return cn;
  
}



PVector ccos(PVector z){
  //complex cos function
  // cos(z) = cos(a)cosh(b) - isin(a) sinh(b)  making use 
  // of cosh x = (e^x + e^-x)/2 for instance.
  //Float a = cos(z.x)*(exp(z.y)+exp(-z.y))/2.0;
  float a = cos(z.x)*(float)Math.cosh(z.y);
  //Float b = 1.0*(sin(z.x)*(exp(-z.y)-exp(z.y))/2.0);
  float b = -1.0*(sin(z.x)*(float)Math.sinh(z.y));
  return new PVector(a,b,0.0);
}

PVector csin(PVector z){
  //complex cos function
  // cos(z) = cos(a)cosh(b) - isin(a) sinh(b)  making use 
  // of cosh x = (e^x + e^-x)/2 for instance.
  Float a = sin(z.x)*(exp(z.y)+exp(-z.y))/2.0;
  Float b = (cos(z.x)*(exp(-z.y)-exp(z.y))/-2.0);
  return new PVector(a,b, 0.0);
}

Float cplaneinc = 9.0;
Float cplanelineinc = 10.0;
Integer dimX = 1080;
Integer dimY = 720;
Float time = 1.0;
Float timeinc = .01;
boolean popm = true;

ArrayList<ArrayList<PVector>> ConformalMap = new ArrayList<ArrayList<PVector>>(); 
ArrayList<ArrayList<PVector>> ConformalMap2 = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> ConformalMap3 = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> ConformalMap4 = new ArrayList<ArrayList<PVector>>();
ArrayList<ArrayList<PVector>> ConformalMap5 = new ArrayList<ArrayList<PVector>>();
void setup(){
  size(1080,720);
  
  PVector a = new PVector(1.0, 2.5, 0.0);
  PVector b = new PVector(1.0, -2.0, 0.0);
  println(cdiv(a,b));
  PVector z = new PVector(-dimX,-dimY,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    
    while(z.y < dimY){
      z.y += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist2.add(creq(z,-2.0));
      
    }
    z.x += cplanelineinc;
    z.y = -dimY;
    ConformalMap.add(nlist);
    ConformalMap5.add(nlist2);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist2.add(creq(z,-2.0));
    }
    z.y += cplanelineinc;
    z.x = -dimX;
    ConformalMap.add(nlist);
    ConformalMap5.add(nlist2);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.y < dimY){
      z.y += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      ArrayList<PVector> roots = nthroot(z, 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));
      //nlist.add(nroot(z, 2.0));
    }
    z.x += cplanelineinc;
    z.y = -dimY;
    ConformalMap2.add(nlist);
    ConformalMap2.add(nlist2);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      ArrayList<PVector> roots = nthroot(z, 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));     
      //nlist.add(nroot(z, 2.0));
    }
    z.y += cplanelineinc;
    z.x = -dimX;
    ConformalMap2.add(nlist);
    ConformalMap2.add(nlist2);
  }
  //println(ConformalMap);
  z = new PVector(0.0,0.0,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    while(z.y < dimY){
      z.y += cplaneinc;
      PVector sz = PVector.mult(z,1.0/50.0);
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist.add(ccos(sz));
    }
    z.x += cplanelineinc;
    z.y = 0.0;
    ConformalMap3.add(nlist);
  }
  z = new PVector(0.0,0.0,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      PVector sz = PVector.mult(z,1.0/50.0);
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      nlist.add(ccos(sz));
    }
    z.y += cplanelineinc;
    z.x = 0.0;
    ConformalMap3.add(nlist);
    //println(ConformalMap3);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.x < dimX){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.y < dimY){
      z.y += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      //nlist.add(cdiv(new PVector(1.0, 0.0, 0.0), z));
      ArrayList<PVector> roots = nthroot(cdiv(new PVector(1.0, 0.0, 0.0), z), 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));
    }
    z.x += cplanelineinc;
    z.y = -dimY;
    ConformalMap4.add(nlist);
    ConformalMap4.add(nlist2);
  }
  z = new PVector(-dimX,-dimY,0.0);
  while(z.y < dimY){
    ArrayList<PVector> nlist = new ArrayList<PVector>();
    ArrayList<PVector> nlist2 = new ArrayList<PVector>();
    while(z.x < dimX){
      z.x += cplaneinc;
      //println(z);
      //nlist.add(new PVector(z.x, z.y, 0.0));
      ArrayList<PVector> roots = nthroot(cdiv(new PVector(1.0, 0.0, 0.0), z), 2.0);
      nlist.add(roots.get(0));
      nlist2.add(roots.get(1));
    }
    z.y += cplanelineinc;
    z.x = -dimX;
    ConformalMap4.add(nlist);
    ConformalMap4.add(nlist2);
  }
}

void draw(){
  background(0);
  translate(dimX/2, dimY/2);
  
  noFill();
  stroke(255);
  strokeWeight(.5*1/time);
  //println(time);
  //strokeWeight(1.5*1/time);
  if (time < 2){
   scale(time);
   for (ArrayList<PVector> curvepts : ConformalMap){
     beginShape();
     for(PVector curvept : curvepts){
       curveVertex(10000*curvept.x, 10000*curvept.y);
     }
     endShape();
   }
   fill(255);
   textSize(12);
   text("1/z",-100.0,-100.0);  
  }
  
  else if (time < 3.0){
   scale(time);
   for (ArrayList<PVector> curvepts : ConformalMap2){
     beginShape();
     for(PVector curvept : curvepts){
       curveVertex(2.5*curvept.x, 2.5*curvept.y);
    }
    endShape();
    }
   fill(255);
   textSize(6);
   text("square root z",-100.0,-100.0);  
  }
  else if (time < 4.0){
   strokeWeight(.1*1/time);
   scale(6.0*time);
   for (ArrayList<PVector> curvepts : ConformalMap3){
     beginShape();
     for(PVector curvept : curvepts){
       curveVertex(curvept.x, curvept.y);
     }
     endShape();
   }
   fill(255);
   textSize(1);
   text("Cos(z)",-5.0,-5.0);  
  }
  else if (time < 5.0){
  
   if (popm){
     pushMatrix();
     popMatrix();
     translate(dimX/2, dimY/2);
     popm = false;
   }
   scale(time);
   for (ArrayList<PVector> curvepts : ConformalMap4){
     beginShape();
     for(PVector curvept : curvepts){
       curveVertex(400.0*curvept.x, 400.0*curvept.y);
     }
     endShape();
   }
   fill(255);
   textSize(5);
   text("1/sqroot(z)",-30.0,-30.0); 
  }
  else{
   scale(time);
   for (ArrayList<PVector> curvepts : ConformalMap5){
     beginShape();
     for(PVector curvept : curvepts){
       curveVertex(1000000.0*curvept.x, 1000000.0*curvept.y);
     }
     endShape();
   }
   fill(255);
   textSize(5);
   text("1/z^2",-30.0,-30.0); 
  }
  //scale(6.0*time);
  //for (ArrayList<PVector> curvepts : ConformalMap3){
  //  beginShape();
  //  for(PVector curvept : curvepts){
  //    curveVertex(curvept.x, curvept.y);
  //  } 
  //  endShape();
  //}
  time += timeinc;
  //saveFrame();
}