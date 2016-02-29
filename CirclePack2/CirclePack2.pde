import java.util.Map;
import java.util.*;

/*
Compute circle packings according to the Koebe-Thurston-Andreev theory,
Following a numerical algorithm by C. R. Collins and K. Stephenson,
"A Circle Packing Algorithm", Comp. Geom. Theory and Appl. 2003.
*/
Integer pass = 1;
float Eradius = 30.0;
float Efactor = 1.3;
float tolerance  = 1.0+1.0e-12;
Integer iterationMax = 8000;
ArrayList<ArrayList<Integer>> ComplexFamily = new ArrayList<ArrayList<Integer>>();
ArrayList<ArrayList<Circle>> CirclesFamily = new ArrayList<ArrayList<Circle>>();
ArrayList<ArrayList<ArrayList<Integer>>> ComplexFamilyPacks = new ArrayList<ArrayList<ArrayList<Integer>>>(); 


class Edge{
  PVector p1;  // p1 is point 1 on the edge pair
  PVector p2;  // p2 is point 2 on the edge pair
  int p1index;  // polygon vertices index
  int p2index;  
  PVector genCenter;  //generator center.  This is the edge coordinate generating center
  // genCenter is the main generating circle for such edge
  float genR; // radius of genrating circle.
  float angle1; //polar angle1 generating p1
  float angle2; //polar angle2 generating p2
  float rangle1; //rendering angles for arc rendering in processing
  float rangle2; //rendering angles for arc rendering in processing
  boolean thick; //Thick edge is true otherwise thin edge is false
  boolean linear;  // Linear if true otherwise arc if false
  Circle circle;
  int interiornp1;  //this is only called for main edge to interior edges mappings
  int interiornp2; //this is only called for main edge to interior edges mappings 
  int pole1;  //poles are the original vertices on the polygon prior to subdivision
              //This is a subdivision index mapping.
  int pole2;  //used in special case edge mappings
  boolean pole;  //flagged edge that wraps around a pole
  Edge(PVector P1, PVector P2, boolean Thick, int P1index, int P2index){
    //This is a assumed linear edge constructor call...use this constructor call
    //only if the edge is linear
    p1 = P1;
    p2 = P2;
    linear = true;
    thick = Thick;
    p1index = P1index;
    p2index = P2index;
  }
  Edge(PVector P1, PVector P2, boolean Thick, PVector GenCenter,
       float GenR, float Angle1, float Angle2, int P1index, int P2index){
    //This is a assumed linear edge constructor call...use this constructor call
    //only if the edge is linear
    p1 = P1;
    p2 = P2;
    linear = false;
    thick = Thick;
    genCenter = GenCenter;
    genR = GenR;
    angle1 = Angle1;
    angle2 = Angle2;
    rangle1 = angle1;
    rangle2 = angle2;
    p1index = P1index;
    p2index = P2index;
    circle = new Circle(GenCenter,GenR);
  }
  Edge(PVector P1, PVector P2, PVector GenCenter,
       float GenR, int P1index, int P2index){
    //This is a special edge constructor call used for edge Circle to Interior
    // Circles mapping. That is in forming the subdivisions of an existing
    // polygon.
    p1 = P1;
    p2 = P2;
    linear = false;
    //thick = Thick;
    genCenter = GenCenter;
    genR = GenR;
    //angle1 = Angle1;
    //angle2 = Angle2;
    p1index = P1index;
    p2index = P2index;
    //interiornp1 = Interiornp1;
    //interiornp2 = Interiornp2;
    circle = new Circle(GenCenter,GenR);
  }
  Edge(Circle circlei){
    linear = false;
    genCenter = circlei.center;
    genR = circlei.radius;
    circle = circlei;
    angle1 = 2.0*PI;
    angle2 = 0.0;
    rangle1 = angle1;
    rangle2 = angle2;
  }
  Edge(Circle circlei, float Angle1, float Angle2){
    linear = false;
    genCenter = circlei.center;
    genR = circlei.radius;
    circle = circlei;
    angle1 = Angle1;
    angle2 = Angle2;
    rangle1 = angle1;
    rangle2 = angle2;
  }
  Edge(){
    p1 = new PVector(0.0,0.0,0.0);
    p2 = new PVector(0.0,0.0,0.0);  // p2 is point 2 on the edge pair
    p1index = 0;  // polygon vertices index
    p2index = 0;  
    genCenter = new PVector(0.0,0.0,0.0);  //generator center.  This is the edge coordinate generating center
  // genCenter is the main generating circle for such edge
    genR = 1.0; // radius of genrating circle.
    angle1 = 0.0; //polar angle1 generating p1
    angle2 = 0.0; //polar angle2 generating p2
    rangle1 = angle1;
    rangle2 = angle2;
    thick = false; //Thick edge is true otherwise thin edge is false
    linear = false;  // Linear if true otherwise arc if false
    circle = new Circle(genCenter,genR);
    interiornp1 = 0;  //this is only called for main edge to interior edges mappings
    interiornp2 = 0; //this is only called for main edge to interior edges mappings 
    pole1 = 0;  //poles are the original vertices on the polygon prior to subdivision
              //This is a subdivision index mapping.
    pole2 = 0;  //used in special case edge mappings
    pole = false;  //flagged edge that wraps around a pole
  }
  void computeAngles(){
    PVector angleGen1 = PVector.sub(p1,genCenter);
    PVector angleGen2 = PVector.sub(p2,genCenter);
    
    angle1 = angleGen1.heading();
    angle2 = angleGen2.heading();
    rangle1 = angle1;
    rangle2 = angle2;
    if (rangle1 < 0.0){
      rangle1 += 2.0*PI;
    }
    if (rangle2 < 0.0){
      rangle2 += 2.0*PI;
    }
    if (angle1 < 0.0){
      angle1 += 2.0*PI;
    }
    if (angle2 < 0.0){
      angle2 += 2.0*PI;
    }

    boolean t1 = rangle2 > rangle1;
    if (t1){
      rangle1 += 2.0*PI;
    }
  }
  void computeRadius(){
    PVector angleGen1 = PVector.sub(p1,genCenter);
    genR = angleGen1.mag();
  }
  float computeAnglediff(){
    float a12 = angle2 - angle1;
    float sign = 1.0;
    if (a12 < 0){
      sign *= -1.0;
    }
    if (abs(a12) >= PI){    
      return -1.0*sign*(2.0*PI - abs(a12));
    }
    else{
      return a12;
    }
  }
  
  PVector getEdgeHalfAnglePt(){
    PVector p1pc = PVector.sub(p1, genCenter);
    float anglediff = angle2 - angle1;
    anglediff *= .5;
    p1pc.rotate(anglediff);
    return PVector.add(p1pc,genCenter);
  }
}

class Circle{
  PVector center;
  float radius;
  Circle(PVector Center, float Radius){
    center = Center;
    radius = Radius;
  }
  Circle(){
    center = new PVector(0.0,0.0);
    radius = 0.0;
  }
}

float acxyz(float x, float y, float z){
    //Angle at a circle of radius x given by two circles of radii y and z"""
    try{
        if ((2.0*(x+y)*(x+z))== 0){
          return PI;
        }
        return acos((pow((x+y),2)+pow((x+z),2)-pow((y+z),2))/(2.0*(x+y)*(x+z)));
    }
    catch(IllegalArgumentException e){
        return PI/3.0;
    }
}

float flower(HashMap<Integer, Float> radius, Integer center, ArrayList<Integer> cycle){
    //"""Compute the angle sum around a given internal circle"""
    float sum = 0.0;
    for (int i = 0; i < cycle.size(); i++){
      int j = i-1;
      if (j == -1){
        j = cycle.size()-1;
      }
      //println(j);
      //println(i);
      //println(center);
      //println(radius);
      sum += acxyz(radius.get(center), radius.get(cycle.get(j)), radius.get(cycle.get(i)));
    }
    return sum;
    /*return sum(acxyz(radius[center],radius[cycle[i-1]],radius[cycle[i]])
               for i in range(len(cycle)))
    */
}

void place(HashMap<Integer,PVector> placements,HashMap<Integer,Float> radii,
           HashMap<Integer,ArrayList<Integer>> internal, Integer center){
    //"""Recursively find centers of all circles surrounding k"""
    if (!internal.containsKey(center)){
        return;
    }
    ArrayList<Integer> cycle = internal.get(center);
    for (int i = -1 * cycle.size(); i < cycle.size(); i++){
    //for i in range(-len(cycle),len(cycle)-1):
        //if cycle[i] in placements and cycle[i+1] not in placements:
        int ii = i;
        if (i < 0){
          ii = i%cycle.size();
          if (ii != 0){
            ii += cycle.size();
          }
        }
        //println("i in place", i);
        int j = (ii+1)%cycle.size();
        if (placements.containsKey(cycle.get(ii)) && !placements.containsKey(cycle.get(j))){
            Integer s = cycle.get(ii);
            Integer t = cycle.get(j);
            float theta = acxyz(radii.get(center),radii.get(s),radii.get(t));
            PVector centTos = PVector.sub(placements.get(s),placements.get(center));
            centTos.normalize();
            centTos.rotate(theta);
            
            //PVector offset = new PVector(0.0,0.0,0.0);
            //offset.x = (placements.get(s).x - placements.get(center).x)/(radii.get(s)+ radii.get(center));
            //offset.y = (placements.get(s).y - placements.get(center).y)/(radii.get(s)+ radii.get(center));
            //offset = (placements[s]-placements[center])/(radii[s]+radii[center]);
            //offset *= e**(-1j*theta);
            //offset.x *= cos(theta);
            //offset.y *= -1.0*sin(theta);
            //println(offset);
            //println(theta);
            PVector tvec = PVector.mult(centTos, radii.get(t)+radii.get(center));
            tvec = PVector.add(tvec, placements.get(center));
            //tvec.x = placements.get(center).x + offset.x*(radii.get(t)+radii.get(center));
            //tvec.y = placements.get(center).y + offset.y*(radii.get(t)+radii.get(center));
            placements.put(t,tvec);
            //placements[t] = placements[center] + offset*(radii[t]+radii[center]);
            place(placements,radii,internal,t);
        }
    }
}

void circlePack(HashMap<Integer,ArrayList<Integer>> internal, 
                HashMap<Integer,Float> external, HashMap<Integer,Circle> out){
  /*
  Find a circle packing for the given data.
    The two arguments should be HashMaps with disjoint keys; the
    keys of the two arguments are identifiers for circles in the packing.
    The internal argument maps each internal circle to its cycle of
    surrounding circles; the external argument maps each external circle
    to its desired radius. The return function is a mapping from circle
    keys to pairs (center,radius) where center is a complex number.
  */
  HashMap<Integer,Float> radii = new HashMap<Integer,Float>(external);
  for (Map.Entry<Integer,ArrayList<Integer>> k : internal.entrySet()){
    if (external.containsKey(k.getKey())){
      println ("CirclePack internal and external keys need to be disjoint.");
    }
    radii.put(k.getKey(), 1.0);
  }
  float lastChange = 2.0;
  Integer iterat = 0;
  while (lastChange > tolerance){
    lastChange = 1;
    if (iterat > iterationMax){
      break;
    }
    for (Map.Entry<Integer,ArrayList<Integer>> k : internal.entrySet()){
      //println(k.getValue());
      float theta = flower(radii, k.getKey(), k.getValue());
      float hat = radii.get(k.getKey())/(1/sin(theta/(2*(k.getValue()).size()))-1);
      float newrad = hat * (1/(sin(PI/(k.getValue()).size())) - 1);
      float kc = max(newrad/radii.get(k.getKey()),radii.get(k.getKey())/newrad);
      lastChange = max(lastChange,kc);
      radii.put(k.getKey(),newrad);
    }
    iterat += 1;
  }
  println("finished iteration");
  HashMap<Integer,PVector> placements = new HashMap<Integer,PVector>();
  Iterator it = internal.entrySet().iterator();
  Map.Entry<Integer,ArrayList<Integer>> k1 = (Map.Entry)it.next();
  placements.put(k1.getKey(),new PVector(0.0,0.0,0.0));
  Integer k2 = k1.getValue().get(0);
  placements.put(k2, new PVector(radii.get(k1.getKey())+radii.get(k2), 0.0, 0.0));
  place(placements, radii, internal, k1.getKey());
  //place(placements, radii, internal, k2);
  for (Map.Entry<Integer, Float> k : radii.entrySet()){
    out.put(k.getKey(), new Circle(placements.get(k.getKey()), k.getValue()));
  }
}

HashMap<Integer,Circle> out = new HashMap<Integer,Circle>();

float computeDet2Dmat(float[][] iMat){
  return iMat[0][0]*iMat[1][1] - (iMat[1][0]*iMat[0][1]);
}

float computeDet3Dmat(float[][] iMat){
  float[][] cofMat1 = {{iMat[1][1], iMat[1][2]},{iMat[2][1],iMat[2][2]}};
  float[][] cofMat2 = {{iMat[1][0], iMat[1][2]},{iMat[2][0],iMat[2][2]}};
  float[][] cofMat3 = {{iMat[1][0], iMat[1][1]},{iMat[2][0],iMat[2][1]}};
  float det1 = computeDet2Dmat(cofMat1);
  float det2 = computeDet2Dmat(cofMat2);
  float det3 = computeDet2Dmat(cofMat3);
  return iMat[0][0]*det1 -1.0*iMat[0][1]*det2 + iMat[0][2]*det3;
}

void CircumCircleCenter(PVector p1, PVector p2, PVector p3, PVector out){
  float[][] mata = {{p1.x,p1.y, 1.0},{p2.x,p2.y,1.0},{p3.x,p3.y,1.0}};
  float[][] matbx = {{pow(p1.x,2.0)+pow(p1.y,2.0), p1.y,1.0},
                     {pow(p2.x,2.0)+pow(p2.y,2.0), p2.y,1.0},
                     {pow(p3.x,2.0)+pow(p3.y,2.0), p3.y,1.0}};
  float[][] matby = {{pow(p1.x,2.0)+pow(p1.y,2.0), p1.x,1.0},
                     {pow(p2.x,2.0)+pow(p2.y,2.0), p2.x,1.0},
                     {pow(p3.x,2.0)+pow(p3.y,2.0), p3.x,1.0}};
  float a = computeDet3Dmat(mata);
  float bx = -1.0*computeDet3Dmat(matbx);
  float by = 1.0*computeDet3Dmat(matby);
  out.x = -1.0*(bx)/(2.0*a);
  out.y = -1.0*(by)/(2.0*a);
}

void CircleCircleIntersection(Circle c1, Circle c2, ArrayList<PVector> ipts){
  PVector c1c2 = PVector.sub(c1.center, c2.center);
  float c1c2angle = c1c2.heading();
  float dc1c2 = c1c2.mag();
  float x = (dc1c2*dc1c2 - c1.radius*c1.radius + c2.radius*c2.radius)/(2*dc1c2);
  float y2 = (4*pow(dc1c2,2)*pow(c2.radius,2) -pow(pow(dc1c2,2)-pow(c1.radius,2)+pow(c2.radius,2),2))/(4*pow(dc1c2,2));
  //ipoint.x = x;
  //ipoint.y = y;
  float yn = -1*pow(y2,.5);
  float yp = pow(y2,.5);
  PVector xyn = PVector.add((new PVector(x,yn,0.0)).rotate(c1c2angle),c2.center);
  PVector xyp = PVector.add((new PVector(x,yp,0.0)).rotate(c1c2angle),c2.center);
  ipts.add(xyn);
  ipts.add(xyp);
}

ArrayList<ArrayList<Integer>> subdivide(ArrayList<ArrayList<Integer>> complexfamily,
               HashMap<Integer, ArrayList<Integer>> internal,
               HashMap<Integer, Float> external, Integer[] maxlabel){
  // 5 point subdivision rule...we track labels only not points or positions.
  ArrayList<ArrayList<Integer>> nout = new ArrayList<ArrayList<Integer>>();
  HashMap<Integer[], Integer> subdivEdges = new HashMap<Integer[],Integer>();
  HashMap<ArrayList<Integer>, Integer> subdivEdges2 = new HashMap<ArrayList<Integer>, Integer>();
  for (ArrayList<Integer> complex : complexfamily){
    //create new 5 complex
    ArrayList<Integer> inner5 = new ArrayList<Integer>();
    Integer[] isubdivlabels = {0,0,0,0,0};
    maxlabel[0] += 1;
    isubdivlabels[0] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[1] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[2] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[3] = maxlabel[0];
    maxlabel[0] += 1;
    isubdivlabels[4] = maxlabel[0];
    maxlabel[0] += 1;
    Integer[] subdivlabels = {0,0,0,0,0};
    Boolean[] internallabel = {false,false,false,false,false};
    //get subdivlabel pts set to subdivlabel array

    //rewrite the original internal cycle.  Done once only...
    for (int i = 0; i < complex.size(); i++){
      int j = i-1;
      if (j == -1){
        j = complex.size()-1;
      }
      Integer nlabel = complex.get((i+1)%complex.size());
      Integer ilabel = complex.get(i);
      Integer[] lpair = {nlabel,ilabel};
      Integer[] lpair2 = {ilabel, nlabel};
      ArrayList<Integer> lpair2al = new ArrayList<Integer>();
      Collections.addAll(lpair2al, lpair2);
      ArrayList<Integer> lpairal = new ArrayList<Integer>();
      Collections.addAll(lpairal, lpair);
      if (subdivEdges2.containsKey(lpairal)){
        
        subdivlabels[j] = subdivEdges2.get(lpairal);

      }
      else if (subdivEdges2.containsKey(lpair2al)){
        
        subdivlabels[j] = subdivEdges2.get(lpair2al);
      }
      else{
        maxlabel[0] += 1;
        subdivlabels[j] = maxlabel[0];
        subdivEdges.put(lpair2, maxlabel[0]);
        subdivEdges2.put(lpair2al, maxlabel[0]);
      }
      if (internal.containsKey(nlabel)){
        ArrayList<Integer> cycle = internal.get(nlabel);
        Integer ilcind = cycle.indexOf(ilabel);
        if (ilcind != -1){
          cycle.set(ilcind, subdivlabels[j]);
        }
      }
      if (internal.containsKey(ilabel)){
        ArrayList<Integer> cycle = internal.get(ilabel);
        Integer nlcind = cycle.indexOf(nlabel);
        if (nlcind != -1){
          cycle.set(nlcind, subdivlabels[j]);
        }
      }

      if (internal.containsKey(subdivlabels[j])){
        internal.get(subdivlabels[j]).add(isubdivlabels[i]);
      }
      else{
        if (!external.containsKey(ilabel)||!external.containsKey(nlabel)){
          //set internal
          ArrayList<Integer> cycle = new ArrayList<Integer>();
          cycle.add(nlabel);
          cycle.add(isubdivlabels[i]);
          cycle.add(ilabel);
          internal.put(subdivlabels[j],cycle);
        }
        else{
          external.put(subdivlabels[j], Eradius);
        }
      }
      
    }
    //pass to write internal points
    for (int i = 0; i < isubdivlabels.length; i++){
      ArrayList<Integer> iComplex = new ArrayList<Integer>();
      Integer label1 = isubdivlabels[(i+1)%isubdivlabels.length];
      int pi = i - 1;
      if (i == 0){
        pi = isubdivlabels.length-1;
      }
      Integer label2 = isubdivlabels[pi];
      Integer label3 = subdivlabels[pi];
      ArrayList<Integer> cycle = new ArrayList<Integer>();
      Collections.addAll(cycle, new Integer[] {label1,label2,label3});
      internal.put(isubdivlabels[i], cycle);
      label1 = isubdivlabels[i];
      label2 = subdivlabels[(pi)%5];
      label3 = complex.get((i+1)%5);
      if (pass == 1){
        label3 = complex.get(i);
      }
      Integer label4 = subdivlabels[(i)%subdivlabels.length];
      Integer label5 = isubdivlabels[(i+1)%isubdivlabels.length];
      ArrayList<Integer> ncomplex = new ArrayList<Integer>();
      Collections.addAll(ncomplex, new Integer[] {label1,label2,label3,label4,label5});
      inner5.add(isubdivlabels[i]);
      nout.add(ncomplex);
      int j = i -1;
      if (j == -1){
        j = 4;
      }

    }
    nout.add(inner5);
    //Integer label6 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label7 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label8 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label9 = maxlabel[0];
    //maxlabel[0] += 1;
    //Integer label10 = maxlabel[0];
    
  }
  return nout;
}

void getEradius(Boolean polytype, Float inputR, Float[] out){
  //if polytype is 4gon then inputR is 1/3 sidelength and not 1/2*sidelength as
  ///in the case of 5gon (false).
  if (polytype){
    out = new Float[3];
    Float sidelength = abs(3.0*inputR/2.0*(1+Efactor)); //on subdiv 4 gon sl
    Float sidelength2 = abs(inputR*Efactor-sidelength);
    Float sidelength3 = abs(abs(3.0*inputR/2.0 - Efactor*inputR/2.0)-sidelength);
    out[0] = sidelength/2.0;
    out[1] = sidelength2/3.0;
    out[2] = sidelength3/2.0;
  }
  else{
    out = new Float[2];
    Float sidelength =  abs(2.0*inputR-2.0*inputR*Efactor);
    Float sidelength2 = inputR*Efactor-sidelength;
    out[0] = sidelength/3.0;
    out[1] = sidelength2/2.0;
  }
}

void triangulateComplex2(HashMap<Integer, ArrayList<Integer>> internal, 
                         HashMap<Integer, Float> external,
                        Integer[] maxlabel){
    Integer tlabelSize = ComplexFamily.size()+1;
    ArrayList<Integer> newLabels = new ArrayList<Integer>();
    ArrayList<ArrayList<Integer>> ncomplexfamily = new ArrayList<ArrayList<Integer>>();
    
    for (int h = 0; h < ComplexFamily.size(); h++){
      if (ComplexFamily.get(h).size() == 5){
        //25 labels plus centroids on 5 gon subdivison 
        Integer[] outersubdiv = {7,8,9,10,11,12,13,14,15,16};
        Integer[] innersubdiv1 = {17,18,19,20,21,22,23,24,25,26};
        Integer[] innersubdiv2 = {27,28,29,30,31};
        Integer[] outercentroiddiv = {32,33,34,35,36,37,38,39,40,41,42,43,44,45,46};
        Integer[] innercentroiddiv = {47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,
                                      62,63,64,65,66,67,68,69,70,71};
        Integer[] outersubdivmid4gon = {72,73,74,75,76,77,78,79,80,81};
        Integer[] innersubdivmidgon = {82,83,84,85,86,87,88,89,90,91};
        Integer[] mapoutertocent = {0,1,3,5,6,8,9,11,12,14};
        Integer[] mapinnertoicent = {0,3,5,8,10,13,15,18,20,23};
        Integer[] mapinnertoinner2 = {0,0,1,1,2,2,3,3,4,4};
        Float[] divEradii = {0.0,0.0,0.0};
        getEradius(false, Eradius/2.0, divEradii);
        for (int i = 0; i < outersubdiv.length; i++){
          ArrayList<Integer> ncycle = new ArrayList<Integer>();
          Integer im = mapoutertocent[i];
          Integer jm = (im+1)%outercentroiddiv.length;
          if (i % 2 == 0){
            Integer[] pol = {0,0,0,0,0,0,0,0};
            //write outer 8 polys for 4gon labeling
            Integer j = (i+1)%outersubdiv.length;
            pol[0] = outersubdiv[i];
            pol[1] = outersubdiv[j];
            pol[2] = innersubdiv1[i];
            pol[3] = innersubdiv1[j];
            pol[4] = outersubdivmid4gon[i];
            pol[5] = outersubdivmid4gon[j];
            pol[6] = outercentroiddiv[im];
            pol[7] = outercentroiddiv[jm];
            Collections.addAll(ncycle,pol);
            ncomplexfamily.add(ncycle);
            Integer[] pol1 = {pol[0],pol[4],pol[6]};
            Integer[] pol2 = {pol[0],pol[6],pol[2]};
            Integer[] pol3 = {pol[2],pol[6],pol[5]};
            Integer[] pol4 = {pol[4],pol[7],pol[6]};
            Integer[] pol5 = {pol[5],pol[6],pol[7]};
            Integer[] pol6 = {pol[4],pol[1],pol[7]};
            Integer[] pol7 = {pol[1],pol[3],pol[7]};
            Integer[] pol8 = {pol[3],pol[5],pol[7]};
            Integer[] icycle1 = {pol[0],pol[6],pol[5]}; // key pol[2]
            Integer[] icycle2 = {pol[0],pol[4],pol[7], pol[5], pol[2]}; //key pol[6]
            Integer[] icycle3 = {pol[1],pol[3],pol[5], pol[6], pol[4]};//key pol[7]
            Integer[] icycle4 = {pol[2],pol[6],pol[7], pol[3]}; //key pol[5]
            Integer[] icycle5 = {pol[5], pol[7], pol[1]}; //key pol[3]
            Integer[][] icycles = {icycle1,icycle2,icycle3,icycle4,icycle5};
            Integer[] cycleskeys = {pol[2],pol[6],pol[7],pol[5],pol[3]};
            for (int k = 0; k < 5; k++){
              ArrayList<Integer> icycle = new ArrayList<Integer>();
              if  (internal.containsKey(cycleskeys[k])){
                ArrayList<Integer> ocycle = internal.get(cycleskeys[k]);
                icycle.addAll(ocycle);
              }
              Collections.addAll(icycle,icycles[k]);
              internal.put(cycleskeys[k],icycle);
            }
            external.put(pol[0], divEradii[0]);
            external.put(pol[4], divEradii[0]/2.0);
            external.put(pol[1], divEradii[0]);
          }
          else{
            ncycle = new ArrayList<Integer>();
            Integer v = i/2;
            Integer[] pol = {0,0,0,0,0,0};
            //write outer 8 polys for 4gon labeling
            Integer j = (i+1)%outersubdiv.length;
            pol[0] = outersubdiv[i];
            pol[1] = v;
            pol[2] = outersubdiv[j];
            pol[3] = innersubdiv1[i];
            pol[4] = innersubdiv1[j];
            pol[5] = outercentroiddiv[im];
            Collections.addAll(ncycle, pol);
            ncomplexfamily.add(ncycle);
            Integer[] icycle1 = {pol[0],pol[1],pol[2],pol[4],pol[3]}; //key pol[5]
            ArrayList<Integer> icycle = new ArrayList<Integer>();
            Collections.addAll(icycle,icycle1);
            internal.put(pol[5], icycle);
            internal.get(pol[3]).add(pol[4]);
            if (internal.containsKey(pol[4])){
              internal.get(pol[4]).add(0,pol[3]);
            }
            else{
              Integer[] icycle2 = {pol[3],pol[2]};
              icycle = new ArrayList<Integer>();
              Collections.addAll(icycle, icycle2);
              internal.put(pol[4],icycle);
            }
            external.put(pol[0], divEradii[1]);
            external.put(pol[1], divEradii[1]);
            external.put(pol[2], divEradii[1]);
          }
        }
        //write next inner level loop of polygons 
        for(int i = 0; i < innersubdiv1.length; i++){
          Integer im = mapinnertoicent[i];
          Integer jm = (im+1)%innercentroiddiv.length;
          if (i % 2 == 0){
            Integer j = (i+1)%innersubdiv1.length;
            Integer km = (im+2)%innercentroiddiv.length;
            Integer in = mapinnertoinner2[i];
            //write 3gon 
            Integer[] pol = {0,0,0,0,0,0,0,0,0};
            //write outer 8 polys for 4gon labeling
           
            pol[0] = innersubdiv1[i];
            pol[1] = innersubdiv1[j];
            pol[2] = innersubdiv2[in];
            pol[3] = outersubdivmid4gon[i];
            pol[4] = innersubdivmidgon[i];
            pol[5] = innersubdivmidgon[j];
            pol[6] = innercentroiddiv[im];
            pol[7] = innercentroiddiv[jm]; 
            pol[8] = innercentroiddiv[km]; //pol[8] -> {pol[4],pol[6],pol[7],pol[5],pol[2]}
            ArrayList<Integer> ncycle = new ArrayList<Integer>();
            Collections.addAll(ncycle, pol);
            ncomplexfamily.add(ncycle);
            ncycle = new ArrayList<Integer>();
            Collections.addAll(internal.get(pol[0]), new Integer[] {pol[6],pol[4]});
            internal.get(pol[1]).add(0,pol[5]);
            internal.get(pol[1]).add(0,pol[7]);
            if (internal.containsKey(pol[2])){
              Collections.addAll(internal.get(pol[2]), new Integer[] {pol[8],pol[5]});
            }
            else{
              ncycle = new ArrayList<Integer>();
              Collections.addAll(ncycle, new Integer[]{pol[4],pol[8],pol[5]});
              internal.put(pol[2],ncycle);
            }
            if (internal.containsKey(pol[3])){
              ArrayList insertcycle = new ArrayList<Integer>();
              Collections.addAll(insertcycle, new Integer[] {pol[7],pol[6]});
              insertcycle.addAll(internal.get(pol[3]));
              internal.put(pol[3], insertcycle);
            }
            else{
              ncycle = new ArrayList<Integer>();
              Collections.addAll(ncycle, new Integer[]{pol[1],pol[7],pol[6],pol[0]});
              internal.put(pol[3],ncycle);
            }
            if (internal.containsKey(pol[5])){
              ArrayList insertcycle = new ArrayList<Integer>();
              Collections.addAll(insertcycle, new Integer[] {pol[8],pol[7]});
              insertcycle.addAll(internal.get(pol[5]));
              internal.put(pol[5], insertcycle);
            }
            else{
              ncycle = new ArrayList<Integer>();
              Collections.addAll(ncycle, new Integer[]{pol[2],pol[8],pol[7],pol[1]});
              internal.put(pol[5],ncycle);
            }
            if (internal.containsKey(pol[4])){
              Collections.addAll(internal.get(pol[4]), new Integer[] {pol[7],pol[6]});
            }
            else{
              ncycle = new ArrayList<Integer>();
              Collections.addAll(ncycle, new Integer[]{pol[1],pol[7],pol[6],pol[0]});
              internal.put(pol[4],ncycle);
            }
            ncycle = new ArrayList<Integer>();
            //pol[6] -> {pol[0],pol[3],pol[7],pol[8],pol[4]}
            Collections.addAll(ncycle, new Integer[]{pol[0],pol[3],pol[7],pol[8],pol[4]});
            internal.put(pol[6], ncycle);
            ncycle = new ArrayList<Integer>();
            //pol[7] -> {pol[3], pol[1], pol[5], pol[8],pol[6]}
            Collections.addAll(ncycle, new Integer[]{pol[3], pol[1], pol[5], pol[8],pol[6]});
            internal.put(pol[7], ncycle);
            ncycle = new ArrayList<Integer>();
            //pol[8] -> {pol[4],pol[6],pol[7],pol[5],pol[2]}
            Collections.addAll(ncycle, new Integer[]{pol[4],pol[6],pol[7],pol[5],pol[2]});
            internal.put(pol[8], ncycle);
          }
          else{
            Integer j = (i+1)%innersubdiv1.length;
            Integer km = (im+2)%innercentroiddiv.length;
            Integer in = mapinnertoinner2[i];
            Integer jn = (in+1)%innersubdiv2.length;
            //write 3gon 
            Integer[] pol = {0,0,0,0,0,0,0,0,0};
            Integer[] pol2 = {0,0,0,0,0,0,0,0,0};
            //write outer 8 polys for 4gon labeling
           
            pol[0] = innersubdiv1[i];
            pol[1] = innersubdiv1[j];
            pol[2] = innersubdiv2[in];
            pol[3] = innersubdiv2[jn];
            pol[4] = innersubdivmidgon[i];
            pol[5] = innersubdivmidgon[j];
            pol[6] = innercentroiddiv[im];
            pol[7] = innercentroiddiv[jm];
            pol2[0] = pol[2];
            pol2[1] = pol[0];
            pol2[2] = pol[3];
            pol2[3] = pol[1];
            pol2[4] = pol[4];
            pol2[5] = pol[5];
            pol2[6] = pol[6];
            pol2[7] = pol[7];
            
            ArrayList<Integer> ncycle = new ArrayList<Integer>();
            Collections.addAll(ncycle, pol2);
            ncomplexfamily.add(ncycle);
            internal.get(pol[2]).add(pol[6]);
            ncycle = new ArrayList<Integer>();
            Integer[] ncyclearr = {pol[2], pol[4], pol[7], pol[5], pol[3]};
            Collections.addAll(ncycle, ncyclearr);
            internal.put(pol[6],ncycle);
            ncycle = new ArrayList<Integer>();
            ncyclearr = new Integer[] {pol[0], pol[1], pol[5], pol[6], pol[4]};
            Collections.addAll(ncycle, ncyclearr);
            internal.put(pol[7], ncycle);
            internal.get(pol[0]).add(pol[7]);
            internal.get(pol[1]).add(0,pol[5]);
            internal.get(pol[1]).add(0,pol[7]);
            Collections.addAll(internal.get(pol[4]), new Integer[]{pol[6],pol[7]});
            ncycle = new ArrayList<Integer>();
            Collections.addAll(ncycle, new Integer[]{pol[3],pol[6],pol[7],pol[1]});
            internal.put(pol[5],ncycle);
            ncycle = new ArrayList<Integer>();
            Collections.addAll(ncycle, new Integer[]{pol[2],pol[6],pol[5]});
            internal.put(pol[3],ncycle);
          }
        }
        ArrayList<Integer> ncycle = new ArrayList<Integer>();
        Collections.addAll(ncycle, innersubdiv2);
        ncycle.add(92);
        ncomplexfamily.add(ncycle);
        for (int i = 0; i < innersubdiv2.length; i++){
          internal.get(innersubdiv2[i]).add(92);
        }
      }

      //maxlabel[0] += 1;
      //newLabels.add(maxlabel[0]);
      //ArrayList<Integer> cycle = new ArrayList<Integer>(ComplexFamily.get(i));
      //internal.put(maxlabel[0], cycle);
      //for (int j = 0; j < cycle.size(); j++){
      //  ArrayList<Integer> ncycle = new ArrayList<Integer>();
      //  ncycle.add(maxlabel[0]);
      //  ncycle.add(cycle.get(j));
      //  ncycle.add(cycle.get((j+1)%cycle.size()));
      //  ncomplexfamily.add(ncycle);
      //  Integer v = cycle.get(j);
      //  if (internal.containsKey(v)){
      //    Integer ni1 = j-1;
      //    if (j == 0){
      //       ni1 = cycle.size()-1;
      //    }
      //    Integer ni2 = (j+1)%cycle.size();
      //    Integer nlabel1 = cycle.get(ni1);
      //    Integer nlabel2 = cycle.get(ni2);
      //    ArrayList<Integer> vpetal = internal.get(v);
      //    //println(nlabel1);
      //    //println(nlabel2);
      //    //println(v);
      //    //println(vpetal);
      //    Integer vpind1 = vpetal.indexOf(nlabel1);
      //    Integer vpind2 = vpetal.indexOf(nlabel2);
      //    boolean t1 = vpind1 == 0 || vpind1 == (vpetal.size()-1);
      //    boolean t2 = vpind2 == 0 || vpind2 == (vpetal.size()-1);
      //    if (t1 && t2){
      //      vpetal.add(maxlabel[0]);
      //    }
      //    else{
      //      Integer insIndex = max(vpind1,vpind2);
      //      vpetal.add(insIndex, maxlabel[0]);
      //    }
      //  }
      //}
    }
    ComplexFamilyPacks.add(ncomplexfamily);
}

void triangulateComplex(HashMap<Integer, ArrayList<Integer>> internal, 
                        Integer[] maxlabel){
    Integer tlabelSize = ComplexFamily.size()+1;
    ArrayList<Integer> newLabels = new ArrayList<Integer>();
    ArrayList<ArrayList<Integer>> ncomplexfamily = new ArrayList<ArrayList<Integer>>();
    
    for (int i = 0; i < ComplexFamily.size(); i++){
      maxlabel[0] += 1;
      newLabels.add(maxlabel[0]);
      ArrayList<Integer> cycle = new ArrayList<Integer>(ComplexFamily.get(i));
      internal.put(maxlabel[0], cycle);
      for (int j = 0; j < cycle.size(); j++){
        ArrayList<Integer> ncycle = new ArrayList<Integer>();
        ncycle.add(maxlabel[0]);
        ncycle.add(cycle.get(j));
        ncycle.add(cycle.get((j+1)%cycle.size()));
        ncomplexfamily.add(ncycle);
        Integer v = cycle.get(j);
        if (internal.containsKey(v)){
          Integer ni1 = j-1;
          if (j == 0){
             ni1 = cycle.size()-1;
          }
          Integer ni2 = (j+1)%cycle.size();
          Integer nlabel1 = cycle.get(ni1);
          Integer nlabel2 = cycle.get(ni2);
          ArrayList<Integer> vpetal = internal.get(v);
          //println(nlabel1);
          //println(nlabel2);
          //println(v);
          //println(vpetal);
          Integer vpind1 = vpetal.indexOf(nlabel1);
          Integer vpind2 = vpetal.indexOf(nlabel2);
          boolean t1 = vpind1 == 0 || vpind1 == (vpetal.size()-1);
          boolean t2 = vpind2 == 0 || vpind2 == (vpetal.size()-1);
          if (t1 && t2){
            vpetal.add(maxlabel[0]);
          }
          else{
            Integer insIndex = max(vpind1,vpind2);
            vpetal.add(insIndex, maxlabel[0]);
          }
        }
      }
    }
    ComplexFamily = new ArrayList<ArrayList<Integer>>(ncomplexfamily);
}

HashMap<Integer, ArrayList<Integer>> internal;

void setup(){
  size(1080,720);
  internal = new HashMap<Integer, ArrayList<Integer>>();
  ArrayList<Integer> ivgraph = new ArrayList<Integer>();
  Integer[] ivgarr = {1,2,3,4,5};
  Collections.addAll(ivgraph,ivgarr);
  println(ivgraph);
  internal.put(6, ivgraph);
  HashMap<Integer, Float> external = new HashMap<Integer,Float>();
  external.put(1,Eradius);
  external.put(2,Eradius);
  external.put(3,Eradius);
  external.put(4,Eradius);
  external.put(5,Eradius);
  ArrayList<ArrayList<Integer>> complexfamily = new ArrayList<ArrayList<Integer>>();
  ComplexFamily.add(ivgraph);
  Integer[] maxlabel = {6};
  //ComplexFamily = subdivide(ComplexFamily, internal, external, maxlabel);
  pass += 1;
  //ComplexFamily = subdivide(ComplexFamily, internal, external, maxlabel);
  //ComplexFamily = subdivide(ComplexFamily, internal, external, maxlabel);
  triangulateComplex2(internal, external, maxlabel);
  println(internal);
  println(external);
  println(ComplexFamily);
  HashMap<Integer, Integer> vcount = new HashMap<Integer,Integer>();
  for(ArrayList<Integer> cycle: ComplexFamily){
    for(Integer v: cycle){
      if (vcount.containsKey(v)){
        vcount.put(v, vcount.get(v)+1);
      }
      else{
        vcount.put(v, 1);
      }
    }
  }
  println(vcount);
  //triangulateComplex(internal, maxlabel);
  
  //triangulateComplex(internal, maxlabel);
  //triangulateComplex(internal, maxlabel);
  println(internal);
  circlePack(internal, external, out);
}

void draw(){
  background(0);
  translate(1080.0/2.0, 720.0/2.0);
  for (Map.Entry<Integer,Circle> circ : out.entrySet()){
    Circle circv = circ.getValue();
    noFill();
    stroke(255);
    ellipse(circv.center.x, circv.center.y, 1.3*2*circv.radius, 1.3*2*circv.radius);
  }
  for (Map.Entry<Integer,ArrayList<Integer>> intern : internal.entrySet()){
    for (Integer av : intern.getValue()){
      line(out.get(intern.getKey()).center.x, out.get(intern.getKey()).center.y,
           out.get(av).center.x, out.get(av).center.y);
    }
  }
  //saveFrame();
}