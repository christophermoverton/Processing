import java.util.Map;
import java.util.*;

/*
Compute circle packings according to the Koebe-Thurston-Andreev theory,
Following a numerical algorithm by C. R. Collins and K. Stephenson,
"A Circle Packing Algorithm", Comp. Geom. Theory and Appl. 2003.
*/
float Eradius = 60.0;
float tolerance  = 1.0+1.0e-12;
Integer iterationMax = 3000;
ArrayList<ArrayList<Integer>> ComplexFamily = new ArrayList<ArrayList<Integer>>();

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
      Integer nlabel = complex.get((i+1)%complex.size());
      Integer ilabel = complex.get(i);
      Integer[] lpair = {nlabel,ilabel};
      Integer[] lpair2 = {ilabel, nlabel};
      ArrayList<Integer> lpair2al = new ArrayList<Integer>();
      Collections.addAll(lpair2al, lpair2);
      ArrayList<Integer> lpairal = new ArrayList<Integer>();
      Collections.addAll(lpairal, lpair);
      if (subdivEdges2.containsKey(lpairal)){
        subdivlabels[i] = subdivEdges2.get(lpairal);

      }
      else if (subdivEdges2.containsKey(lpair2al)){
        subdivlabels[i] = subdivEdges2.get(lpair2al);
      }
      else{
        maxlabel[0] += 1;
        subdivlabels[i] = maxlabel[0];
        subdivEdges.put(lpair2, maxlabel[0]);
        subdivEdges2.put(lpair2al, maxlabel[0]);
      }
      if (internal.containsKey(nlabel)){
        ArrayList<Integer> cycle = internal.get(nlabel);
        Integer ilcind = cycle.indexOf(ilabel);
        if (ilcind != -1){
          cycle.set(ilcind, subdivlabels[i]);
        }
      }
      if (internal.containsKey(ilabel)){
        ArrayList<Integer> cycle = internal.get(ilabel);
        Integer nlcind = cycle.indexOf(nlabel);
        if (nlcind != -1){
          cycle.set(nlcind, subdivlabels[i]);
        }
      }

      if (internal.containsKey(subdivlabels[i])){
        internal.get(subdivlabels[i]).add(isubdivlabels[i]);
      }
      else{
        if (internal.containsKey(ilabel)&&internal.containsKey(nlabel)){
          //set internal
          ArrayList<Integer> cycle = new ArrayList<Integer>();
          cycle.add(nlabel);
          cycle.add(isubdivlabels[i]);
          cycle.add(ilabel);
          internal.put(subdivlabels[i],cycle);
        }
        else{
          external.put(subdivlabels[i], Eradius);
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
      Integer label3 = subdivlabels[i];
      ArrayList<Integer> cycle = new ArrayList<Integer>();
      Collections.addAll(cycle, new Integer[] {label1,label2,label3});
      internal.put(isubdivlabels[i], cycle);
      label1 = isubdivlabels[i];
      label2 = subdivlabels[i];
      label3 = complex.get(i);
      Integer label4 = subdivlabels[(i+1)%subdivlabels.length];
      Integer label5 = isubdivlabels[(i+1)%isubdivlabels.length];
      ArrayList<Integer> ncomplex = new ArrayList<Integer>();
      Collections.addAll(ncomplex, new Integer[] {label1,label2,label3,label4,label5});
      inner5.add(isubdivlabels[i]);
      nout.add(ncomplex);
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

void triangulateComplex(ArrayList<ArrayList<Integer>> complexfamily,
                        HashMap<Integer, ArrayList<Integer>> internal, 
                        Integer[] maxlabel){
    Integer tlabelSize = complexfamily.size()+1;
    ArrayList<Integer> newLabels = new ArrayList<Integer>();
    for (int i = 0; i < complexfamily.size(); i++){
      maxlabel[0] += 1;
      newLabels.add(maxlabel[0]);
      ArrayList<Integer> cycle = new ArrayList<Integer>(complexfamily.get(i));
      internal.put(maxlabel[0], cycle);
      for (int j = 0; j < cycle.size(); j++){
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
          println(nlabel1);
          println(nlabel2);
          println(v);
          println(vpetal);
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
}

void setup(){
  size(1080,720);
  HashMap<Integer, ArrayList<Integer>> internal = new HashMap<Integer, ArrayList<Integer>>();
  ArrayList<Integer> ivgraph = new ArrayList<Integer>();
  Integer[] ivgarr = {1,2,3,4,5};
  Collections.addAll(ivgraph,ivgarr);
  println(ivgraph);
  //internal.put(6, ivgraph);
  HashMap<Integer, Float> external = new HashMap<Integer,Float>();
  external.put(1,Eradius);
  external.put(2,Eradius);
  external.put(3,Eradius);
  external.put(4,Eradius);
  external.put(5,Eradius);
  ArrayList<ArrayList<Integer>> complexfamily = new ArrayList<ArrayList<Integer>>();
  complexfamily.add(ivgraph);
  Integer[] maxlabel = {5};
  complexfamily = subdivide(complexfamily, internal, external, maxlabel);
  complexfamily = subdivide(complexfamily, internal, external, maxlabel);
  println(internal);
  println(external);
  println(complexfamily);
  triangulateComplex(complexfamily,internal, maxlabel);
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
    ellipse(circv.center.x, circv.center.y, 2*circv.radius, 2*circv.radius);
  }
}