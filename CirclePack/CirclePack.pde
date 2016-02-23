import java.util.Map;
import java.util.*;

/*
Compute circle packings according to the Koebe-Thurston-Andreev theory,
Following a numerical algorithm by C. R. Collins and K. Stephenson,
"A Circle Packing Algorithm", Comp. Geom. Theory and Appl. 2003.
*/
float tolerance  = 1.0+1.0e-12;

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
        return acos(( pow((x+y),2)+pow((x+z),2)-pow((y+z),2))/(2.0*(x+y)*(x+z)));
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
      sum += acxyz(radius.get(center), radius.get(j), radius.get(i));
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
    for (int i = -1 * cycle.size()+1; i < cycle.size(); i++){
    //for i in range(-len(cycle),len(cycle)-1):
        //if cycle[i] in placements and cycle[i+1] not in placements:
        int ii = i;
        if (i < 0){
          ii = i%cycle.size();
          ii += cycle.size();
        }
        int j = (ii+1)%cycle.size();
        if (placements.containsKey(cycle.get(ii)) && !placements.containsKey(cycle.get(j))){
            Integer s = cycle.get(i);
            Integer t = cycle.get(j);
            float theta = acxyz(radii.get(center),radii.get(s),radii.get(t));
            PVector offset = new PVector(0.0,0.0,0.0);
            offset.x = (placements.get(s).x - placements.get(center).x)/(radii.get(s)+ radii.get(center));
            offset.y = (placements.get(s).y - placements.get(center).y)/(radii.get(s)+ radii.get(center));
            //offset = (placements[s]-placements[center])/(radii[s]+radii[center]);
            //offset *= e**(-1j*theta);
            offset.x *= cos(theta);
            offset.y *= -1.0*sin(theta);
            PVector tvec = new PVector(0.0,0.0,0.0);
            tvec.x = placements.get(center).x + offset.x*(radii.get(t)+radii.get(center));
            tvec.y = placements.get(center).y + offset.y*(radii.get(t)+radii.get(center));
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
    if (external.containsKey(k)){
      println ("CirclePack internal and external keys need to be disjoint.");
    }
    radii.put(k.getKey(), 1.0);
  }
  float lastChange = 2.0;
  while (lastChange > tolerance){
    lastChange = 1;
    for (Map.Entry<Integer,ArrayList<Integer>> k : internal.entrySet()){
      float theta = flower(radii, k.getKey(), internal.get(k.getValue()));
      float hat = radii.get(k.getKey())/(1/sin(theta/(2*(internal.get(k.getValue())).size()))-1);
      float newrad = hat * (1/(sin(PI/(internal.get(k.getValue())).size())) - 1);
      float kc = max(newrad/radii.get(k.getKey()),radii.get(k.getKey())/newrad);
      lastChange = max(lastChange,kc);
      radii.put(k.getKey(),newrad);
    }
  }
  
  HashMap<Integer,PVector> placements = new HashMap<Integer,PVector>();
  Iterator it = internal.entrySet().iterator();
  Map.Entry<Integer,ArrayList<Integer>> k1 = (Map.Entry)it.next();
  placements.put(k1.getKey(),new PVector(0.0,0.0,0.0));
  Integer k2 = k1.getValue().get(0);
  placements.put(k2, new PVector(radii.get(k1.getKey())+radii.get(k2), 0.0, 0.0));
  place(placements, radii, internal, k1.getKey());
  place(placements, radii, internal, k2);
  for (Map.Entry<Integer, Float> k : radii.entrySet()){
    out.put(k.getKey(), new Circle(placements.get(k.getKey()), k.getValue()));
  }
}