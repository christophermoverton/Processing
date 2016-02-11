import java.util.Collections;
void ptSetinPolygon(ArrayList<PVector> vertices, ArrayList<PVector> pts, ArrayList<PVector> out){
  //point in convex polygon testing
  for (int i = 0; i < pts.size(); i++){
    //get min max y for x test y interior for pt in pts
    boolean test1 = false;
    ArrayList<ArrayList<PVector>> edges = new ArrayList<ArrayList<PVector>>();
    for (int j = 0; j < vertices.size(); j++){
      test1 = false;
      ArrayList<PVector> edge = new ArrayList<PVector>();
      float yi = vertices.get(j).y;
      float nyi = vertices.get((j+1)% vertices.size()).y;
      float py = pts.get(i).y;
      boolean e1 = yi <= py;
      boolean e2 = py <= nyi;
      boolean e3 = yi >= py;
      boolean e4 = py >= nyi;
      if (e1 && e2) {
        test1 = true;
        edge.add(vertices.get(j));
        edge.add(vertices.get((j+1)% vertices.size()));
        edges.add(edge);
      }
      if (e3 && e4 && !test1){
  
        edge.add(vertices.get(j));
        edge.add(vertices.get((j+1)% vertices.size()));
        edges.add(edge);
      }
    }
    if (edges.size() == 2){
      ArrayList<PVector> edg1 = edges.get(0);
      ArrayList<PVector> edg2 = edges.get(1);
      float x1 = xfromLine(edg1.get(0), edg1.get(1), pts.get(i).y);
      float x2 = xfromLine(edg2.get(0), edg2.get(1), pts.get(i).y);
      boolean e1 = min(x1,x2) <= pts.get(i).x;
      boolean e2 = pts.get(i).x <= max(x1,x2); 
      if ( e1 && e2){
        out.add(pts.get(i));
      }
    }
  }
}

float xfromLine(PVector p1, PVector p2, float y){
  return (p2.x-p1.x)*(y-p1.y)/(p2.y-p1.y)+p1.x;
}

void createNGon(ArrayList<PVector> pos, PShape s){
  s.beginShape();
  for (int i = 0; i < pos.size(); i++){
    s.vertex(pos.get(i).x, pos.get(i).y);
  }
  s.stroke(255);
  s.strokeWeight(.5);
  s.noFill();
  s.endShape(CLOSE);
}

ArrayList<PVector> verts = new ArrayList<PVector>();
PShape s1;
void setup(){
  size(1280,720);
  background(0);
  PVector[] vertices = {new PVector(0.0,200.0,0.0), new PVector(200.0,0.0,0.0),
                        new PVector(-200.0,0.0,0.0)};
  Collections.addAll(verts,vertices);
  s1 = createShape();
  createNGon(verts, s1);
}
void draw(){
  fill(0,2);
  float angle = random(PI*2.0);
  float radius = random(300.0);
  PVector t = new PVector(radius, 0.0,0.0);
  t.rotate(angle);
  //println(t);
  ArrayList<PVector> pts = new ArrayList<PVector>();
  pts.add(t);
  ArrayList<PVector> out = new ArrayList<PVector>();
  ptSetinPolygon(verts, pts, out);
  translate(1280.0/2.0, 720.0/2.0);
  stroke(255);
  strokeWeight(.5);
  noFill();
  triangle(verts.get(0).x, verts.get(0).y, verts.get(1).x, verts.get(1).y, verts.get(2).x,
           verts.get(2).y);
  //shape(s1,0.0,0.0);

  if (out.size() != 0){
    fill(255,90);
    ellipse(t.x,t.y, 10.0,10.0);
    //println(t);
  }
}