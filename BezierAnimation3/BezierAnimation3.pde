import java.util.Map;
PVector[] points = {new PVector(0.0, 0.0, 0.0), new PVector(0.6294293999671936, 0.7892796993255615, 0.0), new PVector(-0.38169065117836, 1.6722959280014038, 0.0), new PVector(-1.9889508485794067, 0.9578285217285156, 0.0), new PVector(-2.65193510055542, -1.2771034240722656, 0.0), new PVector(-0.9542266130447388, -4.180739879608154, 0.0), new PVector(3.776573657989502, -4.735679626464844, 0.0), new PVector(7.539822578430176, -1.7976337858272018e-06, 0.0), new PVector(5.035435199737549, 6.314237594604492, 0.0), new PVector(-1.717611312866211, 7.525330066680908, 0.0), new PVector(-6.629836082458496, 3.1927616596221924, 0.0), new PVector(-7.292819499969482, -3.512037515640259, 0.0), new PVector(-2.2901484966278076, -10.033774375915527, 0.0), new PVector(8.18258285522461, -10.260636329650879, 0.0), new PVector(15.079645156860352, 3.5952675716544036e-06, 0.0), new PVector(9.441457748413086, 11.839181900024414, 0.0), new PVector(-3.053513288497925, 13.37837028503418, 0.0), new PVector(-11.270718574523926, 5.427700042724609, 0.0), new PVector(-11.93370532989502, -5.746971130371094, 0.0), new PVector(-3.6260826587677, -15.88680648803711, 0.0), new PVector(12.588571548461914, -15.785606384277344, 0.0), new PVector(22.61946678161621, -1.617870293557644e-05, 0.0), new PVector(13.847456932067871, 17.364147186279297, 0.0), new PVector(-4.389442443847656, 19.231401443481445, 0.0), new PVector(-15.911595344543457, 7.6626505851745605, 0.0), new PVector(-16.57460594177246, -7.98187255859375, 0.0), new PVector(-4.962007522583008, -21.739839553833008, 0.0), new PVector(16.994571685791016, -21.31056785583496, 0.0), new PVector(30.159290313720703, -2.1571606339421123e-05, 0.0), new PVector(18.253509521484375, 22.889068603515625, 0.0), new PVector(-5.725313663482666, 25.084447860717773, 0.0), new PVector(-20.55247688293457, 9.897589683532715, 0.0), new PVector(-21.215496063232422, -10.216796875, 0.0), new PVector(-6.297882080078125, -27.592884063720703, 0.0), new PVector(21.40052032470703, -26.8355712890625, 0.0), new PVector(37.69911193847656, -9.886985935736448e-05, 0.0), new PVector(22.659420013427734, 28.41410255432129, 0.0), new PVector(-7.061220169067383, 30.93748664855957, 0.0), new PVector(-25.193359375, 12.132530212402344, 0.0), new PVector(-25.856359481811523, -12.451769828796387, 0.0), new PVector(-7.633796691894531, -33.445919036865234, 0.0), new PVector(25.80663299560547, -32.360443115234375, 0.0), new PVector(45.23893356323242, 5.3929012210574e-05, 0.0), new PVector(27.06541633605957, 33.93906784057617, 0.0), new PVector(-8.397266387939453, 36.79049301147461, 0.0), new PVector(-29.83429718017578, 14.367356300354004, 0.0), new PVector(-30.49724578857422, -14.6867036819458, 0.0), new PVector(-8.969861030578613, -39.298919677734375, 0.0), new PVector(30.212499618530273, -37.8855094909668, 0.0), new PVector(52.77875518798828, -0.00013841778854839504, 0.0), new PVector(31.471567153930664, 39.46390914916992, 0.0), new PVector(-9.733033180236816, 42.64356231689453, 0.0), new PVector(-34.475120544433594, 16.60240936279297, 0.0), new PVector(-35.13813018798828, -16.9216365814209, 0.0), new PVector(-10.305624961853027, -45.15199279785156, 0.0), new PVector(34.618656158447266, -43.41035079956055, 0.0), new PVector(60.318580627441406, 7.190534961409867e-05, 0.0), new PVector(35.877410888671875, 44.98899459838867, 0.0), new PVector(-11.068939208984375, 48.49660110473633, 0.0), new PVector(-39.116004943847656, 18.83734893798828, 0.0), new PVector(-39.77901840209961, -19.156570434570312, 0.0), new PVector(-11.641539573669434, -51.00502395629883, 0.0), new PVector(39.02466583251953, -48.935302734375, 0.0), new PVector(67.8583984375, 8.08935146778822e-05, 0.0), new PVector(40.283409118652344, 50.51395797729492, 0.0), new PVector(-12.405052185058594, 54.34959411621094, 0.0), new PVector(-43.7569694519043, 21.07212257385254, 0.0), new PVector(-44.41990280151367, -21.391502380371094, 0.0), new PVector(-12.977887153625488, -56.85796356201172, 0.0), new PVector(43.43027114868164, -54.460594177246094, 0.0), new PVector(75.39822387695312, -0.0004853611171711236, 0.0), new PVector(44.6898307800293, 56.03858184814453, 0.0), new PVector(-13.740981101989746, 60.2026252746582, 0.0), new PVector(-48.397857666015625, 23.307044982910156, 0.0), new PVector(-49.06096649169922, -23.626062393188477, 0.0), new PVector(-14.313846588134766, -62.71099090576172, 0.0), new PVector(47.836238861083984, -59.98557662963867, 0.0), new PVector(82.93804931640625, -0.0005338972550816834, 0.0), new PVector(49.09587478637695, 61.563514709472656, 0.0), new PVector(-15.07640266418457, 66.05577850341797, 0.0), new PVector(-53.03855514526367, 25.54237174987793, 0.0), new PVector(-53.7016716003418, -25.86136817932129, 0.0), new PVector(-15.649282455444336, -68.56413269042969, 0.0), new PVector(52.242698669433594, -65.51016235351562, 0.0), new PVector(90.47786712646484, 0.000107858024421148, 0.0), new PVector(53.50140380859375, 67.08885192871094, 0.0), new PVector(-16.412837982177734, 71.90869140625, 0.0), new PVector(-57.67963790893555, 27.776887893676758, 0.0), new PVector(-58.34255599975586, -28.09630012512207, 0.0), new PVector(-16.985197067260742, -74.41716766357422, 0.0), new PVector(56.648712158203125, -71.03512573242188, 0.0), new PVector(98.01769256591797, 0.00011684619676088914, 0.0), new PVector(57.90740203857422, 72.61381530761719, 0.0), new PVector(-17.74876594543457, 77.76171875, 0.0), new PVector(-62.3203010559082, 30.012285232543945, 0.0), new PVector(-62.983680725097656, -30.33075714111328, 0.0), new PVector(-18.32172393798828, -80.27007293701172, 0.0), new PVector(61.05414581298828, -76.56053924560547, 0.0), new PVector(105.55751037597656, -0.0006795055232942104, 0.0), new PVector(62.313987731933594, 78.13829803466797, 0.0)};

PVector[] lpoints = {new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0)};

PVector[] rpoints = {new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0), new PVector(0.0, 0.0, 0.0)};

class CurveDat{
  PVector pos;
  PVector lcpos;
  PVector rcpos;
  int frameID;
  int pointID;
  int curveID;
  int maxPts;
  CurveDat(){
    pos = new PVector(0.0,0.0,0.0);
    lcpos = new PVector(0.0,0.0,0.0);
    rcpos = new PVector(0.0,0.0,0.0);
    frameID = 0;
    pointID = 0;
    curveID = 0;
    maxPts = 0;
  }
}
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

float tstep = .001;  //should be 0<tstep< 1
int curvelen = 10;//in point steps;
boolean showTrack = false;
boolean showTrack2 = true;
boolean multiCurve = true;
int maxFrames = 0;
HashMap<Integer,HashMap<Integer,PVector[]>> curvetoframepoints;
HashMap<Integer,HashMap<Integer,PVector[]>> curvetoframelpoints;
HashMap<Integer,HashMap<Integer,PVector[]>> curvetoframerpoints;
HashMap<Integer,HashMap<Integer, ArrayList<PVector>>> cfpoints;

ArrayList<PVector> curvepoints = new ArrayList<PVector>();
void setup(){
  size(1080,720);
  background(255);
  smooth(8);
  Table table = loadTable("bezierpoints.csv","header");
  points = new PVector[table.getRowCount()];
  lpoints = new PVector[table.getRowCount()];
  rpoints = new PVector[table.getRowCount()];
  int h = 0;
  HashMap<Integer,ArrayList<Integer>> frametoCDats = new HashMap<Integer,ArrayList<Integer>>();
  ArrayList<CurveDat> CurveDats = new ArrayList<CurveDat>();
  HashMap<Integer,ArrayList<Integer>> curvetoCDats = new HashMap<Integer,ArrayList<Integer>>();
  for (TableRow row : table.rows()) {
    CurveDat curvdat = new CurveDat();
    float x = row.getFloat("x");
    float y = row.getFloat("y");
    float z = row.getFloat("z");
    float lpx = row.getFloat("cplx");
    float lpy = row.getFloat("cply");
    float lpz = row.getFloat("cplz");
    float rpx = row.getFloat("cprx");
    float rpy = row.getFloat("cpry");
    float rpz = row.getFloat("cprz");
    Integer curveID = row.getInt("curveID");
    Integer pointID = row.getInt("pointID");
    Integer frameID = row.getInt("frameID");
    Integer maxPts = row.getInt("maxPts");
    curvdat.pos = new PVector(x,y,z);
    curvdat.lcpos = new PVector(lpx,lpy,lpz);
    curvdat.rcpos = new PVector(rpx,rpy,rpz);
    curvdat.curveID = curveID;
    curvdat.frameID = frameID;
    curvdat.pointID = pointID;
    curvdat.maxPts = maxPts;
    if (curvetoCDats.containsKey(curveID)){
      curvetoCDats.get(curveID).add(CurveDats.size());
    }
    else{
      ArrayList<Integer> cdatslist = new ArrayList<Integer>();
      cdatslist.add(CurveDats.size());
      curvetoCDats.put(curveID, cdatslist);
    }
    if (frametoCDats.containsKey(frameID)){
      frametoCDats.get(frameID).add(CurveDats.size());
    }
    else{
      ArrayList<Integer> cdatslist = new ArrayList<Integer>();
      cdatslist.add(CurveDats.size());
      frametoCDats.put(frameID, cdatslist);
    }
    CurveDats.add(curvdat);
    points[h] = new PVector(x,y,z);
    lpoints[h] = new PVector(lpx,lpy,lpz);
    rpoints[h] = new PVector(rpx,rpy,rpz);
    h += 1;
  }
  maxFrames = frametoCDats.size();
  HashMap<Integer,PVector[]> framepoints = new HashMap<Integer,PVector[]>();
  HashMap<Integer,PVector[]> framelpoints = new HashMap<Integer,PVector[]>();
  HashMap<Integer,PVector[]> framerpoints = new HashMap<Integer,PVector[]>();
  curvetoframepoints = new HashMap<Integer,HashMap<Integer,PVector[]>>();
  curvetoframelpoints = new HashMap<Integer,HashMap<Integer,PVector[]>>();
  curvetoframerpoints = new HashMap<Integer,HashMap<Integer,PVector[]>>();
  for (Map.Entry<Integer,ArrayList<Integer>> me : frametoCDats.entrySet()) {
    Integer frameID = me.getKey();
    ArrayList<Integer> cDats = me.getValue();
    for (Integer cDat : cDats){
      Integer curveID = CurveDats.get(cDat).curveID;
      Integer pointID = CurveDats.get(cDat).pointID;
      PVector pos = CurveDats.get(cDat).pos;
      PVector lcpos = CurveDats.get(cDat).lcpos;
      PVector rcpos = CurveDats.get(cDat).rcpos;
      Integer maxPts = CurveDats.get(cDat).maxPts;
      if (curvetoframepoints.containsKey(curveID)){
        
        framepoints = curvetoframepoints.get(curveID);
        if (framepoints.containsKey(frameID)){
          PVector[] points = framepoints.get(frameID);
          points[pointID] = pos;
        }
        else{
          PVector[] points = new PVector[maxPts];
          points[pointID] = pos;
          framepoints.put(frameID, points);
        }
      }
      else{
         framepoints = new HashMap<Integer,PVector[]>();
         PVector[] points = new PVector[maxPts];
         points[pointID] = pos;
         framepoints.put(frameID, points);
         curvetoframepoints.put(curveID, framepoints);
      }
      if (curvetoframelpoints.containsKey(curveID)){
        
        framelpoints = curvetoframelpoints.get(curveID);
        if (framelpoints.containsKey(frameID)){
          PVector[] lpoints = framelpoints.get(frameID);
          lpoints[pointID] = lcpos;
        }
        else{
          PVector[] lpoints = new PVector[maxPts];
          lpoints[pointID] = lcpos;
          framelpoints.put(frameID, lpoints);
        }
      }
      else{
         framelpoints = new HashMap<Integer,PVector[]>();
         PVector[] lpoints = new PVector[maxPts];
         lpoints[pointID] = lcpos;
         framelpoints.put(frameID, lpoints);
         curvetoframelpoints.put(curveID, framelpoints);
      }
      if (curvetoframerpoints.containsKey(curveID)){
        
        framerpoints = curvetoframerpoints.get(curveID);
        if (framerpoints.containsKey(frameID)){
          PVector[] rpoints = framerpoints.get(frameID);
          rpoints[pointID] = rcpos;
        }
        else{
          PVector[] rpoints = new PVector[maxPts];
          rpoints[pointID] = rcpos;
          framerpoints.put(frameID, rpoints);
        }
      }
      else{
         framerpoints = new HashMap<Integer,PVector[]>();
         PVector[] rpoints = new PVector[maxPts];
         rpoints[pointID] = rcpos;
         framerpoints.put(frameID, rpoints);
         curvetoframerpoints.put(curveID, framerpoints);
      }
    }
  }
  cfpoints = new HashMap<Integer,HashMap<Integer, ArrayList<PVector>>>();
  int steps = (int)(1/tstep);
  for (Map.Entry<Integer,HashMap<Integer,PVector[]>> me : curvetoframepoints.entrySet()){
    Integer curveID = me.getKey();
    framepoints = me.getValue();
    for (Map.Entry<Integer,PVector[]> me2 : framepoints.entrySet()){
      Integer frameID = me2.getKey();
      PVector[] mappoints = me2.getValue();
      PVector[] maplpoints = curvetoframelpoints.get(curveID).get(frameID);
      PVector[] maprpoints = curvetoframerpoints.get(curveID).get(frameID);
      ArrayList<PVector> newpoints = new ArrayList<PVector>();
      if (cfpoints.containsKey(curveID)){
        HashMap<Integer,ArrayList<PVector>> fpoints = cfpoints.get(curveID);
        if (fpoints.containsKey(frameID)){
          newpoints = fpoints.get(frameID);
        }
        else{
          fpoints.put(frameID, newpoints);
        }
      }
      else{
        HashMap<Integer, ArrayList<PVector>> fpoints = new HashMap<Integer,ArrayList<PVector>>();
        fpoints.put(frameID, newpoints);
        cfpoints.put(curveID,fpoints);
      }
      for (int i = 1; i < mappoints.length; i++){
        int j = (i-1)%mappoints.length;
        PVector P1 = mappoints[j];
        PVector P4 = mappoints[i];
        PVector P3 = maplpoints[i];
        PVector P2 = maprpoints[j];  

        float pstep = 0.0;
        for (int k = 0; k < steps; k++){
            PVector out = new PVector(0.0,0.0,0.0);
            CubicBezier(pstep, P1, P2, P3, P4, out);
            newpoints.add(out);
            pstep += tstep;
        }
      }
      PVector P1 = mappoints[mappoints.length-1];
      PVector P4 = mappoints[0];
      PVector P3 = maplpoints[0];
      PVector P2 = maprpoints[mappoints.length-1];
      float pstep = 0.0;
      for (int k = 0; k < steps; k++){
        PVector out = new PVector(0.0,0.0,0.0);
        CubicBezier(pstep, P1, P2, P3, P4, out);
        newpoints.add(out);
        pstep += tstep;
      }        
    }
    
  }
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
  //println(curvepoints.size());
  println(curvetoframepoints);
}
float astep = 0;
int astep2 = 0;
float astep3 = 0;

void draw(){
  fill(255,0);
  rect(0,0, 1920,1080);
  if (astep > (curvepoints.size()-curvelen-1)){
    astep = 0;
    println("hit");
  }
  translate(1080.0/2.0, 720.0/2.0);
  //rotate(.1*pow(astep2, .5));
  scale(50.0);
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
  if (showTrack2){
    Integer frameID = ((int) astep3)%maxFrames + 1;
    for (Map.Entry<Integer,HashMap<Integer,PVector[]>> me : curvetoframepoints.entrySet()){
      Integer curveID = me.getKey();
      HashMap<Integer, PVector[]> framepoints = me.getValue();
      
      //println(frameID);
      //for (Map.Entry<Integer,PVector[]> me2 : framepoints.entrySet()){
        //Integer frameID = me2.getKey();
        //PVector[] mappoints = me2.getValue();
        PVector[] mappoints = curvetoframepoints.get(curveID).get(frameID);
        PVector[] maplpoints = curvetoframelpoints.get(curveID).get(frameID);
        PVector[] maprpoints = curvetoframerpoints.get(curveID).get(frameID);
        print(mappoints);
        beginShape();
        vertex(mappoints[0].x,mappoints[0].y);
        for(int i = 1; i < mappoints.length; i++){
          int j = (i-1)%mappoints.length;
          i = i % mappoints.length;
          float x1 = mappoints[i].x;
          float y1 = mappoints[i].y;
          float x2 = maplpoints[i].x;
          float y2 = maplpoints[i].y;
          float x3 = maprpoints[j].x;
          float y3 = maprpoints[j].y;
          float x4 = mappoints[j].x;
          float y4 = mappoints[j].y;
          bezierVertex(x3, y3, x2, y2, x1, y1);
          //bezier(x1, y1, x2, y2, x3, y3, x4, y4);
        }
        endShape();
        PVector lp1 = mappoints[mappoints.length-1];
        PVector lp2 = mappoints[0];
        PVector lpc1 = maplpoints[0];
        PVector lpc2 = maprpoints[mappoints.length-1];
        bezier(lp2.x,lp2.y, lpc1.x, lpc1.y, lpc2.x, lpc2.y, lp1.x, lp1.y);
      //}
    }
  }
  float astepi = astep;
  stroke(0);
  strokeWeight(.01);
  beginShape();
  if (multiCurve){
    Integer frameID = ((int) astep3)%maxFrames + 1;
    for (Map.Entry<Integer,HashMap<Integer,ArrayList<PVector>>> me : cfpoints.entrySet()){
      Integer curveID = me.getKey();
      Integer cfsize = cfpoints.get(curveID).get(frameID).size();
      int k = 0;
      while(k < (curvelen)){
        Integer poskey = ((int)astepi)%cfsize;
        curveVertex(cfpoints.get(curveID).get(frameID).get(poskey).x, 
                    cfpoints.get(curveID).get(frameID).get(poskey).y);
        astepi += .5;
        k += 1;
      }
      endShape();   
    }
  }
  else{
    //int k = 0;
    //while(k < (curvelen)){
    //  curveVertex(curvepoints.get(astepi).x, curvepoints.get(astepi).y);
    //  astepi += 1;
    //  k += 1;
    //}
    //endShape();
  }
  astep += 100.0;
  astep2 += 1000;
  astep3 += 1;

}