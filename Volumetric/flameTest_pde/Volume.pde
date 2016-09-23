interface Volume
{
  vRes sample( PVector p );
  float getRange();
  int numSamples();
  public void updateCount();
  public float getCount();
}
class Sphere implements Volume
{
  float count = .01f;
  vRes  sample( PVector p )
  {
    float a = smoothstep( 1.,0.9,p.mag())*0.1;
    return new vRes( new PVector(4.,4.,4.), a);
  }
  void updateCount(){
    count += .01f;
  }
  float getCount(){
    return count;
  }
  float getRange() { return 1.;}
  int numSamples() { return 20;}
}
 
class Cloud implements Volume
{
  float count = .01f;
  void updateCount(){
    count += .01;
  }
  float getCount(){
    return count;
  }
  vRes  sample( PVector p )
  {
    noiseDetail(6,0.6);
    float d = p.mag();
     float nf = 1.f;
     float namp = 3.f;
    d += (noise(p.x*nf+57.f,p.y*nf+513.f,p.z*nf+13.f+count)-.5)*namp;
     
    float a = smoothstep( 1.,0.9, d)*0.3;
     
    // do approx lighting on p.y
    float lgt = constrain( (-p.y * d ) *.5+.5,0.5,1.);
    return new vRes( new PVector(4.*lgt,4.*lgt,4.*lgt), a);
  }
  float getRange() { return 1.0;}
  int numSamples() { return 80;}
}
// does planet amtosphere
class EdgedSphere implements Volume
{
  float count = .01f;
  void updateCount(){
    count += count;
  }
  float getCount(){
    return count;
  }
  vRes  sample( PVector p )
  {
    float d =p.mag();
     
    float s = 1.;
    float e = 0.9;
    if ( d < e )
    {
      return  new vRes( new PVector(0.8,0.4,0.1),1.);
    }
    float a = constrain((d - s)/(e-s),0.,1);
    a*=a;
    float den = 0.8;
    return new vRes( new PVector(1.,1.,4.), a*den);
  }
  float getRange() { return 1.;}
  int numSamples() { return 30;}
}
// does planet amtosphere
class CrabNebula implements Volume
{
  float count = .01f;
  void updateCount(){
    count += .01;
  }
  float getCount(){
    return count;
  }
  vRes  sample( PVector p )
  {
    noiseDetail(6,0.6);
    float d =p.mag();
    d *=2.;
     float nf = 1.f;
     float namp = 6.f;
    d += (abs(noise(p.x*nf+57.f,p.y*nf+513.f+.05f*count,p.z*nf+13.f+count)-.5)-.25)*namp;
    d = max(d,0.);
    float s = 1.;
    float e = 0.9;
   float a = constrain((d - s)/(e-s),0.,1);
    
   float a2 = constrain(1. -d *10.,0.,1);
    
   float a3 = constrain(1.-abs((d-0.3)*10.),0.,1);
   a3 = a3*5.;
    
   float c = a;
    a*=a;
    float da = max(a,a2);
    da= max(a3,da);
    float den = 0.1;
    PVector cExterior = new PVector(4.+a2*2.,4.+a2*2.,8.);
    PVector cInterior = new PVector(4.,1.,1.);
    PVector cLine = new PVector(4.,4.,1.);
    PVector cr = a > a2 ? cInterior : cExterior;
    cr = a3 > 0 ? cLine : cr;
     
    return new vRes( cr, da*den);
  }
  float getRange() { return 1.5;}
  int numSamples() { return 40;}
}
class Nebula2 implements Volume
{
  float count = .01f;
  void updateCount(){
    count += .01f;
  }
  float getCount(){
    return count;
  }
  vRes  sample( PVector p )
  {
    noiseDetail(6,0.6);
    float d =p.mag();
    d *=2.;
     float nf = 1.f;
     float namp = 6.f;
    d += (abs(noise(p.x*nf+357.f,p.y*nf+53.f,p.z*nf+13.f+count)-.5)-.25)*namp;
    d = max(d,0.);
   float a3 = 1. - min(abs(d-0.5)*2.,1.);
   float ba = a3*a3;
   ba = ba*ba;
   a3 = ba + a3*.5f;
    
    PVector cLine = new PVector(2.+ba*2.,2.+ba*2.,4.);   
    return new vRes( cLine, a3*0.3);
  }
  float getRange() { return 1.5;}
  int numSamples() { return 80;}
}
 
class Eskimo implements Volume
{
  float count = .01f;
  vRes  sample( PVector p )
  {
   // a hollow sphere with slight wobble
    float d = p.mag();
     
    float c= smoothstep( 0.2f,0.0f,d);
     
     noiseDetail(4,0.4);
    float nf = 1.f;
     float namp = 1.f;
    d += (abs(noise(p.x*nf+57.f+count,p.y*nf+513.f,p.z*nf+13.f)-.5)-.25)*namp;
 
    float c3= smoothstep( 0.,1.1,d) - smoothstep( 1.1,1.3,d);   
    d = smoothstep( 0.68f,0.7f,d)-smoothstep(0.7f,0.8f,d);
   d +=c*2.;
   PVector col = new PVector(2.+c*8.,2.+c*3.,24.);
   col.mult(d);
   col.x += c3*5.;
   d += c3*0.1;
   return new vRes( col,d *.1);
  }
  void updateCount(){
    count += .01f;
  }   
  float getCount(){
    return count;
  }
  float getRange() { return 1.5;}
  int numSamples() { return 40;}
};
 
float smoothstep (float edge0, float edge1, float x)
{
  x =  min(max((x - edge0) / (edge1 - edge0),0.0f),1.0f);
  return x*x*(3-2*x);
}
class HieghtField implements Volume
{
  float count = .1f;
  void updateCount(){
    count += .1f;
  }
  PVector col( float d )
  {
     return new PVector(1.+d*8.,6.,1.);
  }
  vRes  sample( PVector p )
  {
    float d = abs(p.y);
    float nf = 1.f;
    noiseDetail(5,0.6f);
    d += (abs(noise(p.x*nf+57.f,p.z*nf+13.f+count)-.5)-.25)*1.f;
    float a = smoothstep( 0.45f,0.5,d) - smoothstep(0.5,0.8f,d);
     
    a += smoothstep( 0.2,0.0,abs(p.x)) * smoothstep( 0.2,0.0,abs(p.z)) * 2.;
    
    //print("Hit count");
    return new vRes( col(d), a*.1);
  }
  float getCount(){
    return count;
  }
  int numSamples() { return 30;}
  float getRange() { return 2.0;};
}
class ColouredSphere implements Volume
{
  float count = .0001f;
  float getCount(){
    return count;
  }
  void updateCount(){
    count += .0001f;
  }
  PVector col( float d )
  {
    float id = constrain( 2.0f - d*2., 0.,1.);
    float r = id;
    float b = 1.-abs(id-.5f)*2.;
    float g = id*id*id;
    float i = 2.f;
     return new PVector(r*i,g*i,b*i);
  }
  vRes  sample( PVector p )
  {
    float d = p.mag();
    float nf = 3.f;
    d += (abs(noise(p.x*nf+57.f+count,p.y*nf+513.f,p.z*nf+13.f)-.5)-.25)*1.f;
    float a = smoothstep( 1.f,0.9,d);
    return new vRes( col(d), a);
  }
   
  int numSamples() { return 30;}
  float getRange() { return 1.5;};
}
Volume Volumes[]={
    new HieghtField(),
    new CrabNebula(),
    new Eskimo(),
    new ColouredSphere(),
    new Nebula2(),
    new EdgedSphere(),
    new Cloud(),
    new Sphere(),
};