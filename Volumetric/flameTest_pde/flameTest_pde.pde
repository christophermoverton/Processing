// ftp://69.31.121.43/developer/presentations/2005/GDC/Sponsored_Day/GDC_2005_VolumeRenderingForGames.pdf
int px = 0;
int py = 0;
void setup()
{
  size(640,480);
   
  st = 0;
}
class vRes
{
  PVector c;
  float   a;
  vRes( PVector _c, float _a )
  {
      c = _c;
      a = _a;
  }
}
 
PVector intSphereBound( PVector sp, float r2, PVector ro, PVector dir)
{
  PVector  d = new PVector();
  d.set( ro);  //ray origin
  d.sub( sp);  //sphere origin  subtracting rayorigin with sphere origin
  float c = d.dot(d) - r2*r2;  //distance between the rayotosphereorigin to radius of the sphere
  //hmm...the vector of rayotosphereorigin would eminate technically in the direction from the ray origin towards
  //the sphere origin, so the difference of distance (radius of the sphere) simply 
  //measures the relative squared distance of rayotosphereorigin relative the sphere radius squared.
  //The ray origin I believe should be picked to reside outside of the bounding volume.
  //so that c is always positive.  That is, the distance squared of the ray origin 
  //to the boundary of the surface of the volume.
  float b = d.dot(dir);  //distance of the direction vector in the direction of rayotosphereorigin
  float t = b*b-c;  //square of distance of ray direction vector in the direction of rayotosphereorigin relative to 
  //the distance of the bounding volume's surface.  It should be positive (>0) to count anything for 
  //for a given volumetric measurement.
  if ( t > 0.0f) {  
    float st =sqrt(t);  //distanced difference between surface distance to dir vec distance
    t = -b-st;
    return new PVector( max( -b-st, 0.f), max( -b+st, 0.f),0.f);
    //returns an upper and lower bound of distance (dir) +/- the (distance(dir)^2 - dist(vol_surf)^2)^(1/2) 
  }
  return new PVector( -1.f, -1.f);
}
 
 
void doRay( Raymarcher marcher )
{
  PVector origin = new PVector(0.f,0.f,-4.f);
    PVector dir = new PVector( -width/2,-height/2,width);
    dir.add( new PVector( px, py,0.0));
    dir.normalize();
     
    PVector c =marcher.integrate( origin, dir);
     
    c.mult( 255.);
    // should do gamma
    set( (int)px,(int)py,color(c.x,c.y,c.z));
    px++;
}
int g_CurrentVolume=0;
float oldtime = 0;
float st;
int count = 0;
int count2 = 0;
void mousePressed()
{
  g_CurrentVolume= (g_CurrentVolume+1)%Volumes.length;
  py=0;
  st = 0;
   
}
void draw()
{
   //background(0);
   if ( py == 0)
     background(0);
   
    // lets do a ray
    Raymarcher marcher = new Raymarcher( Volumes[g_CurrentVolume], 0.04f);
   for (int i =0; i < 16; i++)
    {
      px=0;
      for (int j =0; j < width; j++)
         doRay(marcher);
         
       if ( py < height)
         py+=2;
       else if (st == 0)
       {
         st = 1;
         py=1;
       }        
    }
    if (count%80 == 0){
       marcher.updateCount();
       saveFrame("screen_"+String.format("%04d",count2)+".tif");
       background(0);
       
       py = 0;
       st=0;
       count2 += 1;
    }
    count += 1;
    //print("completed cycle");
}