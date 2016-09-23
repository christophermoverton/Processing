class Raymarcher
{
  Volume m_vol;
  float m_stepsize;
   
  Raymarcher( Volume vol, float stepsize )
  {
   m_vol = vol;
    m_stepsize = stepsize;
  }
  PVector integrate( PVector ro, PVector rd )
  {
    float absorption=1.;//2.;
    float spRad = m_vol.getRange();
   PVector tvals = intSphereBound( new PVector(0.f,0.f,0.f),spRad, ro, rd);
    int numsteps = (int)(ceil(tvals.y-tvals.x)/m_stepsize);
    float ds = (tvals.y-tvals.x)/numsteps;
     
    PVector stepdir = new PVector();
    stepdir.set(rd);
    stepdir.mult(ds);
     
    PVector raypos = new PVector();
    raypos.set(rd);
    raypos.mult(tvals.x + random(0,ds));
    raypos.add(ro);
    raypos.add(stepdir);
     
     
    float rhomult = -absorption*ds;
    float T = 1.f;
    PVector ct = new PVector(0.,0.,0.);
     
    for (int step=0; step < numsteps; ++step) {
       vRes r =  m_vol.sample( raypos);
       float rho = r.a;
       T *= exp( rhomult* rho);
       if ( T < 0.0001f)
         break;
       PVector ci = r.c;
       // do lighting
        
       ci.mult( T * ds*rho);
       ct.add(ci);
       raypos.add(stepdir);
    }
    return ct;
  }
};