# Processing

My latest is a flame animation which is actually an animation that attempts that makes use of a particle emitter using
bezier curve data.  Thus controlled curve deflection.  Variables as related to Glow radius are finnicky, so it may require 
tweaking for desired specifications.  At the outset, default exhaust approximation has some translated random motions (but not 
vortical turbulence).  Adding more particle (bezier strands) means greater particle density and hence a more consistent
density appearance, but this also increases the thrust opaqueness which consequently means on Glow radius applying dampening 
factors.  As it turns out increasing the Glow radius dampening factor also means adjusting Glowraddiv, but there is often
because of glow aggregation, very fine thresholds between lack and too much glow (i.e., literally a difference of .1f or less in
some cases even while the factor ranges from 300 to 400 for Glowraddiv).  

minParticleSize and maxParticleSize control the random min max range of the strands size in terms of bezier stroke weight...
glow radius is, however, independent of this value.

MaximumParticles controls the maximum number of particles active at any given time in the particle system in preventing
both undesired computations/ density issues and so forth.  

diffParticleVelocity and minParticleVelocity are deprecated in version 2 of the script but applicable in version 1 of the 
of the script.  Again, this applies in so far as the random assignment of the t step velocity parameter of the bezier 
function.  

turbulence  is a bit of a misnomer in script since this applies more so in terms of translated motions of the exhaust flame,
this doesn't produce voritical turbulence for instance.  

glowdampening factor scales alongside the ith iteration of strokeWeight or in other words strokeWeight(i*glowdampening) on 
the bezier curve glow curve.  

minstrength,particleStrength are applicable to the particle's main stroke alpha channel value in so far as random 
min,max possible assignments.  This is independent of particle glow.  

randomP3 is a boolean value controlling whether or not a supplied bezier curve control point system is parametrically drawn 
from 0 to 1 or to some value less than 1.  
minRandomt3 forces the random assignment of the third control point by random(minRandomt3,1)  which is to say that a value
minRandomt3 = .95 picks a random value from .95 to 1 for the t3 parameter...randomP3 is not to be confused also with 
rendering a substring bezier curve from (0,1) but instead approximates a new curve randomly by computing a new 
control point from the old bezier curve data at parameter t3 and then randomly assigns a value random (0,t3) for t2
which is then applied yet again in computing t1 by random(0,t2)  and so forth.  This produces a distinct bezier curve
from the original bezier curve data and thus adds more randomness into the system.  
