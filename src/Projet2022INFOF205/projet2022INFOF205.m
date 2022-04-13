% This version of the solution works with yhatE(step) and slopehatE(step)
% throughout the iteration steps
if exist('tmax') ~= 1, tmax = 6; end
if exist('nsteps') ~= 1, nsteps = 4000; end
if exist('stepsize') ~= 1, stepsize = tmax/nsteps; end
if exist('playvideo') ~= 1, playvideo = false; end
nsteps = tmax/stepsize;
wheelcolor = [6 3 2]/10;
rand('state',matricule)

t0 = 0;
h = stepsize;
t = t0+h*(0:nsteps);
x0 = 0;
v0 = rand;
y0 = roadprofile(x0);
g = 9.81;

xhatE = zeros(1,nsteps+1); xhatE(1) = x0;
yhatE = zeros(1,nsteps+1); slopehatE = zeros(1,nsteps+1);
vxhatE = zeros(1,nsteps+1); vxhatE(1) = v0;
direction = sign(v0);
for step = 1:nsteps,
   [yhatE(step) slopehatE(step)] = roadprofile(xhatE(step)); 
   % à compléter ...
   if vxhatE(step)==0,
      axhat = % à compléter ...
      xhatE(step+1) = xhatE(step) + h^2*direction*axhat/2;
   end
end
vxhatE(nsteps+1) = 0;


if playvideo==true, figure(1)
   % L'objectif des calculs suivants est l'identification de la surface sur
   % laquelle la roue roule. La function f(x) décrivant les positions du centre
   % de la roue, il faut d'abord vérifier que le rayon de la roue ne soit pas
   % trop grand. (voir le calcul de maxr = rayon maximal).
   % Ensuite on calcule les positions des points de tangente depuis les
   % positions du centre et les pentes
   % (f(x) et f'(x), voir yhatE(xhatE) et slopehatE)
   dslope = initmom(slopehatE); dx = initmom(xhatE);
   secondderivative = row(dslope./dx);
   curvature = secondderivative./(1+slopehatE.^2).^(3/2);
   maxcurvature = max(-curvature);
   maxr = 1/maxcurvature;
   r0 = 0.03; r = min(r0,maxr);
   screensize = get(0,'screensize'); screensize = screensize(3:4);
   figsize = [560*r0/r 420]; figsize = min(figsize,screensize);
   figpaperpos = [0 0 figsize/96];
   figpos = get(gcf,'position'); figpos(3:4) = figsize;
   set(gcf,'position',figpos,'paperposition',figpaperpos)

   xx = xhatE+r*slopehatE./sqrt(1+slopehatE.^2);
   yy = yhatE-r./sqrt(1+slopehatE.^2);
   rollonhills(xx,yy,vhatE,0.03,wheelcolor)
end

