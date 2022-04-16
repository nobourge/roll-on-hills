% This version of the solution works with yhatE(step) and slopehatE(step)
% throughout the iteration steps
%if exist('tmax','var') ~= 1, tmax = 6; end %original
if exist('tmax','var') ~= 1, tmax = 666; end

%if exist('nsteps','var') ~= 1, nsteps = 4000; end %original
if exist('nsteps','var') ~= 1, nsteps = 40000; end
nsteps
if exist('stepsize','var') ~= 1, stepsize = tmax/nsteps; end
%stepsize = stepsize*100; %Requested array exceeds the maximum possible variable size.
if exist('playvideo','var') ~= 1, playvideo = true; end

nsteps = tmax/stepsize; %original
wheelcolor = [6 3 2]/10;
%%%rand('state',matricule) %Unrecognized function or variable 'matricule'.
rand('state', 496667)
%rand('state', 9)
%rand('state', 999)
%rand('state', 100000)

%%%rng('state',496667) % First input must be a nonnegative integer seed less than 2^32, 'shuffle', 'default', or generator settings captured previously using S = RNG.
%rng(496667)
%rng(seed,'state',496667)
%%%rng(seed, 'state') %Unrecognized function or variable 'seed'.

t0 = 0;
h = stepsize;
t = t0+h*(0:nsteps);
x0 = 0;
v0 = rand;
y0 = roadprofile(x0);
gT = 9.81;
gL = 1.62;
g = gL;

%% EDO resolution Euler explicit method (progressive) 
nsteps
xhatE = zeros(1,nsteps+1); 
xhatE(1) = x0;
yhatE = zeros(1,nsteps+1); 
slopehatE = zeros(1,nsteps+1);
vxhatE = zeros(1,nsteps+1); 
vxhatE(1) = v0;
direction = sign(v0);
for step = 1:nsteps
   %%% position, pente
   [yhatE(step), slopehatE(step)] = roadprofile(xhatE(step)); 
   %%% à compléter ...
   %%%|x^'(t)|²
   %%% vxhat2 = (2.*g.*(y0-yhatE(step)+v0.^2))./(1+(slopehatE(step)).^2);
   vxhat2 = (2*g*(y0-yhatE(step)+v0^2))/(1+(slopehatE(step))^2);
   if vxhat2 < 0
       %la methode numerique entre dans une “zone interdite” 
       %par la physique (et la mathematique).s 
       %Il faut alors inverser la direction du mouvement en evaluant
      direction = -sign(slopehatE(step));
   end
   if vxhat2 < 0
       %%% la valeur est remplacee par zero.
   end
   vxhat = sqrt(vxhat2);
   
   %%% xhat(t+h)= xhat(t) + h * vxhat = xhat(t) + h * sign(vxhat)*vxhat
   %xhatE(step+1)= xhatE(step) + h * vxhat;
   xhatE(step+1)= xhatE(step) + h*sign(vxhat)*vxhat;

   % completion fin
   
   %%%Par contre, si x^'(t) = 0
   if vxhatE(step)==0
      %%% il faut calculer l’acceleration,
      %%% à compléter ...
      axhat = g*(abs(slopehatE(step))/(1+(slopehatE(step))^2));
      %%% completion fin

      xhatE(step+1) = xhatE(step) + h^2*direction*axhat/2;
   end
end
vxhatE(nsteps+1) = 0;

if playvideo==true, figure(1)
   % L'objectif des calculs suivants est 
   % l'identification de la surface sur laquelle la roue roule. 
   % La function f(x) décrivant les positions du centre
   % de la roue, 
   % il faut d'abord vérifier que le rayon de la roue ne soit pas
   % trop grand. (voir le calcul de maxr = rayon maximal).
   % Ensuite on calcule les positions des points de tangente depuis les
   % positions du centre et les pentes
   % (f(x) et f'(x), voir yhatE(xhatE) et slopehatE)
   dslope = initmom(slopehatE); 
   dx = initmom(xhatE);
   secondderivative = row(dslope./dx);
   curvature = secondderivative./(1+slopehatE.^2).^(3/2);
   %curvature %Column <  4 000 : NaN; 4 000 & 4 001 : Inf
   maxcurvature = max(-curvature); %original
   %maxcurvature = max(curvature);
   %maxcurvature
   
   % il faut d'abord vérifier que le rayon de la roue ne soit pas
   % trop grand. (voir le calcul de maxr = rayon maximal).

   maxr = 1/maxcurvature;
   %maxr
   r0 = 0.03; 
   r = min(r0,maxr);
   %r
   screensize = get(0,'screensize'); 
   screensize = screensize(3:4);
   figsize = [560*r0/r, 420]; 
   %figsize = [999, 999]; 

   %figsize
   figsize = min(figsize,screensize);
   %figsize
   figpaperpos = [0, 0, figsize/96];
   figpos = get(gcf,'position'); 
   %figpos
   figpos(3:4) = figsize;
   %figpos
   %get current figure
 
   set(gcf,'position',figpos,'paperposition',figpaperpos)
   
   % Ensuite on calcule les positions des points de tangente depuis les
   % positions du centre et les pentes
   % (f(x) et f'(x), voir yhatE(xhatE) et slopehatE)

   xx = xhatE+r*slopehatE./sqrt(1+slopehatE.^2);
   %xx
   yy = yhatE-r./sqrt(1+slopehatE.^2);
   %yy
   vhatE = 1;
   %rollonhills(x,y, v,    r,   color,    frames)
   %rollonhills(xx,yy,vhatE,0.03,wheelcolor) %original
   rollonhills(xx,yy,vhatE,r,wheelcolor)

end

