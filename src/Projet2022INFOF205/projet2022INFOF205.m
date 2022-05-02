% This version of the solution works with yhatE(step) and slopehatE(step)
% throughout the iteration steps

% *  is matrix multiplication 
% .* is elementwise multiplication

clear variables;

if exist('tmax','var') ~= 1
    tmax = 6;    %%%original
    %tmax = 66;
    %tmax = 666; 
end

if exist('nsteps','var') ~= 1 
    nsteps = 4000; %max   %%%
end 
if exist('stepsize','var') ~= 1 
    stepsize = tmax/nsteps; 
end
%stepsize = stepsize*100;

if exist('playvideo','var') ~= 1
    playvideo = true; 
end

nsteps = tmax/stepsize;     %%%original
%stepsize = stepsize*100;
wheelcolor = [6 3 2]/10;    %%%
%%%rand('state',matricule) %Unrecognized function or variable 'matricule'.
%rand('state', 1)
%rand('state', 9)
%rand('state', 999)
%rand('state', 100000)
%rand('state', 123456)

rand('state', 496667) % matricule

%%%rng('state',496667) % First input must be a nonnegative integer seed less than 2^32, 'shuffle', 'default', or generator settings captured previously using S = RNG.
%rng(496667)
%rng(seed,'state',496667)
%%%rng(seed, 'state') %Unrecognized function or variable 'seed'.

h0 = 0;
%h0 = 2;
h = stepsize;   %%%
t0 = 0;     %%%
t = t0+h*(0:nsteps);    %%%

x0 = 0;     %%%
v0 = rand;  %%%
y0 = roadprofile(x0);   %%%

gL = 1.62;
gT = 9.81;  %%%
gS = 274;
g = gT;
%g = gL;
%g = gS;
m = 1;  %masse roue

%% EDO resolution Euler explicit method (progressive) 
Ecin = zeros(1,nsteps+1);
Epot = zeros(1,nsteps+1);
Etot = zeros(1,nsteps+1);

xhatE = zeros(1,nsteps+1);  %%%
xhatE(1) = x0;  %%%
yhatE = zeros(1,nsteps+1);  %%%
slopehatE = zeros(1,nsteps+1);  %%%
vxhatE = zeros(1,nsteps+1); %%%
vxhatE(1) = v0; %%%
vhatE = zeros(1,nsteps+1);
direction = sign(v0);   %%%


for step = 1:nsteps
   [yhatE(step), slopehatE(step)] = roadprofile(xhatE(step));   %%%
   %%% à compléter ...
   %%%|x^'(t)|²
   vxhat2 = (2*g*(y0-yhatE(step)+v0.^2))/(1+(slopehatE(step)).^2); 
  
   if vxhat2 < 0
       %la methode numerique entre dans une “zone interdite” 
       %par la physique (et la mathematique)
       %Il faut alors inverser la direction du mouvement en evaluant
      direction = -sign(slopehatE(step));
   end
  
   vxhatE(step)= abs(sqrt(vxhat2));
   
   %%% xhat(t+h)= xhat(t) + h * vxhat 
   %%%          = xhat(t) + h * sign(vxhat)*vxhat
   xhatE(step+1)= xhatE(step) + h*direction*vxhatE(step);

   % completion fin 
   
   if vxhatE(step)==0   %%%
      %%% x^'(t) = 0
      %%% il faut calculer l’acceleration,
      %%% axhat = à compléter ... 
      axhat = g*(abs(slopehatE(step))/(1+(slopehatE(step)).^2));    
      %%% completion fin 
      xhatE(step+1) = xhatE(step) + h^2*direction*axhat/2;  %%%
      %xhatE(step+1)
      %disp('---------')
   end
   %% vitesse en fonction de l'acceleration
   vhatE(step) = vxhatE(step)*sqrt(1+slopehatE(step).^2);
   %% energies
   Ecin(step) = (m*(vhatE(step).^2))/2; %
   Epot(step) = m*g*(yhatE(step)-h0);
   Etot(step) = Ecin(step) + Epot(step);
end

vxhatE(nsteps+1) = 0;   %%%

%% video
if playvideo==true, figure(1)   %%%
   % L'objectif des calculs suivants est 
   % l'identification de la surface sur laquelle la roue roule. 
   % La function f(x) décrivant les positions du centre
   % de la roue, 
  
   dslope = initmom(slopehatE); %%%
   dx = initmom(xhatE); %%%
   secondderivative = row(dslope./dx);  %%%
   
   %% road curvature
   curvature = secondderivative./(1+slopehatE.^2).^(3/2);   %%%
   maxcurvature = max(-curvature); %%%
 
   %% wheel radius
   % il faut d'abord vérifier que le rayon de la roue ne soit pas
   % trop grand. (voir le calcul de maxr (= rayon maximal)).
   maxr = 1/maxcurvature;   %%%
   r0 = 0.03;   %%%
   r = min(r0,maxr);    %%%
   screensize = get(0,'screensize');    %%%
   screensize = screensize(3:4);    %%%
   figsize = [560*r0/r, 420];   %%%
   %figsize = [960*r0/r, 520];  
   figsize = min(figsize,screensize);   %%%
   figpaperpos = [0, 0, figsize/96];    %%%
   figpos = get(gcf,'position');    %%%
   figpos(3:4) = figsize; %%%
   %GetCurrentFigure
   set(gcf,'position',figpos,'paperposition',figpaperpos)   %%%
   % Ensuite on calcule les positions des points de tangente %%%
   % depuis les %%%
   % positions du centre et les pentes %%%
   % (f(x) et f'(x), voir yhatE(xhatE) et slopehatE) %%%

   xx = xhatE + r * slopehatE  ./  sqrt(1+slopehatE.^2);    %%%
   yy = yhatE - r              ./  sqrt(1+slopehatE.^2);  %%%
   
   rollonhills(xx, yy,   vhatE,0.03,wheelcolor) %%%
   %rollonhillsenergiesviz(nsteps, xhatE, vxhatE, Epot, Ecin, Etot, xx,yy,vhatE,0.03,wheelcolor)
end

plot_energies(nsteps, xhatE, vxhatE, Epot, Ecin, Etot);


