
function [x00, y00, xoo, yoo] = rollonhills(x,y,v,r,color,frames)

% rollonhills; displays (x,y) as wheel rolling on hill at given speed v
% Usage
%    [x00 y00 xoo yoo] = rollonhills(x,y,vo,r,color,frames)
% Inputs
%    x,y,vo  real vectors of equal lengths n; 
%            (x,y) positions of point of tangency; 
%            vo    speed at wheel center
%    r       positive real number; 
%            radius of the wheel
%    color   (default 'k','black',[0 0 0]) Matlab color specicfication
%    frames  time points in which position is displayed (default is equidistant
%            vector at a rate of 1/24 seconds)
% Outputs (all optional)
%    (x00,y00) shift lower left corner of figure towards origin
%    (xoo,yoo) displayed positions of wheel center
% Description
%
% Note
%    (1)
%    The framerate (i.e., the time gaps between the frames) affects the video
%    quality of the simulation. Too low framerates (i.e., too large gaps, too
%    small vector of frames) may cause aliasing. Too high framerates may slow
%    down Matlab's calculations, leding to unrealistically slow motion.
%    (2)
%    Calculations are approximative, ignoring 2nd order effects. This results
%    in slipping effects in highly curved valleys.
% See Also
%          

%% arg in
if nargin<6, frames = NaN; end
if nargin<5, color = 'k'; end
if nargin<4, r = NaN; end
if nargin<3, v = 1; end
if ischar(r) || length(r)==3
   color = r; 
   r = NaN; 
%elseif length(r)==3
 %  color = r; r = NaN;
end
wheelopt = {'color',color,'linewidth',2};

if size(x,1)==1, x=x'; end 
if size(y,1)==1, y=y'; end
if size(v,1)==1, v=v'; end

n = length(x);
if length(y) == 2
   y = (x-x(1))*(y(2)-y(1))/(x(n)-x(1));
elseif length(y)~=n
   error('lengths x and y dismatch')
end
if length(v)==1, v = v*ones(size(x)); end

%% Find elapsed time and total distance
dx = x(2:n)-x(1:n-1); 
dy = y(2:n)-y(1:n-1);
ds = sqrt(dx.^2+dy.^2).*sign(dx);
s = [0; cumsum(ds)];
vbar = (v(1:n-1)+v(2:n))/2;
dt = abs(ds)./vbar; 
dt(isnan(dt))=0; dt(isinf(dt))=0;
t = [0; cumsum(dt)];

if isnan(r)
    r = s(n)/40;
end
framerate = 24;
framerate = framerate*8;
frametime = 1/framerate;
nframes = floor(framerate*t(n));
if isnan(frames)
    frames = (0:nframes-1)'/(nframes-1);
elseif max(frames)>1
    frames = min(1,frames/t(n)); 
end
nframes = length(frames);
tt = frames*t(n);
Phi = Bsplinefunctions(tt,t,1);
xx = Phi*x;
yy = Phi*y;
vv = Phi*v;
%
ss = Phi*s;
[xxx, i] = unique([xx;x]);
yyy = [yy;y]; yyy = yyy(i);
xxx = [xxx(1)-r;xxx;xxx(end)+r]; 
yyy = [yyy(1);yyy;yyy(end)];
ymin = min(yyy); 
ymax = max(yyy); 
y0 = ymin*1.1-0.1*ymax; 
y1 = ymax+2*r-y0;
%% shift plot so that left lower corner is at the origin 
x0 = min(xxx); 
xxx=xxx-x0;
yyy=yyy-y0;
xx=xx-x0;
yy=yy-y0;
PhiD = Bsplinederivatives(tt,t,1);
dxx = PhiD*x; dyy = PhiD*y;
dxx = initmom(xx); 
dyy = initmom(yy);
xxo = xx+cos(atan(-dxx./dyy)).*sign(-dxx./dyy)*r;
yyo = yy+sin(atan(-dxx./dyy)).*sign(-dxx./dyy)*r;

alfa = -ss/r;

w = 100; 
wheel = (0:w)/w*2*pi; 
wheelx = r*cos(wheel); 
wheely = r*sin(wheel);
for f=1:nframes
   xf = xx(f); 
   yf = yy(f); 
   xo = xxo(f); 
   yo = yyo(f);
   spikes = alfa(f)+(0:5)'*pi/3;
   spikes = [xo+r*cos(spikes) yo+r*sin(spikes)];
   fill([xxx;xxx(end);xxx(1)],[yyy;0;0],[1 1 1]*0.4)
   hold on
   %figure(3)
   plot(xxx,yyy,'k','linewidth',2)
   plot(xo+wheelx,yo+wheely,wheelopt{:})
   fill(xo+wheelx/3,yo+wheely/3,color)   
   plot(xo+wheelx/3,yo+wheely/3,'color',color)
   for s=1:3 
       ss = [s s+3]; 
       plot(spikes(ss,1),spikes(ss,2),wheelopt{:}); 
   end
   hold off
   grid
   axis equal
   axis([xxx(1) xxx(end) 0 y1])
   axis off
   pause(frametime)
end
%% arg out
if nargout>0
   x00 = x0;
end
if nargout>0
   y00 = y0;
end
if nargout>2
   xoo = xxo;
end
if nargout>3
   yoo = yyo;
end
