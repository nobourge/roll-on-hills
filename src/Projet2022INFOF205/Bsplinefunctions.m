
function Phi = Bsplinefunctions(x,knots,degree,ext)

% Bsplinefunctions - version independent from software package
% Usage
%    Phi = Bsplinefunctions(x,knots,degree,ext);
% Inputs
%    x       points in which splines are evaluated, vector of length n
%    knots   knots of the spline basis, vector of length nk
%    degree  degree of polynomials
%    ext     (optional) binary value (default true) 
%            By setting ext=false one generates the interior basis functions
%            only. The number of basis functions then equals the number of
%            knots (= nk)
%            ext=true means that the basis will be extended so that the basis
%            functions sum up to one. The number of basis functions then equals
%            nk+degree-1
% Outputs
%     Phi    matrix of size n x (nk+degree-1) (if ext==true)
%                           n x nk            (if ext==false)
% Description
%     Applies recursion on parameter degree for B splines
% Note
%     It is assumed that the boundaries contain coinciding knots
% See Also
%          

if nargin<4, ext = true; end
if size(x,1)<size(x,2), x=x'; end
if size(knots,1)<size(knots,2), knots=knots'; end
n = length(x);
nk = length(knots);
if ext==true, nPhi = nk-1; else, nPhi = nk; end 
% at the end, we will have for ext==true that nPhi = nk+degree-1;
if degree < 1
   Phi = zeros(n,nPhi); 
   knots = [knots; Inf];
   for k=1:nPhi
      Phi(1:n,k) = (x>=knots(k) & x<knots(k+1));
   end
else
   Phi = zeros(n,nk);
   for k=2:nk 
      i=find(x>knots(k-1)&x<=knots(k));
      Phi(i,k) = (x(i)-knots(k-1))/(knots(k)-knots(k-1));
   end
   for k=1:nk-1
      i=find(x>=knots(k)&x<knots(k+1));
      Phi(i,k) = (knots(k+1)-x(i))/(knots(k+1)-knots(k));
   end
   nPhi = nk;
   for s = 3:degree+1
      Phi0 = Phi;
      t = floor(s/2); r = s-2*t;
      for k = 2-r:t
         Phi(1:n,k) = ...
            (x-knots(1))/(knots(k-t+s-1)-knots(1)).*Phi0(1:n,k+r-1) + ...
            (knots(k-t+s)-x)/(knots(k-t+s)-knots(1)).*Phi0(1:n,k+r);
      end
      for k = t+1:nPhi-s+t
         Phi(1:n,k) = ...
            (x-knots(k-t))/(knots(k-t+s-1)-knots(k-t)).*Phi0(1:n,k+r-1) + ...
            (knots(k-t+s)-x)/(knots(k-t+s)-knots(k-t+1)).*Phi0(1:n,k+r);
      end
      for k = nPhi-s+t+1:nPhi-r
         Phi(1:n,k) = ...
            (x-knots(k-t))/(knots(nPhi)-knots(k-t)).*Phi0(1:n,k+r-1) + ...
            (knots(nPhi)-x)/(knots(nPhi)-knots(k-t+1)).*Phi0(1:n,k+r);
      end
      if r==1
         Phi(1:n,nPhi) = ...
         (x-knots(nPhi-t))/(knots(nPhi)-knots(nPhi-t)).*Phi0(1:n,nPhi);
         if ext==true
            % (k = 0)
            phi0 = (knots(-t+s)-x)/(knots(-t+s)-knots(1)).*Phi0(1:n,1);
            knots = [knots(1); knots];
            Phi = [phi0 Phi];
            nPhi = nPhi+1;
         end
      end
      if r==0
         Phi(1:n,1) = (knots(1-t+s)-x)/(knots(1-t+s)-knots(1)).*Phi0(1:n,1);
         if ext==true
            k = nPhi+1;
            phi1 = (x-knots(k-t))/(knots(nPhi)-knots(k-t)).*Phi0(1:n,nPhi);
            knots = [knots; knots(nPhi)];
            Phi = [Phi phi1];
            nPhi = nPhi+1;
         end
      end
   end
end
% Phi(xorder,1:nk) = Phi;
