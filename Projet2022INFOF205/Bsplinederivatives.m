
function DPhi = Bsplinederivatives(x,knots,degree,ext);

% Bsplinederivative - version independent from software package
% Usage
%    DPhi = Bsplinederivatives(x,knots,degree,ext);
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
%    DPhi    matrix of size n x nk 
% Description
% Note
%    The routine assumes that the first and last knots have a multiplicity to
%    deal with all definitions involving knots with numbers outside {1,...,n}
% See Also
%     help Bsplinefunctions

if size(x,1)<size(x,2), x=x'; end
if size(knots,1)<size(knots,2), knots=knots'; end
n = length(x);
nk = length(knots);
Phi = Bsplinefunctions(x,knots,degree-1,true);
q = ceil(degree/2); r = degree+1-2*q;
knots = [ones(degree-1,1)*knots(1);knots;ones(degree-1,1)*knots(nk)];
nPhi = nk+degree-2;
Delta = knots(degree+1:degree+nPhi)-knots(1:nPhi);
D = [zeros(nPhi,1) diag(1./Delta)]-[diag(1./Delta) zeros(nPhi,1)];
DPhi = degree*Phi*D;
if nargin<4, ext=true; end
if ext==false, DPhi = DPhi(1:n,q+1:q+nk); end
