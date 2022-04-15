
function [y, Dy] = roadprofile(x)
%x=-1:0.1:2;
%y = cos(3*pi*x.^2).*(1-x-x.^2); %original
y = x.^2.*cos(x);
Dy = cos(3*pi*x.^2).*(-1-2*x)-sin(3*pi*x.^2).*(6*pi*x).*(1-x-x.^2);
y = y/10;
Dy = Dy/10;

%plot(x,y);