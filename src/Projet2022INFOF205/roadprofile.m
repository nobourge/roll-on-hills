
function [y, Dy] = roadprofile(x)

y = cos(3*pi*x.^2).*(1-x-x.^2);
Dy = cos(3*pi*x.^2).*(-1-2*x)-sin(3*pi*x.^2).*(6*pi*x).*(1-x-x.^2);
y = y/10;
Dy = Dy/10;
