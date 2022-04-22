
function [y, Dy] = roadprofile(x)

%y = cos(3*pi*x.^2).*(1-x-x.^2); %original
%y = 1;
%y = -x;
y = sin(x);
%y = cos(x);

%y = x.^2.*sin(x);
%y = x.^2.*cos(x);

%Dy = cos(3*pi*x.^2).*(-1-2*x)-sin(3*pi*x.^2).*(6*pi*x).*(1-x-x.^2); %original
%%%Dy = 0;
%%%Dy = -1;
Dy = cos(x);
%Dy = -sin(x);

%Dy = 2.*x.*sin(x) + x.^2.*cos(x);
%Dy = 2.*x.*cos(x) + x.^2.*(-sin(x));

%y = y/10;
%Dy = Dy/10;
end
