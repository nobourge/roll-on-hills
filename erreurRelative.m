d = -1:0.1:3;
d_hat=d-0.01;



function r=err_relat(x, x_hat)
    if x ~= 0
        r = abs(x_hat-x)/abs(x);
        disp(r);
    end
end 



function Fd = F(d)
    Fd=d-1.^2;
end

x=F(d);
x_hat=F(d_hat);

y = err_relat(x, x_hat);
plot(d,y);
hold