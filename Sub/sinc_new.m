function f_val = sinc_new(t,f,type)
%SINC_NEW Summary of this function goes here
%   Detailed explanation goes here

w = pi*f;
x = w*(t+eps);

s = w^2/3;

switch type
    case 0
        if x == 0
            f_val = 1;
        else
            f_val = sin(x)/x;
        end
    case 1
        if t == 0
            f_val = 0;
        else
            f_val = (cos(x)*x-sin(x))/x^2*w/sqrt(s);
        end
    case 2
        if t == 0
            f_val = 1;
        else
            f_val = (2*x*cos(x)-(2-x^2)*sin(x))/x^3/s*w^2;
        end
end


end
