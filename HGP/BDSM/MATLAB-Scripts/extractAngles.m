function [res] = extractAngles(a,b,c) 

res(1) = acos((b^2 + c^2 - a^2)/(2*b*c));
res(2) = acos((-b^2 + c^2 + a^2)/(2*a*c));
res(3) = pi - res(1) - res(2);



    