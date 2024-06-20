function aux = rhsVFE(s, XTnb, c0, t)
% X = XTnb(1:3);
T = XTnb(4:6);
n= XTnb(7:9);
b = XTnb(10:12);
aux = [T;
    c0/sqrt(t)*n;
    -c0/sqrt(t)*T+s/(2*t)*b;
    -s/(2*t)*n];