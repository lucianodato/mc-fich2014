cant = 5;
A= zeros(cant,cant);
b= zeros(cant,1);

A= [132 27.5 0 0 0;
    -27.5 105 22.5 0 0;
    0 -27.5 105 22.5 0;
    0 0 -27.5 105 22.5;
    0 0 0 -27.5 127.5];

%Phi tiempo 0.
% phi0=0;
% phi1= 0;
% phi2= 0;
% phi3= 0;
% phi4= 0;
% phi5= 0;
% dif_phi6 = 0;

% %phi tiempo 1
% phi0=0;
% phi1=0.0010;
% phi2=0.0024;
% phi3=-0.0009;
% phi4=0.0157;
% phi5=-0.0656;
% dif_phi6 = 1;

%phi tiempo 2
% phi0=0;
% phi1=0.0010;
% phi2=0.0070;
% phi3=-0.0114;
% phi4=0.0541;
% phi5=-0.1266;
% dif_phi6 = 1;

% %phi tiempo 3
% phi0=0;
% phi1=-0.0033;
% phi2=0.0174;
% phi3= -0.0343;
% phi4=0.0968;
% phi5= -0.1438;
% dif_phi6 = 1;


%phi tiempo 4
phi0=0;
phi1=-0.0117;
phi2=0.0346;
phi3=-0.0629;
phi4=0.0969;
phi5= -0.1384;
dif_phi6 = 1;


%phi tiempo 5
% phi0=0;
% phi1=-0.0243;
% phi2=0.0537;
% phi3=-0.0767;
% phi4=0.1122;
% phi5= -0.1379;
% dif_phi6 = 1;


st1=100*phi1 - 0.5*(-100*phi0 - 10*(phi0 - phi1) + 50*(phi1 + phi2) - 5*(phi3 - phi2));
st2=100*phi2 - 0.5*(-50*(phi1+phi2) - 5*(phi1 - phi2) + 50*(phi2 + phi3) - 5*(phi3 - phi2));
st3=100*phi3 - 0.5*(-50*(phi2+phi3) - 5*(phi2 - phi3) + 50*(phi3 + phi4) - 5*(phi4 - phi3));
st4=100*phi4 - 0.5*(-50*(phi3+phi4) - 5*(phi3 - phi4) + 50*(phi4 + phi5) - 5*(phi5 - phi4));
st5=100*phi5 - 0.5*(-50*(phi4+phi5) - 5*(phi4 - phi5) + 10*dif_phi6 + 100*phi5 -dif_phi6);

b=[(1/5 + st1) (1/5 + st2) (1/5 + st3) (1/5 + st4) (-8.8 + st5)];

a= A\b';

plot(a,'--or');

