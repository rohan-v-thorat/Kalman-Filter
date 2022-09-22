function [M,P] = ekf_update(M,P,y,u,R)

  
%   S = (C*P*C'+ R);
%   K = P*C'/S;
%   M = M + K * (y-C*M);
%   P = P - K*S*K';
x_multiply_I = [M(4) M(5)];
der_Cd_phi = [0 -2*exp(M(2)) 0; 0 -2*exp(M(2))*exp(M(3)) 0];
u_multipy_I = u;
der_Dd_phi = [-exp(M(1)) 0 0];

% Cd, Dd, der_h_phi, der_h_x
Cd = [-exp(M(2))^2  -2*exp(M(3))*exp(M(2))];
Dd = [exp(-M(1))];
der_h_phi = x_multiply_I * der_Cd_phi + u_multipy_I * der_Dd_phi;
der_h_x = Cd; % how we can get 1*2 dimensions % lets adjust for now(take Cd instead of Cd')

Hz = [der_h_phi der_h_x]; % here we need dim = 1*(3+2) = 1*5
 
v =  y - (Cd * [M(4); M(5)] + Dd * u);
s = Hz * P * Hz' + R;
k = P * Hz' /s;
M = M + k*v;
P = P - k * s * k';


