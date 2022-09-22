function [M,P] = ekf_predict(M,P,u,dt,Q)
    % Read  'der_Ad_phi' as derivative of Ad w.r.t phi

    X_multiply_I = [ M(4) 0 M(5) 0; 0 M(4) 0 M(5)];
    der_Ad_phi = [ 0 0 0; 0 -2*exp(M(2)) 0;...
        0 0 0; 0 -2*exp(M(3))*exp(M(2)) -2*exp(M(3))*exp(M(2))]*dt;
    u_multipy_I = [u 0; 0 u];
    der_Bd_phi = [0; exp(-M(1))]*dt;

% Ad, Bd, der_f2_phi, der_f2_x
    Ad = [1 0 ; 0 1] + [0 1; -exp(M(2))^2 -2*exp(M(3))*exp(M(2))]*dt;
    Bd = [0 ; exp(-M(1))]*dt;
    der_f2_phi = X_multiply_I * der_Ad_phi + u_multipy_I * der_Bd_phi;
    der_f2_x = Ad;

    Fz = [eye(3) zeros([3,2]);  der_f2_phi der_f2_x];

    % update
    M = [M(1); M(2); M(3); Ad * [M(4); M(5)] + Bd * u];
    P = Fz * P * Fz' + Q;

    % Here can we/should we use predicted 'M' to obtain 'Fz' value?
  
