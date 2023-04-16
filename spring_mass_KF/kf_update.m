function [M,P] = kf_update(M,P,y,C,D,u,R)

  
  S = (C*P*C'+ R);
  K = P*C'/S;
  M = M + K * (y-C*M-D*u);
  P = P - K*S*K';
  

