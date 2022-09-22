function [M,P] = kf_update(M,P,y,C,R)

  
  S = (C*P*C'+ R);
  K = P*C'/S;
  M = M + K * (y-C*M);
  P = P - K*S*K';
  

