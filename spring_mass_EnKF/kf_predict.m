function [M,P] = kf_predict(M,P,Ad,Bd,u,Q)

    M = Ad * M + Bd * u;
    P = Ad * P * Ad' + Q;
  
