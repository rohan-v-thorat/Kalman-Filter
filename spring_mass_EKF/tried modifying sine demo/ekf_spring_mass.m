% Demonstration for EKF using a random sine signal model. 
%
% A Very simple demonstration for extended Kalman filter (EKF), which is
% used to track a random single-component sinusoid signal,
% which is modelled as x_k = a_k*sin(\theta_k), dtheta/dt = omega_k.
% The signal is also filtered with unscented Kalman filter (UKF) for
% comparison.
%
% Copyright (C) 2007 Jouni Hartikainen
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

clc;
disp('Filtering the signal with EKF...');

save_plots = 1;

% Measurement model and it's derivative
h_func = @ekf_sine_h;
dh_dx_func = @ekf_sine_dh_dx;


% Initial values for the signal.
f = 0;
w = 10;
a = 1;
  
% Number of samples and stepsize.
d = 5;
n = 500;
dt = d/n;
x = 1:n;

% Check the derivative of the measurement function.
der_check(h_func, dh_dx_func, 1, [f w a]');

% Dynamic state transition matrix in continous-time domain.
F = [0 1 0;
     0 0 0;
     0 0 0];
  
% Noise effect matrix in continous-time domain.
L = [0 0;
     1 0;
     0 1];
  
% Spectral power density of the white noise.
q1 = 0.2;
q2 = 0.1;
Qc = diag([q1 q2]);
  
% Discretize the plant equation.
[A,Q] = lti_disc(F,L,Qc,dt);
  
% Generate the real signal.
X = zeros(3, n);
X(:,1) = [f w a]';
for i = 2:n   % it is process model ?
   X(:,i) = A*X(:,i-1) + gauss_rnd([0 0 0]', Q);
end  
  
% Generate the observations with Gaussian noise.
sd = 1;
R = sd^2;

Y = zeros(1,n);
Y_real = feval(h_func,X);     
Y = Y_real + gauss_rnd(0,R,n);
  
plot(x,Y,'.',x,Y_real)
  
% Initial guesses for the state mean and covariance.
M = [f w a]';
P = diag([3 3 3]);    
  
% Reserve space for estimates.
MM = zeros(size(M,1),size(Y,2));
PP = zeros(size(M,1),size(M,1),size(Y,2));

% Estimate with EKF
for k=1:size(Y,2)
   [M,P] = ekf_predict1(M,P,A,Q);
   [M,P] = ekf_update1(M,P,Y(:,k),dh_dx_func,R*eye(1),h_func);
   MM(:,k)   = M;
   PP(:,:,k) = P;
end
clf; clc;
disp('The filtering results using the 1st order EKF is now displayed')

% Project the estimates to measurement space 
Y_m = feval(h_func, MM);

plot(x,Y,'.', x,Y_real,'--',x,Y_m);
legend('Measurements','Real signal', 'Filtered estimate');
xlim([0 ceil(max(x))]);
title('Estimating a random Sine signal with extended Kalman filter.');

%if save_plots
    %print -dpsc demo2_f1.ps 
%end

clc;
disp('The filtering result using the 1st order EKF is now displayed.');



