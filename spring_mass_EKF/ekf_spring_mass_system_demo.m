clear;clc;
% close all;

% w = omega, Z = zetha
ndof = 1;
dt = 0.01;    % sampling period
fs = 1/dt;    % 
T = 20;       % final time
t = 0:dt:T;   % time samples 

%% Generate input
rng(0);
utmp = 10*randn(1,length(t));
[bb,aa] = butter(5, 5/(fs/2),'low');
u = filter(bb,aa,utmp);
rng('shuffle');
%% System definition
% Continuous time

true_m = 1;
k = 100;
c = 0.4;

true_omg = sqrt(k/true_m);
true_zeta = c/(2*true_omg*true_m);

true_theta = [true_m; true_omg; true_zeta];

phi = log(true_theta);

[Ac, Bc, C, D] = gen_ss(phi);

%% Convert continuous to discrete
Cfull = [eye(2*ndof); C];
Dfull = [ 0 ; 0 ; 1/true_m];
sysc = ss(Ac,Bc,Cfull,Dfull);
sysd = c2d(sysc,dt);
Ad = sysd.A;
Bd = sysd.B;

%% True signal
z0 = [0;0];
outputt = lsim(sysd,u',t,z0);

qdisp = outputt(:,1);
qvel = outputt(:,2);
qddot = outputt(:,3);
qmass = zeros([size(t,2),1]) + 1;
qomega = zeros([size(t,2),1]) + 10;
qzeta = zeros([size(t,2),1]) + 0.02;



% Add noise to measurements
noise_percentage = 10;
noise_std = noise_percentage/100*rms(qddot);
noise = noise_std*randn(size(qddot));
y = qddot + noise;

%% discretized the noise Q
% q = 0.1;
% Q = diag([0 q]);
% L = eye(2);
% n   = size(Ac,1);
% Phi = [Ac L*Q*L'; zeros(n,n) -Ac'];
% AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
% Q   = AB(1:n,:)/AB((n+1):(2*n),:);



%initialized

%M = [0;0];
%P = diag([0 0.1]);

m = 2; % true m = 1
omega = 15; % true omega = 10
zeta = 0.03; % true zeta = 0.02 

phi = log([m; omega; zeta]);

% Ad = [1 0 ; 0 1] + [0 0; -exp(phi(2))^2 -2*exp(phi(3))*exp(2)]*dt;
% Bd = [0 ; exp(-phi(1))]*dt;
% Cd = [-exp(phi(2))^2  -2*exp(phi(3))*exp(2) ];
% Dd = [exp(-phi(1))];


x = [0; 0];
M = [ phi; x];
P = diag([0.1; 0.5; 0.1; 0.2; 0.1]);
% Q =  1*diag([0.001; 0.001; 0.001; 0.001; 0.001]);
Q =  1*diag([0.01; 0.001; 0.005; 0.001; 0.001]);
R = 100;

% Track and animate 
MM = zeros(size(M,1),size(y',2));
PP = zeros(size(M,1),size(M,1),size(y',2));
% clf;
clc;
disp('In this demonstration we estimate a mass, omega, zeta, displacement and velocity from noisy measurements by using the classical Kalman filter.');
disp(' ');
disp('The filtering results are now displayed sequantially for 10 time step at a time.');
disp(' ');

%rng('shuffle');
for k=1:size(y',2)
    [M,P] = ekf_predict(M,P,u(k),dt,Q);
    [M,P] = ekf_update(M,P,y(k),u(k),R);

    MM(:,k) = M;
    PP(:,:,k) = P;
    CC(:,k) = [-exp(M(2))^2  -2 * exp(M(3)) * exp(M(2))];
    DD(:,k) = [exp(-M(1))];

 end
z = CC .* [MM(4,:); MM(5,:)];
z1 = z(1,:);
z2 = z(2,:);

est = (z1 + z2 + (DD.*u))';
 %% plots
 figure(1); clf;
 subplot(3,1,1);
 plot(t,qdisp,'b',t,MM(4,:),'r','linewidth',1.2)
 title('displacement')
 legend('True displacement','Estimated displacement')
grid on;
 subplot(3,1,2);
 plot(t,qvel,'b',t,MM(5,:),'r', 'linewidth',1.2)
 title('velocity')
 legend('True velocity','Estimated velocity')
 grid on;

subplot(3,1,3);
plot(t, qddot, 'k', t, est , 'r', 'linewidth', 1.2);
title('acceleration')
legend('True','Estimated')
grid on;
xlabel('Time ( in seconds )')

 


 figure(2); clf;
 subplot(3,1,1);
 plot(t,qmass,'b',t,exp(MM(1,:)),'r', 'linewidth',1.2)
 title('mass')
 legend('True mass','Estimated mass')
 grid on;
 ylim([0,2]);

 subplot(3,1,2);
 plot(t,qomega,'b',t,exp(MM(2,:)),'r','linewidth',1.2)
 title('omega')
 legend('True omega','Estimated omega')
 grid on;
 ylim([0,20]);

 subplot(3,1,3);
 plot(t,qzeta,'b',t,exp(MM(3,:)),'r','linewidth',1.2)
 title('zeta')
 legend('True zeta','Estimated zeta')
 grid on;
 ylim([0,0.1]);
 xlabel('Time ( in seconds )')

% check whether you are comparing same true and estimated state or with
% next state !!




 %% For understanding purpose
%  figure(2)
%  plot(t,utmp)
%  title('utmp')
%  figure(3)
%  plot(t,u)
%  title('u')