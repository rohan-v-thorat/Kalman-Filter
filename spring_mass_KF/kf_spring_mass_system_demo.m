clear;clc;
% close all;

% w = omega, Z = zetha
ndof = 1;
dt = 0.01;    % sampling period
fs = 1/dt;    % 
T = 10;       % final time
t = 0:dt:T;   % time samples 

%% Generate input
rng(0);
utmp = 10*randn(1,length(t));
[bb,aa] = butter(5, 5/(fs/2),'low');
u = filter(bb,aa,utmp);

%% System definition
% Continuous time

true_m = 1;
k = 100;
c = 0.4;

true_omg = sqrt(k/true_m);
true_zeta = c/(2*true_omg*true_m);

true_theta = [true_m; true_omg; true_zeta];

[Ac, Bc, C, D] = gen_ss(log(true_theta));

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
 

% Add noise to measurements
noise_percentage = 10;
noise_std = noise_percentage/100*rms(qddot);
noise = noise_std*randn(size(qddot));
y = qddot + noise;

%% discretized the noise Q
q = 0.1;
Q = diag([0 q]);
L = eye(2);
n   = size(Ac,1);
Phi = [Ac L*Q*L'; zeros(n,n) -Ac'];
AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
Q   = AB(1:n,:)/AB((n+1):(2*n),:);

Q = 10000*diag([0.001; 0.01]);
R = 10;

%initialized
M = [0;0];
P = diag([0 0.1]);

% Track and animate 
MM = zeros(size(M,1),size(y',2));
PP = zeros(size(M,1),size(M,1),size(y',2));
% clf;
clc;
disp('In this demonstration we estimate a displacement and velocity from noisy measurements by using the classical Kalman filter.');
disp(' ');
disp('The filtering results are now displayed sequantially for 10 time step at a time.');
disp(' ');

for k=1:size(y',2)
    [M,P] = kf_predict(M,P,Ad,Bd,u(k),Q);
    [M,P] = kf_update(M,P,y(k),C,D,u(k),R);

    MM(:,k) = M;
    PP(:,:,k) = P;

 end
 
 %% plots
 figure(2); clf;
 subplot(2,1,2);
 plot(t,qvel,'b',t,MM(2,:),'r', 'linewidth',1.2)
 title('velocity')
 legend('True velocity','Estimated velocity')
 grid on;

 subplot(2,1,1);
 plot(t,qdisp,'b',t,MM(1,:),'r','linewidth',1.2)
 title('displacement')
 legend('True displacement','Estimated displacement')
grid on;

 figure(3); clf;
 plot(t, y, 'k', t, (C*MM+D*u)', 'r', 'linewidth', 1.2);
 title('acceleration')
 legend('observed','Estimated')
grid on;

figure(4); clf;
 plot(t, qddot, 'k', t, (C*MM+D*u)', 'r', 'linewidth', 1.2);
 title('acceleration')
 legend('True','Estimated')
grid on;


 %% For understanding purpose
%  figure(2)
%  plot(t,utmp)
%  title('utmp')
%  figure(3)
%  plot(t,u)
%  title('u')
