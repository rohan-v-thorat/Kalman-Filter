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
ufull = filter(bb,aa,utmp);

%% System definition
% Continuous time

true_m = 1;
k = 100;
c = 0.4;

true_omg = sqrt(k/true_m);
true_zeta = c/(2*true_omg*true_m);

true_theta = [true_m; true_omg; true_zeta];

[Ac, Bc, C, D] = gen_ss(true_theta);

%% Convert continuous to discrete
Cfull = [eye(2*ndof); C];
Dfull = [ 0 ; 0 ; 1/true_m];
sysc = ss(Ac,Bc,Cfull,Dfull);
sysd = c2d(sysc,dt);
Ad = sysd.A;
Bd = sysd.B;

%% True signal
z0 = [0;0];
outputt = lsim(sysd,ufull',t,z0);

qdisp = outputt(:,1);
qvel = outputt(:,2);
qddot = outputt(:,3);
 

% Add noise to measurements
noise_percentage = 1;
noise_std = noise_percentage/100*rms(qddot);
noise = noise_std*randn(size(qddot));
y = qddot + noise;

%% discretized the noise Q
Q = 1e-3*eye(2);
R = 10;
Q_chol = chol(Q);
R_chol = chol(R);

%initialized
M = [0;0];
P = 1e-3*diag([1 1]);

% Track and animate 
MM = zeros(size(M,1),size(y',2));
PP = zeros(size(M,1),size(M,1),size(y',2));
% clf;
clc;
disp('In this demonstration we estimate a displacement and velocity from noisy measurements by using the classical Kalman filter.');
disp(' ');
disp('The filtering results are now displayed sequantially for 10 time step at a time.');
disp(' ');

n = 2; % size of state
m = 1; % size of measurement y
q = 100; % number of samples
x_samples = M + Q_chol*randn(n,q);

for k = 1:size(y',2)
    u = ufull(k);
    y_meas = y(k);

    % prediction/forecast step
    for j = 1:q
        x_samples(:,j) = f(x_samples(:,j),u);
    end
    
    x_mean = mean(x_samples,2);
    
    E_xx = x_samples - x_mean;

    % update/anaylsis step
    for j = 1  
        y_samples = g(x_samples(:,j),u);
    end    

    y_meas_samples = y_meas + R_chol*randn(m,q);
    y_meas_mean = mean(y_samples,2); % check whether to take mean of y_samples or y_meas_samples
    
    E_yy = y_meas_samples - y_meas_mean;
    P_yy = 1/(q-1)*(E_yy*E_yy');
    P_xy = 1/(q-1)*(E_xx*E_yy');

    K = P_xy/P_yy;
    
    x_samples = x_samples + K*(y_meas_samples-y_samples);
    x_mean = mean(x_samples,2);

    E_xx = x_samples - x_mean;
    P_xx = 1/(q-1)*(E_xx*E_xx');

    % store the results
    MM(:,k) = x_mean;
    PP(:,:,k) = P_xx;

 end
 
 %% plots
 figure(2); clf;
 subplot(3,1,1);
 plot(t,qdisp,'b',t,MM(1,:),'r','linewidth',1.2)
 title('displacement')
 legend('True','Estimated')
 grid on;

 subplot(3,1,2);
 plot(t,qvel,'b',t,MM(2,:),'r', 'linewidth',1.2)
 title('velocity')
 grid on;

 subplot(3,1,3);
 plot(t, qddot, 'b', t, (C*MM+D*u)', 'r', 'linewidth', 1.2);
 title('acceleration')

 grid on;

%% For understanding purpose
%  figure(2)
%  plot(t,utmp)
%  title('utmp')
%  figure(3)
%  plot(t,u)
%  title('u')

%% Functions

function output = f(x,u)
    
    dt = 0.01;
    m = 1;
    c = 0.4;
    k = 100;
    Ac = [0 1;-k/m -c/m];
    Bc = [0;1/m];
    Ad = expm(Ac*dt);
    Bd = Ac\(Ad-eye(2))*Bc;
    
    output = Ad*x + Bd*u;

end  

function output = g(x,u)

    m = 1;
    c = 0.4;
    k = 100;
    Ac = [0 1;-k/m -c/m];
    Bc = [0;1/m];
           
    output = Ac(2,:)*x + Bc(2,:)*u;
end