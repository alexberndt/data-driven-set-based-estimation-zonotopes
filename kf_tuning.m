%% Comparison of Matrix Zonotopes
%
% We are considering the system 
%  x(k+1) = Ax(k) + w(k)
%  z(k)   = x(k) + gamma(k)
%
% A = [0.9455   -0.2426;
%      0.2486    0.9455];
%     
% Result:
%  
%   x(k+1) = Ax(k) + w(k)
%   
%   and adhere to the data 
%
% author:  Alexander Berndt
% contact: alberndt@kth.se
%
%% Training phase to identify the system
clc
clear
%cd '/home/alberndt/Documents/research/data_driven/code/data_driven_set_based_estimation_zonotopes'
rng("default");
addpath('./functions/');

% Training phase

% Run the 'training' phase
A_true  = [0.9455   -0.2426;
           0.2486    0.9455];
B_true  = [0.1; 0];
T       = 500;
Ts      = 1.0;

ctrb(A_true,B_true)

useRandPointExtreme = false;

% Init matrices
x       = zeros(2,T+1);
z       = zeros(2,T);
t       = zeros(1,T+1);
gam     = zeros(2,T);
w       = zeros(2,T);
u       = zeros(1,T);

% Define initial condition
Z_X_0   = zonotope([0;0], blkdiag(15,15));
x(:,1)  = randPoint(Z_X_0);
x(:,1)  = [10.0;10.0];

% Define input zonotope
Z_u         = zonotope( 0, 10);

% Define noise zonotopes
c_gam       = 0.01;
c_w         = 0.01;
Z_gamma     = zonotope( [0.0;0.0], blkdiag(c_gam, c_gam) );
Z_w         = zonotope( [0.0;0.0], [c_w 0 0.2*c_w; 0 0.9*c_w -0.05*c_w] ); 

% Generate training data
for k = 1:T
    % random input sequence
    u(k)        = 1.0*randPoint(Z_u);
    % random bounded noise
    if useRandPointExtreme
        gam(:,k)    = randPointExtreme(Z_gamma);
        w(:,k)      = randPointExtreme(Z_w);
    else
        gam(:,k)    = randPoint(Z_gamma);
        w(:,k)      = randPoint(Z_w);
    end
    % system evolution
    x(:,k+1)    = A_true*x(:,k) + B_true*u(k) + w(:,k);
    z(:,k)      = x(:,k) + gam(:,k);
    t(k+1)      = t(k) + Ts;
end

% N4SID model
sysid_data = iddata(z',u',Ts);
m = n4sid(sysid_data,2,'ssp','can'); % canonical form for C = eye(2)
[pvec,pvec_sd] = getpvec(m);

unc_stddev = reshape(pvec_sd,[2,8]);

A_3sigma = intervalMatrix(m.A, 3*unc_stddev(:,1:2) );
B_3sigma = intervalMatrix(m.B, 3*unc_stddev(:,3)   );
M_3sigma = intervalMatrix([m.A m.B], 3*unc_stddev(:,1:3)   ); 

x_start = [-10; 10];
z_start = zonotope(x_start,[0.1 0 0.04; 0 0.1 -0.12]);

% Matrix zonotope identification
U_minus         = u(1:T-1);
Z_minus         = z(:,1:T-1);
Z_plus          = z(:,2:T);

Z_U_minus       = [Z_minus;
                   U_minus];
               
[~,S,~] = svd(Z_U_minus,'econ');
svd_ratio = S(3,3)/S(2,2)

% construct M_v - matrix zonotope of measurement noise
C_gam         = repmat(Z_gamma.center,1,T-1); %zeros(2,T-1);
G_gam         = cell(1,2*(T-1));
Gen_gam     = Z_gamma.generators;
for i = 1:T-1
    G_gam{i}            = zeros(2,T-1);
    G_gam{i}(:,i)       = Gen_gam(:,1); 
    G_gam{i+T-1}        = zeros(2,T-1);
    G_gam{i+T-1}(:,i)   = Gen_gam(:,2);
end
M_gamma         = matZonotope(C_gam,G_gam);

% construct M_w - matrix zonotope of process noise
C_w         = repmat(Z_w.center,1,T-1);
G_w         = cell(1,2*(T-1));
Gen_w       = Z_w.generators;
for i = 1:T-1
    G_w{i}            = zeros(2,T-1);
    G_w{i}(:,i)       = Gen_w(:,1); 
    G_w{i+T-1}        = zeros(2,T-1);
    G_w{i+T-1}(:,i)   = Gen_w(:,2);
end
M_w         = matZonotope(C_w,G_w);

% determine propogations matrices
M_dash      = (Z_plus - C_gam - C_w)*pinv(Z_U_minus);
M_AV        = Z_plus - M_dash*Z_U_minus + (-1)*M_gamma + (-1)*M_w;
Int_Mat_AV  = intervalMatrix(M_AV);
M_v_sup     = Int_Mat_AV.Sup;
M_v_inf     = Int_Mat_AV.Inf;
Z_max       = max(M_v_sup,[],2);
Z_min       = min(M_v_inf,[],2);
Z_AV        = zonotope(interval(Z_min, Z_max));

M_dash_A = M_dash(:,1:2);
M_dash_B = M_dash(:,3);

% determine M_Sigma using AV bound assumption
M_Sigma = (Z_plus + (-1)*M_gamma + A_true*M_gamma + (-1)*M_w)*pinv(Z_U_minus);

% plot the results
% figure(1);
% clf;
% hold on
% grid on
% 
% % z_zon       = M_dash*(z_start + Z_gamma) + Z_AV + Z_w;
% % z_zon_AV    = M_Sigma*(z_start) + Z_w;
% % 
% % z_zon       = reduce(z_zon,'girard',5);
% % z_zon_AV    = reduce(z_zon_AV,'girard',5);
% 
% plot(z_start,[1 2],'k-');
% 
% u_zon = cell(1,5);
% z_true = z_start;
% z_n4sid = z_start;
% z_zon = z_start;
% z_zon_AV = z_start;
% 
% for i = 1:2
%     
%     % Generate random input
%     u_zon{i} = zonotope(randPoint(Z_u),0.01);
%     
%     % True
%     z_true = [A_true B_true]*cartProd(z_true, u_zon{i}) + Z_w;
%     z_true = reduce(z_true,'girard',3);
%     plot(z_true,[1 2],'b-');
%     
%     % N4SID method
%     z_n4sid = M_3sigma*cartProd(z_n4sid, u_zon{i}) + Z_w;
%     z_n4sid = reduce(z_n4sid,'girard',3);
%     plot(z_n4sid,[1 2],'r*-');
%     
%     % Zonotope
%     z_zon = M_dash*cartProd((z_zon + Z_gamma), u_zon{i}) + Z_AV + Z_w;
%     z_zon = reduce(z_zon,'girard',3);
%     plot(z_zon,[1 2],'g*-');
%     
%     % AV Zonotope
%     plot(z_zon_AV,[1 2],'m+-');
% end
% 
% if useRandPointExtreme
%     samplemethod = "RandPointExtreme";
% else
%     samplemethod = "RandPoint";
% end
% 
% xlabel("state 1");
% ylabel("state 2");
% xlim([-15 -8]);
% ylim([2 12]);
% title({"Parameters",""," Noise bounds: c_w = " + c_w + ", c_{\gamma} = " + c_gam, "sampling: " + samplemethod});
% legend('Start','True','N4SID','M dash','M Sigma');

%% Run set-based observer

N   = 45;
nx  = 2;
nu  = 1;
q   = 3;

% define observability matrices
C       = cell(1,q);
C{1}    = [1 0.4];
C{2}    = [0.9 -1.2];
C{3}    = [-0.8 0.2;
           0    0.7];
% define bounds on v(k) using zonotopes
c_v_meas     = cell(1,q);
c_v_meas{1}  = 1.5;
c_v_meas{2}  = 1.5;
c_v_meas{3}  = 1.5;

Z_v_meas     = cell(1,q);
Z_v_meas{1}  = zonotope(0,c_v_meas{1});
Z_v_meas{2}  = zonotope(0,c_v_meas{2});
Z_v_meas{3}  = zonotope([0;0],blkdiag(c_v_meas{3},c_v_meas{3}));
% Z_v_meas{4}  = zonotope(2,c_v_meas{4});

x       = zeros(nx,N+1);
u       = zeros(nu,N+1);
x(:,1)  = [-9.9;-7]; %randPointExtreme(Z_X_0);   %           % true initial condition
w       = zeros(nx,N);
y       = cell(1,q);
v       = cell(1,q);

x_est_svd_z         = cell(1,N);
x_est_svd_prev_z    = cell(1,N);

x_est_svd_AV_z      = cell(1,N);
x_est_svd_AV_prev_z = cell(1,N);

x_est_svd_con_z     = cell(1,N);
x_est_svd_con_prev_z = cell(1,N);

x_est_opt_z         = cell(1,N);
x_est_opt_prev_z    = cell(1,N);

x_est_opt_con_z     = cell(1,N);
x_est_opt_con_prev_z = cell(1,N);

x_est_int_z         = cell(1,N);

for j = 1:q
    [C_sz,~]   = size(C{j});
    y{j}    = zeros(C_sz,N);
    v{j}    = zeros(C_sz,N);
end
t       = zeros(1,N);

% Kalman filter setup

% Kalman filter equations
%   [P,K,~] = idare(m.A',m.C',Q,R);
P = zeros(2,2,N+1);
P(:,:,1) = [12 1; 1 12];

x_hat(:,1)  = zeros(nx,1); % initial state estimate
% x_hat(:,1)  = x(:,1);
y_meas      = zeros(q,N);
x_hat_zon   = cell(1,N+1);

C_n4sid = [];
R = [];
for j = 1:q
   C_n4sid      = [C_n4sid;C{j}];
   [p_i,~]      = size(C{j});
   R = blkdiag(R,c_v_meas{j}*3*eye(p_i));
end
[sz,~]      = size(C_n4sid);
K           = zeros(nx,sz,N);
K(:,:,1)    = rand(nx,sz);
e_k         = cell(q,N);

% Simulation

% state estimate zonotopes - initial set estimate
x_est_svd       = Z_X_0;
x_est_svd_AV    = Z_X_0;
x_est_svd_con   = conZonotope(Z_X_0);
x_est_opt       = Z_X_0;
x_est_opt_con   = conZonotope(Z_X_0);

% zonotopes Z_{x|y^i}
Z_x_y_i         = cell(q,N);


Q = 45*c_w*eye(2);

[P_static,K,~] = idare(m.A',C_n4sid',Q,R);
K = K';

for k = 1:N
    %% Simulate true system
    disp("Timestep: " + t(k) )
    
    u(k)        = randPoint(Z_u);
    u_zon       = zonotope(u(k),0.01);
    
    w(:,k)      = randPoint(Z_w);
    x(:,k+1)    = A_true*x(:,k)  + w(:,k); % + B_true*u(k)
    for i = 1:q
        % v{i}(k) = randPoint(Z_v_meas{i});
        y{i}(:,k) = C{i}*x(:,k) + randPoint(Z_v_meas{i});
    end
    
    %% SYSID APPROACH
    
    disp(' - KF');
%     % Time-varying Kalman filter
%     R_ek = C_n4sid*P(:,:,k)*C_n4sid' + R;
%     K(:,:,k) = m.A*P(:,:,k)*C_n4sid'/R_ek;
%     P(:,:,k+1) = m.A*P(:,:,k)*m.A' + m.K*Q*m.K' - K(:,:,k)*R_ek*K(:,:,k)';
%     e_k_vec = [];
%     for i = 1:q
%         e_k_vec = [e_k_vec; y{i}(k) - C{i}*x_hat(:,k)];
%     end
%     x_hat(:,k+1) = m.A*x_hat(:,k) + m.B*u(k) + K(:,:,k)*e_k_vec;
%     x_hat_zon{k} = M_3sigma*cartProd(x_hat(:,k), u_zon);

    %Static KF
    
    e_k_vec = [];
    for i = 1:q
        e_k_vec = [e_k_vec; y{i}(:,k) - C{i}*x_hat(:,k)];
    end
    x_hat(:,k+1) = m.A*x_hat(:,k) + m.B*u(k) + K*e_k_vec;
    
    %% Advance time
    t(k+1)      = t(k) + Ts;
end

%% Plot the results

figure(2);
clf;
bd = 25;
xlim([-bd bd]);
ylim([-bd bd]);
grid on
hold on

for idx = (N-5):N
    % plot the true state and measurement regions
    plot(x(1,idx),x(2,idx),'k+','MarkerSize',20);
    
    % KALMAN FILTER 
    
%     [U,S,V] = svd(squeeze(P(:,:,idx)));
    [U,S,V] = svd(P_static);
    angle_vec = U(:,1);
    angle = atan2(angle_vec(2),angle_vec(1));
    ra = S(1,1);
    rb = S(2,2); 
    ang = angle; 
    x0  = x_hat(1,idx);
    y0  = x_hat(2,idx);
    ellipse(3*ra,3*rb,ang,x0,y0,'b-');
    plot(x_hat(1,idx),x_hat(2,idx),'b*');
    
    % ESTIMATE ZONOTOPES    
    % for j = 1:q
    %     plot(Z_x_y_i{j,idx},[1 2],'r-');
    % end
end
% legend('x(k)','SVD','SVD with AV assumption','SVD CON','OPT','OPT CON','N4SID 3 std dev');






















































%.
