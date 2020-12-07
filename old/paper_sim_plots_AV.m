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
rng("default");
addpath('./functions/');

% Run the 'training' phase
A_true  = [0.9455   -0.2426;
           0.2486    0.9455];
% B_true  = [1; 0];
T       = 350;
Ts      = 1.0;

useRandPointExtreme = false;

% init matrices
x       = zeros(2,T+1);
z       = zeros(2,0.2*T);
t       = zeros(1,T+1);
gam     = zeros(2,T);
w       = zeros(2,T);
u       = zeros(1,T);

% define initial condition
Z_X_0  = zonotope([0;0], blkdiag(15,15));
x(:,1)  = randPoint(Z_X_0);
x(:,1)  = [10.0;10.0];

% define noise zonotopes
c_gam       = 0.01;
c_w         = 0.01;
Z_gamma     = zonotope( [0.0;0.0], blkdiag(c_gam, c_gam) );
Z_w         = zonotope( [0.0;0.0] ,blkdiag(c_w, c_w) ); % [c_w 0 0.2*c_w; 0 0.9*c_w -0.05*c_w] ); 

% generate training data
for k = 1:T
    % random input sequence
    % u(k)        = 0.1*rand;
    % random bounded noise
    if useRandPointExtreme
        gam(:,k)    = randPointExtreme(Z_gamma);
        w(:,k)      = randPointExtreme(Z_w);
    else
        gam(:,k)    = randPoint(Z_gamma);
        w(:,k)      = randPoint(Z_w);
    end
    % system evolution
    x(:,k+1)    = A_true*x(:,k)  + w(:,k); % + B*u(k); 
    z(:,k)      = x(:,k) + gam(:,k);
    t(k+1)      = t(k) + Ts;
end

% N4SID model
sysid_data = iddata(z',[],Ts);
m = n4sid(sysid_data,2,'ssp','can'); % canonical form for C = eye(2)
[pvec,pvec_sd] = getpvec(m);
unc_stddev = reshape(pvec_sd,[2,6]);
M_n4sid = intervalMatrix(m.A,3*unc_stddev(:,1:2));

x_start = [-10; 10];
z_start = zonotope(x_start,[0.1 0 0.04; 0 0.1 -0.12]);
% z_start = zonotope([-10; 10],0.1*eye(2));
z_n4sid = M_n4sid*z_start + Z_w; 
z_true = A_true*z_start + Z_w;

% Matrix zonotope identification
Z_minus         = z(:,1:T-1);
Z_plus          = z(:,2:T);

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
M_dash      = (Z_plus - C_gam - C_w)*pinv(Z_minus);
M_AV        = Z_plus - M_dash*Z_minus + (-1)*M_gamma + (-1)*M_w;
Int_Mat_AV  = intervalMatrix(M_AV);
M_v_sup     = Int_Mat_AV.Sup;
M_v_inf     = Int_Mat_AV.Inf;
Z_max       = max(M_v_sup,[],2);
Z_min       = min(M_v_inf,[],2);
Z_AV        = zonotope(interval(Z_min, Z_max));

M_Sigma = (Z_plus + (-1)*M_gamma + A_true*M_gamma + (-1)*M_w)*pinv(Z_minus);

%% Run set-based observer

N   = 15;
nx  = 2;
q   = 3;

% define observability matrices
C       = cell(1,q);
C{1}    = [1 0];
C{2}    = [1 -1];
C{3}    = [1 0;
           0 1];
% C{4}    = [1 1];
% define bounds on v(k) using zonotopes
c_v_meas     = cell(1,q);
c_v_meas{1}  = 1.6;
c_v_meas{2}  = 1.7;
c_v_meas{3}  = 2.8;
% c_v_meas{4}  = 10.8;
Z_v_meas     = cell(1,q);
Z_v_meas{1}  = zonotope(0,c_v_meas{1});
Z_v_meas{2}  = zonotope(0,c_v_meas{2});
Z_v_meas{3}  = zonotope([0;0],blkdiag(c_v_meas{3},c_v_meas{3}));
% Z_v_meas{4}  = zonotope(0,c_v_meas{4});
% 
x       = zeros(nx,N+1);
x(:,1)  = [-9.9;-7]; %randPointExtreme(Z_X_0);   %           % true initial condition
w       = zeros(nx,N);
y       = cell(1,q);
v       = cell(1,q);

x_est_svd_z         = cell(1,N);
x_est_svd_prev_z    = cell(1,N);

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

%% Kalman filter setup

Q = 15*c_w*eye(2);

% Kalman filter equations
%   [P,K,~] = idare(m.A',m.C',Q,R);
P = zeros(2,2,N+1);
P(:,:,1) = [10 1; 1 10];

x_hat(:,1)  = zeros(nx,1); % initial state estimate
y_meas      = zeros(q,N);
x_hat_zon   = cell(1,N+1);

C_n4sid = [];
R = [];
for j = 1:q
   C_n4sid      = [C_n4sid;C{j}];
   [p_i,~]      = size(C{j});
   R = blkdiag(R,c_v_meas{j}*9*eye(p_i));
end
[sz,~]      = size(C_n4sid);
K           = zeros(nx,sz,N);
K(:,:,1)    = rand(nx,sz);
e_k         = cell(q,N);

% Simulation

% state estimate zonotopes - initial set estimate
x_est_svd   = Z_X_0;
x_est_svd_con = conZonotope(Z_X_0);
x_est_opt   = Z_X_0;
x_est_opt_con = conZonotope(Z_X_0);

% zonotopes Z_{x|y^i}
Z_x_y_i   = cell(q,N);

for k = 1:N
    %% Simulate true system
    disp("Timestep: " + t(k) )
    w(:,k)      = randPoint(Z_w);
    x(:,k+1)    = A_true*x(:,k) + w(:,k);
    for i = 1:q
        % v{i}(k) = randPoint(Z_v_meas{i});
        y{i}(:,k) = C{i}*x(:,k) + randPoint(Z_v_meas{i});
    end
    
    %% ZONOTOPE SVD METHOD
    disp(' - Zonotope SVD method');
    x_est_svd_prev_z{k} = x_est_svd;
    x_est_svd_prev      = x_est_svd;
    
    for i = 1:q
        Z_x_y_i{i,k}   = measurement_zonotope(y{i}(:,k), C{i}, Z_v_meas{i});
        meas_zonotope  = Z_x_y_i{i,k};
        x_est_svd   = and(Z_x_y_i{i,k},x_est_svd);
        plot(meas_zonotope,[1 2],'r-');
    end
    x_est_svd_z{k}  = x_est_svd;

    % propogate current state estimate
    x_est_svd_kp1   = M_Sigma * (x_est_svd ) + Z_w;
    x_est_svd       = reduce(x_est_svd_kp1,'girard',5);

    %% CONZONOTOPE SVD METHOD
    disp(' - ConZonotope SVD method');
    x_est_svd_con_prev_z{k} = x_est_svd_con;
    x_est_svd_con_prev      = x_est_svd_con;
    
    for i = 1:q
        conZono_measurement = conZonotope(Z_x_y_i{i,k});
        x_est_svd_con   = and(conZono_measurement,x_est_svd_con);
    end
    x_est_svd_con_z{k}  = x_est_svd_con;
    x_est_svd_con_kp1   = M_dash * (x_est_svd_con + Z_gamma ) + Z_AV + Z_w;    
    x_est_svd_con       = reduce(x_est_svd_con_kp1,'girard',5);
    
    plot(x_est_svd_con,[1 2],'k--');
    
    %% ZONOTOPE OPT METHOD
    disp(' - Zonotope OPT method');
    x_est_opt_prev_z{k} = x_est_opt;
    x_est_opt_prev      = x_est_opt;
    
    yl = {y{1}(:,k), y{2}(:,k), y{3}(:,k) };    
    x_est_opt = intersectZonoZono(x_est_opt_prev,C,Z_v_meas,yl,'frobenius');
    x_est_opt_z{k}  = x_est_opt;

    % propogate current state estimate 
    x_est_opt_kp1   = M_Sigma * (x_est_opt) + Z_w;
    x_est_opt       = reduce(x_est_opt_kp1,'girard',5);

    %% CONZONOTOPE OPT METHOD
    disp(' - ConZonotope OPT method');
    x_est_opt_con_prev_z{k} = x_est_opt_con;
    x_est_opt_con_prev      = x_est_opt_con;

    yl = {y{1}(:,k), y{2}(:,k), y{3}(:,k) };
%     plot(x_est_opt_con,[1 2],'g*-');
    x_est_opt_con = intersectConZonoZono2(x_est_opt_con, C, Z_v_meas, yl, 'frobenius');
%     plot(x_est_opt_con,[1 2],'b--');
    x_est_opt_con_z{k}  = x_est_opt_con;
    
    % propogate current state estimate 
    x_est_opt_con_kp1   = M_dash * (x_est_opt_con + Z_gamma ) + Z_AV + Z_w; 
    % x_est_opt_con       = x_est_opt_con_kp1;
    x_est_opt_con       = reduce(x_est_opt_con_kp1,'girard',5);

    %% INTERSECTION OF BOTH SVD and OPT METHOD
    
%     figure(5);
%     clf;
%    
%     hold on;
%     plot(x_est_svd_con,[1 2],'r*-');
%     plot(x_est_opt_con,[1 2],'b-');
% 
%     x_est_inter = and( x_est_svd_con, x_est_opt_con );
% 
%     x_est_int_z{k} = x_est_inter;
%     
%     plot(x_est_inter,[1 2],'k--');
    
    %% SYSID APPROACH
    
    disp(' - KF');
    % Kalman filter
    R_ek = C_n4sid*P(:,:,k)*C_n4sid' + R;
    K(:,:,k) = m.A*P(:,:,k)*C_n4sid'/R_ek;
    P(:,:,k+1) = m.A*P(:,:,k)*m.A' + m.K*Q*m.K' - K(:,:,k)*R_ek*K(:,:,k)';
    e_k_vec = [];
    for i = 1:q
        e_k_vec = [e_k_vec; y{i}(k) - C{i}*x_hat(:,k)];
    end
    x_hat(:,k+1) = m.A*x_hat(:,k) + K(:,:,k)*e_k_vec;
    x_hat_zon{k} = M_n4sid*x_hat(:,k);
    
    %% Advance time
    t(k+1)      = t(k) + Ts;
    
%     pause(1.0);
end

%% FIGURE 

width = 600;
height = 250;

%% Plot the results

gcf2 = figure(2);
set(gcf2, 'Position',  [100, 600, width, height])
clf;
% subplot(1,2,1);
bd = 25;
xlim([4.6 11.9]);
ylim([-0.1 8.1]);
grid on
hold on

for idx = 13
    % plot the true state and measurement regions
    plot(x(1,idx),x(2,idx),'k+','MarkerSize',20);
    
    % ZONTOPE

    % METHOD 1
    % plot(x_est_svd_prev_z{idx},[1 2],'k--'); 
    plot(x_est_svd_z{idx},[1 2],'k*-');
    plot(x_est_svd_con_z{idx},[1 2],'m*-');  
    
    % METHOD 2
    % plot the state-estimate zonotope 
    plot(x_est_opt_z{idx},[1 2],'g--'); 
    plot(x_est_opt_con_z{idx},[1 2],'b--');
    
%     conZon1 = x_est_svd_con_z{idx};
%     conZon2 = x_est_opt_con_z{idx};
%     x_inter = and(conZon1, conZon2);
    % plot(x_est_int_z{idx},[1 2],'k*--');
%     plot(x_inter,[1 2],'k*--');

    % KALMAN FILTER 
    
    [U,S,V] = svd(squeeze(P(:,:,idx)));
    angle_vec = U(:,1);
    angle = atan2(angle_vec(2),angle_vec(1));
    ra = S(1,1);
    rb = S(2,2); 
    ang = angle; 
    x0  = x_hat(1,idx);
    y0  = x_hat(2,idx);
    ellipse(3*ra,3*rb,ang,x0,y0,'b-');
    plot(x_hat(1,idx),x_hat(2,idx),'b*');
    
end
legend('$x(k)$','Method 1','Method 1 CON','Method 2','Method 2 CON','N4SID 3 std dev','Interpreter','latex','Location','southeast');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');


saveas(gcf2,'figures/sim_2.eps','epsc');

%%
gcf3 = figure(3);
set(gcf3, 'Position',  [100, 250, width, height])
clf;
xlim([1.9 8.3]);
ylim([2.3 10.2]);
grid on
hold on

for idx = 15
    % plot the true state and measurement regions
    plot(x(1,idx),x(2,idx),'k+','MarkerSize',20);
    
    % ZONTOPE

    % METHOD 1
    % plot(x_est_svd_prev_z{idx},[1 2],'k--'); 
    plot(x_est_svd_z{idx},[1 2],'k*-');
    plot(x_est_svd_con_z{idx},[1 2],'m*-');  
    
    % METHOD 2
    % plot the state-estimate zonotope 
    plot(x_est_opt_z{idx},[1 2],'g--'); 
    plot(x_est_opt_con_z{idx},[1 2],'b--');
    
%     conZon1 = x_est_svd_con_z{idx};
%     conZon2 = x_est_opt_con_z{idx};
%     x_inter = and(conZon1, conZon2);
    % plot(x_est_int_z{idx},[1 2],'k*--');
%     plot(x_inter,[1 2],'k*--');

    % KALMAN FILTER 
    
    [U,S,V] = svd(squeeze(P(:,:,idx)));
    angle_vec = U(:,1);
    angle = atan2(angle_vec(2),angle_vec(1));
    ra = S(1,1);
    rb = S(2,2); 
    ang = angle; 
    x0  = x_hat(1,idx);
    y0  = x_hat(2,idx);
    ellipse(3*ra,3*rb,ang,x0,y0,'b-');
    plot(x_hat(1,idx),x_hat(2,idx),'b*');
    
end
legend('$x(k)$','Method 1','Method 1 CON','Method 2','Method 2 CON','N4SID 3 std dev','Interpreter','latex','Location','southeast');
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');

saveas(gcf3,'figures/sim_3.eps','epsc');

%% Measure radius 

% x_est_svd_radius        = zeros(1,N);
% x_est_opt_radius        = zeros(1,N);
% x_est_svd_con_radius    = zeros(1,N);
% x_est_opt_con_radius    = zeros(1,N);
% x_n4sid_max_radius      = zeros(1,N);
% 
% for idx = 1:N
%     
%     zon_to_tes = x_est_svd_z{idx};
%     zon_to_test = zonotope([0;0], zon_to_tes.generators);
%     
% %     x_est_svd_radius(idx)       = norm(x_est_svd_z{idx});
%     x_est_svd_radius(idx)       = radius(zon_to_test);
%     x_est_opt_radius(idx)       = norm(x_est_opt_z{idx});
% %     x_est_svd_con_radius(idx)   = volume(x_est_svd_con_z{idx});
% %     x_est_opt_con_radius(idx)   = volume(x_est_opt_con_z{idx});
% 
%     [U,S,V] = svd(squeeze(P(:,:,idx)));
%     ra = S(1,1);
%     rb = S(2,2); 
%     x_n4sid_max_radius(idx) = max(ra,rb);
%     
% end
% 
% figure(10);
% clf;
% hold on
% grid on
% stairs(x_est_svd_radius,'k-');
% % stairs(x_est_opt_radius,'r');
% stairs(x_n4sid_max_radius,'b');
% ylim([-1 15]);

%%

% figure(6);
% grid on
% hold on
% 
% plot(x_est_svd_con_z{3},[1 2],'m*-');
% plot(x_est_opt_con_z{3},[1 2],'r*-');
% 
% conZon1 = x_est_svd_con_z{3};
% conZon2 = x_est_opt_con_z{3};
% 
% x_inter = and(conZon1, conZon2);
% 
% plot(x_inter,[1 2],'b*-');































































%.
