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
% author:  Alexander Berndt, Amr Alanwar
% contact: alberndt@kth.se, amr.alanwar@gmail.com
%
%% Training phase to identify the system
clc
clear
% cd '/home/alberndt/Documents/research/data_driven/code/data_driven_set_based_estimation_zonotopes'
% rng(1);
 addpath('functions/');

% Training phase

% Run the 'training' phase
A_true  = [0.9455   -0.2426;
           0.2486    0.9455];
B_true  = [0.1; 0];
T       = 100;
Ts      = 1.0;
N   = 50;
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
c_v_meas{1}  = 1.0;
c_v_meas{2}  = 1.0;
c_v_meas{3}  = 1.0;

Z_v_meas     = cell(1,q);
Z_v_meas{1}  = zonotope(0,c_v_meas{1});
Z_v_meas{2}  = zonotope(0,c_v_meas{2});
Z_v_meas{3}  = zonotope([0;0],diag([c_v_meas{3},c_v_meas{3}]));

useRandPointExtreme = false;

% Init matrices
x       = zeros(2,T+1);
%z       = zeros(2,T);
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
c_gam       = 0.02;
c_w         = 0.02;
Z_gamma     = zonotope( [0.0;0.0], blkdiag(c_gam, c_gam) );
Z_w         = zonotope( [0.0;0.0], blkdiag(c_w, c_w) ); % [c_w 0 0.2*c_w; 0 0.9*c_w -0.05*c_w] ); 

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
    for i = 1:q
        % v{i}(k) = randPoint(Z_v_meas{i});
        z{i}(:,k) = C{i}*x(:,k) + randPoint(Z_v_meas{i});
    end
    z_ssyid(:,k)      = x(:,k) + gam(:,k);
    t(k+1)      = t(k) + Ts;
end

% N4SID model
sysid_data = iddata(z_ssyid',u',Ts);
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
for i = 1:q
Z_minus{i}         = z{i}(:,1:T-1);
Z_plus{i}          = z{i}(:,2:T);
Z_U_minus{i}       = [Z_minus{i};
                   U_minus];
end

               

% construct M_v - matrix zonotope of measurement noise
for j=1:q
    C_gam         = repmat(Z_v_meas{j}.center,1,T-1); %zeros(2,T-1);
    G_gam = {};
    Gen_gam     = Z_v_meas{j}.generators;
    for i = 1:T-1
        G_gam{i}            = zeros(length(Z_v_meas{j}.center),T-1);
        G_gam{i}(:,i)       = Gen_gam(:,1);
        if size(Gen_gam,2)==2
        G_gam{i+T-1}        = zeros(length(Z_v_meas{j}.center),T-1);
        G_gam{i+T-1}(:,i)   = Gen_gam(:,2);
        end
    end
    M_gamma{j}         = matZonotope(C_gam,G_gam);
end
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


for i = 1:q
Z_plus_zono{i} = measurement_zonotope_matzono(Z_plus{i},C{i},M_gamma{i});
end

for i = 1:q
Z_minus_zono{i} = measurement_zonotope_matzono(Z_minus{i},C{i},M_gamma{i});
end

selectedNode =1;
XU_cent = [Z_minus_zono{selectedNode}.center; U_minus];
for i =1:length(Z_minus_zono{selectedNode}.generator)
   XU_gen{i} =  [Z_minus_zono{selectedNode}.generator{i};zeros(size(U_minus))];
end
XU_matzono = intervalMatrix(matZonotope(XU_cent,XU_gen));
M_Sigma = intervalMatrix((Z_plus_zono{selectedNode} + (-1)*M_w))*invIntMatrixAutoRinv(XU_matzono);

intAB1 = M_Sigma.int;
intAB1.sup >= [A_true,B_true]
intAB1.inf <= [A_true,B_true]
% % plot the results
% % figure(1);
% % clf;
% % hold on
% % grid on
% % 
% % % z_zon       = M_dash*(z_start + Z_gamma) + Z_AV + Z_w;
% % % z_zon_AV    = M_Sigma*(z_start) + Z_w;
% % % 
% % % z_zon       = reduce(z_zon,'girard',5);
% % % z_zon_AV    = reduce(z_zon_AV,'girard',5);
% % 
% % plot(z_start,[1 2],'k-');
% % 
% % u_zon = cell(1,5);
% % z_true = z_start;
% % z_n4sid = z_start;
% % z_zon = z_start;
% % z_zon_AV = z_start;
% % 
% % for i = 1:2
% %     
% %     % Generate random input
% %     u_zon{i} = zonotope(randPoint(Z_u),0.01);
% %     
% %     % True
% %     z_true = [A_true B_true]*cartProd(z_true, u_zon{i}) + Z_w;
% %     z_true = reduce(z_true,'girard',3);
% %     plot(z_true,[1 2],'b-');
% %     
% %     % N4SID method
% %     z_n4sid = M_3sigma*cartProd(z_n4sid, u_zon{i}) + Z_w;
% %     z_n4sid = reduce(z_n4sid,'girard',3);
% %     plot(z_n4sid,[1 2],'r*-');
% %     
% %     % Zonotope
% %     z_zon = M_Sigma*cartProd(z_zon, u_zon{i}) +  Z_w;
% %     z_zon = reduce(z_zon,'girard',3);
% %     plot(z_zon,[1 2],'g*-');
% %     
% %     % AV Zonotope
% %     plot(z_zon_AV,[1 2],'m+-');
% % end
% % 
% % if useRandPointExtreme
% %     samplemethod = "RandPointExtreme";
% % else
% %     samplemethod = "RandPoint";
% % end
% % 
% % xlabel("state 1");
% % ylabel("state 2");
% % xlim([-15 -8]);
% % ylim([2 12]);
% % title({"Parameters",""," Noise bounds: c_w = " + c_w + ", c_{\gamma} = " + c_gam, "sampling: " + samplemethod});
% % legend('Start','True','N4SID','M dash','M Sigma');

%% Run set-based observer


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

Q = 5*c_w*eye(2);

% Kalman filter equations
%   [P,K,~] = idare(m.A',m.C',Q,R);
P = zeros(2,2,N+1);
P(:,:,1) = [2 1; 1 2];

x_hat(:,1)  = zeros(nx,1); % initial state estimate
% x_hat(:,1)  = x(:,1);
y_meas      = zeros(q,N);
x_hat_zon   = cell(1,N+1);

C_n4sid = [];
R = [];
for j = 1:q
   C_n4sid      = [C_n4sid;C{j}];
   [p_i,~]      = size(C{j});
   R = blkdiag(R,c_v_meas{j}*eye(p_i));
end
[sz,~]      = size(C_n4sid);
K           = zeros(nx,sz,N);
K(:,:,1)    = rand(nx,sz);
e_k         = cell(q,N);

[P_static,K,~] = idare(m.A',C_n4sid',Q,R);
K = K';

% Simulation

% state estimate zonotopes - initial set estimate
x_est_svd       = Z_X_0;
x_est_svd_AV    = Z_X_0;
x_est_svd_con   = conZonotope(Z_X_0);
x_est_opt       = Z_X_0;
x_est_opt_con   = conZonotope(Z_X_0);

% zonotopes Z_{x|y^i}
Z_x_y_i         = cell(q,N);

comp_times = struct();
comp_times.svd = zeros(1,N);
comp_times.opt = zeros(1,N);
comp_times.svdc = zeros(1,N);
comp_times.optc = zeros(1,N);

for k = 1:N
    %% Simulate true system
    disp("Timestep: " + t(k) )
    
    u(k)        = randPoint(Z_u);
    u_zon       = zonotope(u(k),0.01);
    
    w(:,k)      = randPoint(Z_w);
    x(:,k+1)    = A_true*x(:,k) + B_true*u(k) + w(:,k);
    for i = 1:q
        % v{i}(k) = randPoint(Z_v_meas{i});
        y{i}(:,k) = C{i}*x(:,k) + randPoint(Z_v_meas{i});
    end
    
    %% ZONOTOPE SVD METHOD
    disp(' - Zonotope SVD method');
    
    x_est_svd_prev_z{k} = x_est_svd;
    x_est_svd_prev      = x_est_svd;
    
    % meas_list = cell(1,q+1);
    tic;
    for i = 1:q
        Z_x_y_i{i,k}   = measurement_zonotope(y{i}(:,k), C{i}, Z_v_meas{i});
        meas_zonotope  = Z_x_y_i{i,k};
        % meas_list{i}   = Z_x_y_i{i,k};
        x_est_svd = and(meas_zonotope, x_est_svd);
    end
    % meas_list{q+1} = x_est_svd;
    % x_est_svd   = andAveraging1(meas_list);
    
    x_est_svd_z{k}  = x_est_svd;

    % propogate current state estimate
    x_est_svd_kp1   = M_Sigma*cartProd(x_est_svd,u_zon) + Z_w;
    x_est_svd       = reduce(x_est_svd_kp1,'girard',5);
    comp_times.svd(k) = toc;
    


    %% CONZONOTOPE SVD METHOD
    disp(' - ConZonotope SVD method');
    tic;
    x_est_svd_con_prev_z{k} = x_est_svd_con;
    x_est_svd_con_prev      = x_est_svd_con;
    
    for i = 1:q
        conZono_measurement = conZonotope(Z_x_y_i{i,k});
        x_est_svd_con   = and(conZono_measurement,x_est_svd_con);
    end
    x_est_svd_con_z{k}  = x_est_svd_con;

    % x_est_svd_con_kp1   = M_dash * (x_est_svd_con + Z_gamma) + Z_AV + Z_w;
    x_est_svd_con_kp1   = M_Sigma*zonotope(cartProd(x_est_svd_con,u_zon)) + Z_w;
    x_est_svd_con       = reduce(x_est_svd_con_kp1,'girard',5);
    comp_times.svdc(k) = toc;
    
    %% ZONOTOPE OPT METHOD
    disp(' - Zonotope OPT method');
    tic;
    x_est_opt_prev_z{k} = x_est_opt;
    x_est_opt_prev      = x_est_opt;
    
    yl = {y{1}(:,k), y{2}(:,k), y{3}(:,k) };    
    x_est_opt = intersectZonoZono(x_est_opt_prev,C,Z_v_meas,yl,'frobenius');
    x_est_opt_z{k}  = x_est_opt;

    % propogate current state estimate 
    % x_est_opt_kp1   = M_dash * (x_est_opt + Z_gamma) + Z_AV + Z_w;
    x_est_opt_kp1   = M_Sigma*cartProd(x_est_opt,u_zon) + Z_w;
    x_est_opt       = reduce(x_est_opt_kp1,'girard',5);
    comp_times.opt(k) = toc;

    %% CONZONOTOPE OPT METHOD
    disp(' - ConZonotope OPT method');
    tic;
    x_est_opt_con_prev_z{k} = x_est_opt_con;
    x_est_opt_con_prev      = x_est_opt_con;
    
    if isempty(x_est_opt_con)
        disp("no zonotope found in last step.");
        x_est_opt_con = x_est_svd_con;
    end

    yl = {y{1}(:,k), y{2}(:,k), y{3}(:,k) }; 
    x_est_opt_con = intersectConZonoZono2(conZonotope(x_est_opt_con), C, Z_v_meas,yl,'frobenius');
    x_est_opt_con_z{k}  = x_est_opt_con;
    
    % propogate current state estimate 
    % x_est_opt_con_kp1   = M_dash * (x_est_opt_con + Z_gamma) + Z_AV + Z_w;
    x_est_opt_con_kp1   = M_Sigma*zonotope(cartProd(x_est_opt_con,u_zon))  + Z_w;
    
    % x_est_opt_con       = x_est_opt_con_kp1;
    x_est_opt_con       = reduce(x_est_opt_con_kp1,'girard',5);
    comp_times.optc(k) = toc;

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

mean(comp_times.svd)
mean(comp_times.svdc)
mean(comp_times.opt)
mean(comp_times.optc)

% gcf1 = figure(1);
% clf;
% grid on 
% hold on 
% plot(comp_times.svd);
% plot(comp_times.opt);
% plot(comp_times.svdc);
% plot(comp_times.optc);
% legend('Approach 1','Approach 2','Approach 1 CZ','Approach 2 CZ','Interpreter','latex');

%%

gcf2 = figure('Renderer', 'painters', 'Position', [10 10 700 700]);
clf;
bd = 25;
% xlim([-bd bd]);
% ylim([-bd bd]);
grid on
hold on

for idx = 34
    % plot the true state and measurement regions
    plot(x(1,idx),x(2,idx),'k+','MarkerSize',20);
    
    % ZONTOPE

    % METHOD 1
    % plot(x_est_svd_prev_z{idx},[1 2],'k--');
    plot(x_est_svd_z{idx},[1 2],'k*-');
    % plot(x_est_svd_AV_z{idx},[1 2],'r-');
    plot(x_est_svd_con_z{idx},[1 2],'m*-');  
    
    % METHOD 2
    % plot the state-estimate zonotope 
    plot(x_est_opt_z{idx},[1 2],'g--'); 
    plot(x_est_opt_con_z{idx},[1 2],'b--');
    
    % INTERSECTION OF METHOD 1 AND 2
    % conZon1 = x_est_svd_con_z{idx};
    % conZon2 = x_est_opt_con_z{idx};
    % x_inter = and(conZon1, conZon2);
    % plot(x_est_int_z{idx},[1 2],'k*--');
    % plot(x_inter,[1 2],'k*--');

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
  %  plot(x_hat(1,idx),x_hat(2,idx),'b*');
    
    % ESTIMATE ZONOTOPES    
    % for j = 1:q
    %     plot(Z_x_y_i{j,idx},[1 2],'r-');
    % end
end
x0 = 500;
y0 = 300;
width = 450;
height = 150;

%xlim([3.5 7]);
%ylim([-5.5, -2]);
% xlim([-10 10]);
% ylim([-14 0]);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
warOrig = warning; warning('off','all');
%set(gcf2,'units','points','position',[x0,y0,width,height])

legend('$x(k)$','Approach 1','Approach 1 CZ','Approach 2','Approach 2 CZ','KF $3\sigma$ bounds','Interpreter','latex');
   box on
warning(warOrig);
    ax = gca;
    ax.FontSize = 18;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
%% Measure radius 
% 
coord = 1; % or 2 for x_1 or x_2 comparison

bounds = struct();
bounds.true      = zeros(1,N);
bounds.n4sid.sup = zeros(1,N);
bounds.n4sid.inf = zeros(1,N);
bounds.svd.sup   = zeros(1,N);
bounds.svd.inf   = zeros(1,N);
bounds.svdc.sup  = zeros(1,N);
bounds.svdc.inf  = zeros(1,N);
bounds.opt.sup   = zeros(1,N);
bounds.opt.inf   = zeros(1,N);
bounds.optc.sup  = zeros(1,N);
bounds.optc.inf  = zeros(1,N);

for idx = 1:N
    % plot the true state and measurement regions
    % plot(x(1,idx),x(2,idx),'k+','MarkerSize',20);
    
    bounds.true(idx) = x(coord,idx);
    
    % ZONOTOPE

    % METHOD 1    
    %int_svd     = interval(x_est_svd_prev_z{idx});
    int_svd     = interval(x_est_svd_z{idx});
    bounds.svd.sup(idx) = int_svd.sup(coord);
    bounds.svd.inf(idx) = int_svd.inf(coord);
    
    %int_svdc    = interval(x_est_svd_con_prev_z{idx});
    int_svdc    = interval(x_est_svd_con_z{idx});
    bounds.svdc.sup(idx) = int_svdc.sup(coord);
    bounds.svdc.inf(idx) = int_svdc.inf(coord);
    
    % METHOD 2    
    %int_opt     = interval(x_est_opt_prev_z{idx});
    int_opt     = interval(x_est_opt_z{idx});
    bounds.opt.sup(idx) = int_opt.sup(coord);
    bounds.opt.inf(idx) = int_opt.inf(coord);
    
    %int_optc    = interval(x_est_opt_con_prev_z{idx});
    int_optc    = interval(x_est_opt_con_z{idx});
    bounds.optc.sup(idx) = int_optc.sup(coord);
    bounds.optc.inf(idx) = int_optc.inf(coord);
  
    % KALMAN FILTER N4SID 3
    [U,S,V] = svd(squeeze(P(:,:,idx)));
    angle_vec = U(:,1);
    angle = atan2(angle_vec(2),angle_vec(1));
    ra = S(1,1);
    rb = S(2,2); 
    ang = angle; 
    x0  = x_hat(1,idx);
    y0  = x_hat(2,idx);
end 
%%
gcf6 = figure('Renderer', 'painters', 'Position', [10 10 700 300]);;
clf;
hold on
grid on

h1 = plot(t(1:N),bounds.svdc.sup,'b-');
h2 = plot(t(1:N),bounds.svdc.inf,'b-');

% h3 = plot(t(1:N),bounds.svd.sup,'g-');
% h4 = plot(t(1:N),bounds.svd.inf,'g-');

h5 = plot(t(1:N),bounds.optc.sup,'m--');
h6 = plot(t(1:N),bounds.optc.inf,'m--');

% h7 = plot(t(1:N),bounds.opt.sup,'r-');
% h8 = plot(t(1:N),bounds.opt.inf,'r-');

h9 = plot(t(1:N),bounds.true,'k-.','Linewidth',2);

xlabel("time step $k$",'Interpreter','latex');
ylabel("$x_1$",'Interpreter','latex');

hlegend = [h9,h1,h5];
x0 = 10;
y0 = 10;
width = 450;
height = 150;

%ylim([-15 15]);

%set(gcf6,'units','points','position',[x0,y0,width,height])
legend(hlegend,'$x(k)$','Approach 1 CZ','Approach 2 CZ','Location','southeast','Interpreter','latex');
   box on
warning(warOrig);
    ax = gca;
    ax.FontSize = 18;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

%%
gcf7 = figure('Renderer', 'painters', 'Position', [10 10 700 300]);;
clf;
hold on
grid on

% h1 = plot(t(1:N),bounds.svdc.sup,'b-');
% h2 = plot(t(1:N),bounds.svdc.inf,'b-');

h3 = plot(t(1:N),bounds.svd.sup,'b-');
h4 = plot(t(1:N),bounds.svd.inf,'b-');

% h5 = plot(t(1:N),bounds.optc.sup,'m--');
% h6 = plot(t(1:N),bounds.optc.inf,'m--');

h7 = plot(t(1:N),bounds.opt.sup,'m--');
h8 = plot(t(1:N),bounds.opt.inf,'m--');

h9 = plot(t(1:N),bounds.true,'k-.','Linewidth',2);

xlabel("time step $k$",'Interpreter','latex');
ylabel("$x_1$",'Interpreter','latex');

hlegend = [h9,h3,h7];
x0 = 10;
y0 = 10;
width = 450;
height = 150;

%ylim([-14.8 15.2]);

%set(gcf7,'units','points','position',[x0,y0,width,height])
legend(hlegend,'$x(k)$','Approach 1','Approach 2','Location','southeast','Interpreter','latex');
   box on
warning(warOrig);
    ax = gca;
    ax.FontSize = 18;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

