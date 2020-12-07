%% Rank condition test
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
cd '/home/alberndt/Documents/research/data_driven/code/data_driven_set_based_estimation_zonotopes'
rng("default");
addpath('./functions/');

% Training phase

% Run the 'training' phase
A_true  = [0.9455   -0.2426;
           0.2486    0.9455];
B_true  = [0.1; 0];
T       = 50;
Ts      = 1.0;

useRandPointExtreme = false;

% Init constant trajectories
gam_const   = zeros(2,T);
w_const     = zeros(2,T);
u_const     = zeros(1,T);


% Define input zonotope
Z_u         = zonotope( 0, 10);

% Define noise zonotopes
c_gam       = 0.02;
c_w         = 0.02;
Z_gamma     = zonotope( [0.0;0.0], blkdiag(c_gam, c_gam) );
Z_w         = zonotope( [0.0;0.0], [c_w 0 0.2*c_w; 0 0.9*c_w -0.05*c_w] ); 

%% GENERATE A-PRIORI DATA

for k = 1:T
    u_const(k)        = randPoint(Z_u);
    if useRandPointExtreme
        gam_const(:,k)    = randPointExtreme(Z_gamma);
        w_const(:,k)      = randPointExtreme(Z_w);
    else
        gam_const(:,k)    = randPoint(Z_gamma);
        w_const(:,k)      = randPoint(Z_w);
    end
end

%% A-PRIORI ESTIMATOR DATA

N   = 3;
nx  = 2;
nu  = 1;
q   = 3;

Z_X_0   = zonotope([0;0], blkdiag(15,15));

% define observability matrices
C       = cell(1,q);
C{1}    = [1 0.4];
C{2}    = [0.9 -1.2];
C{3}    = [-0.8 0.2;
           0    0.7];
       
% define bounds on v(k) using zonotopes
c_v_meas     = cell(1,q);
c_v_meas{1}  = 15.5;
c_v_meas{2}  = 15.7;
c_v_meas{3}  = 15.1;

Z_v_meas     = cell(1,q);
Z_v_meas{1}  = zonotope(0,c_v_meas{1});
Z_v_meas{2}  = zonotope(0,c_v_meas{2});
Z_v_meas{3}  = zonotope([0;0],blkdiag(c_v_meas{3},c_v_meas{3}));

% 
x_est       = zeros(nx,N+1);
u_est       = zeros(nu,N+1);
x_est(:,1)  = [-9.9;-7]; %randPointExtreme(Z_X_0);   %           % true initial condition
w_est       = zeros(nx,N);
y_est       = cell(1,q);
v_est       = cell(1,q);
t_est       = zeros(1,N+1);

for j = 1:q
    [C_sz,~]    = size(C{j});
    y_est{j}    = zeros(C_sz,N);
    v_est{j}    = zeros(C_sz,N);
end


for k_est = 1:N
    
    % generate measurement noise signal
    for j = 1:q
        v_est{j}(:,k_est) = randPoint(Z_v_meas{j});
    end
    
    % generate input sequence
    u_est(k_est) = randPoint(Z_u);
    w_est(:,k_est) = randPoint(Z_w);
   
end


%%

u_bounds = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0];
u_bounds = logspace(-2,0.1,20);
% u_bounds = [0.01, 0.5, 2.0, 5.0];

u_bnds_len = length(u_bounds);
svd_ratios = zeros(1,u_bnds_len);
results = struct();
results.true = cell(1,u_bnds_len);
results.svd = cell(1,u_bnds_len);

u_idx = 1;

for u_bnd = u_bounds
    
    %% IDENTIFICATION STEP
    disp("u bound: " + num2str(u_bnd));
    
    x_est(:,1)  = [-9.9;-8.9];
    
    % init vectors
    t           = zeros(1,T+1);
    x_learn     = zeros(2,T+1);
    z_learn     = zeros(2,T);
    u_learn     = zeros(1,T);
    
    % Define initial condition
    % Z_X_0   = zonotope([0;0], blkdiag(15,15));
    % x(:,1)  = randPoint(Z_X_0);
    x_learn(:,1)  = [-10.0;10.0];
    
    % Generate training data
    for k = 1:T
        % random input sequence
        u_learn(k) = u_bnd*u_const(k);

        % system evolution
        x_learn(:,k+1)    = A_true*x_learn(:,k) + B_true*u_learn(k) + w_const(:,k);
        z_learn(:,k)      = x_learn(:,k) + gam_const(:,k);
        t(k+1)      = t(k) + Ts; 
    end
    
    % N4SID model
    sysid_data = iddata(z_learn',u_learn',Ts);
    m = n4sid(sysid_data,2,'ssp','can'); % canonical form for C = eye(2)
    [pvec,pvec_sd] = getpvec(m);

    unc_stddev = reshape(pvec_sd,[2,8]);

    A_3sigma = intervalMatrix(m.A, 3*unc_stddev(:,1:2) );
    B_3sigma = intervalMatrix(m.B, 3*unc_stddev(:,3)   );
    M_3sigma = intervalMatrix([m.A m.B], 3*unc_stddev(:,1:3)   ); 

    x_start = [-10; 10];
    z_start = zonotope(x_start,[0.1 0 0.04; 0 0.1 -0.12]);

    % Matrix zonotope identification
    U_minus         = u_learn(1:T-1);
    Z_minus         = z_learn(:,1:T-1);
    Z_plus          = z_learn(:,2:T);

    Z_U_minus       = [Z_minus;
                       U_minus];

    [~,S,~] = svd(Z_U_minus,'econ');
    svd_ratio = S(3,3)/S(2,2);
    svd_ratios(u_idx) = svd_ratio;
    disp("  - svd_ratio: " + num2str(svd_ratio));

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
    
    
    %% ESTIMATION STEP
    % state estimate zonotopes - initial set estimate
    x_est_svd       = Z_X_0;
    x_est_svd_AV    = Z_X_0;
    x_est_svd_con   = conZonotope(Z_X_0);
    x_est_opt       = Z_X_0;
    x_est_opt_con   = conZonotope(Z_X_0);
    
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

    % zonotopes Z_{x|y^i}
    Z_x_y_i         = cell(q,N);
    
    for k = 1:N
        %% Simulate true system
        disp("Timestep: " + t_est(k) )

        x_est(:,k+1)    = A_true*x_est(:,k) + B_true*u_est(k) + w_est(:,k);
        for i = 1:q
            % v{i}(k) = randPoint(Z_v_meas{i});
            y_est{i}(:,k) = C{i}*x_est(:,k) + v_est{i}(:,k);
        end
        
        yl = {y_est{1}(:,k), y_est{2}(:,k), y_est{3}(:,k) };   
        
        %% ZONOTOPE SVD METHOD
        disp(' - Zonotope SVD method');
        x_est_svd_prev_z{k} = x_est_svd;
        x_est_svd_prev      = x_est_svd;

        % meas_list = cell(1,q+1);
        for i = 1:q
            Z_x_y_i{i,k}   = measurement_zonotope(y_est{i}(:,k), C{i}, Z_v_meas{i});
            meas_zonotope  = Z_x_y_i{i,k};
            % meas_list{i}   = Z_x_y_i{i,k};
            x_est_svd = and(meas_zonotope, x_est_svd);
        end
        % meas_list{q+1} = x_est_svd;
        % x_est_svd   = andAveraging1(meas_list);

        x_est_svd_z{k}  = x_est_svd;

        % propogate current state estimate
        
        u_zon = zonotope(u_est(k),0.01);

        x_est_svd_kp1   = M_Sigma*cartProd(x_est_svd,u_zon) + Z_w;
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

        % x_est_svd_con_kp1   = M_dash * (x_est_svd_con + Z_gamma) + Z_AV + Z_w;
        x_est_svd_con_kp1   = M_Sigma*cartProd(x_est_svd_con,u_zon) + Z_w;
        x_est_svd_con       = reduce(x_est_svd_con_kp1,'girard',5);
        
        %% ZONOTOPE OPT METHOD
        disp(' - Zonotope OPT method');
        x_est_opt_prev_z{k} = x_est_opt;
        x_est_opt_prev      = x_est_opt;
 
        x_est_opt = intersectZonoZono(x_est_opt_prev,C,Z_v_meas,yl,'frobenius');
        x_est_opt_z{k}  = x_est_opt;

        % propogate current state estimate 
        % x_est_opt_kp1   = M_dash * (x_est_opt + Z_gamma) + Z_AV + Z_w;
        x_est_opt_kp1   = M_Sigma*cartProd(x_est_opt,u_zon) + Z_w;
        x_est_opt       = reduce(x_est_opt_kp1,'girard',5);
        
        %% CONZONOTOPE OPT METHOD
        disp(' - ConZonotope OPT method');
        x_est_opt_con_prev_z{k} = x_est_opt_con;
        x_est_opt_con_prev      = x_est_opt_con;

        if isempty(x_est_opt_con)
            disp("no zonotope found in last step.");
            x_est_opt_con = x_est_svd_con;
        end

        x_est_opt_con = intersectConZonoZono2(x_est_opt_con, C, Z_v_meas,yl,'frobenius');
        x_est_opt_con_z{k}  = x_est_opt_con;

        % propogate current state estimate 
        % x_est_opt_con_kp1   = M_dash * (x_est_opt_con + Z_gamma) + Z_AV + Z_w;
        x_est_opt_con_kp1   = M_Sigma*cartProd(x_est_opt_con,u_zon)  + Z_w;

        % x_est_opt_con       = x_est_opt_con_kp1;
        x_est_opt_con       = reduce(x_est_opt_con_kp1,'girard',5);
        
        %% Advance time
        t_est(k+1)      = t_est(k) + Ts;
        
    end
    
    results.svd{u_idx} = x_est_svd_z{N};    
    results.opt{u_idx} = x_est_opt_z{N}; 
    results.svdc{u_idx} = x_est_svd_con_z{N}; 
    results.optc{u_idx} = x_est_opt_con_z{N}; 

    resuts.true{u_idx} = x_est(:,N);
    
    if u_idx == 5
       figure(9);
       hold on
       plot(x_est(1,k),x_est(2,k),'k*');
       plot(x_est_svd_z{N},[1 2],'b-');
       plot(x_est_opt_z{N},[1 2],'g-');
       plot(x_est_svd_con_z{N},[1 2],'r-');
       plot(x_est_opt_con_z{N},[1 2],'m-');
    end
    
    u_idx = u_idx + 1;
    
end
   
%% Plot the bounds
coord = 1; % or 2 for x_1 or x_2 comparison

bounds = struct();
len_u_bnds = length(u_bounds);
bounds.true      = zeros(1,len_u_bnds);
bounds.svd.sup   = zeros(1,len_u_bnds);
bounds.svd.inf   = zeros(1,len_u_bnds);
bounds.svdc.sup  = zeros(1,len_u_bnds);
bounds.svdc.inf  = zeros(1,len_u_bnds);
bounds.opt.sup   = zeros(1,len_u_bnds);
bounds.opt.inf   = zeros(1,len_u_bnds);
bounds.optc.sup  = zeros(1,len_u_bnds);
bounds.optc.inf  = zeros(1,len_u_bnds);

u_idx = 1;
for u_bnd = u_bounds
    % plot the true state and measurement regions
%     plot(x(1,idx),x(2,idx),'k+','MarkerSize',20);
    
    bounds.true(u_idx) = resuts.true{u_idx}(coord);
    
    % ZONOTOPE

    % METHOD 1    
    int_svd     = interval(results.svd{u_idx});
    bounds.svd.sup(u_idx) = int_svd.sup(coord);
    bounds.svd.inf(u_idx) = int_svd.inf(coord);
    
    int_svdc     = interval(results.svdc{u_idx});
    bounds.svdc.sup(u_idx) = int_svdc.sup(coord);
    bounds.svdc.inf(u_idx) = int_svdc.inf(coord);
    
    % METHOD 2    
    int_opt     = interval(results.opt{u_idx});
    bounds.opt.sup(u_idx) = int_opt.sup(coord);
    bounds.opt.inf(u_idx) = int_opt.inf(coord);
    
    int_optc    = interval(results.optc{u_idx});
    bounds.optc.sup(u_idx) = int_optc.sup(coord);
    bounds.optc.inf(u_idx) = int_optc.inf(coord);
  
%     % KALMAN FILTER 
%     [U,S,V] = svd(squeeze(P(:,:,idx)));
%     angle_vec = U(:,1);
%     angle = atan2(angle_vec(2),angle_vec(1));
%     ra = S(1,1);
%     rb = S(2,2); 
%     ang = angle; 
%     x0  = x_hat(1,idx);
%     y0  = x_hat(2,idx);
    
    u_idx = u_idx + 1;
end 

figure(1);
clf;
h1 = semilogx(svd_ratios, bounds.true, 'k*-');
hold on

h2  = semilogx(svd_ratios, bounds.svd.sup, 'b+-');
h2b = semilogx(svd_ratios, bounds.svd.inf, 'b+-');

h3  = semilogx(svd_ratios, bounds.svdc.sup, 'r+-');
h3b = semilogx(svd_ratios, bounds.svdc.inf, 'r+-');

h4  = semilogx(svd_ratios, bounds.opt.sup, 'g+-');
h4b = semilogx(svd_ratios, bounds.opt.inf, 'g+-');

h5  = semilogx(svd_ratios, bounds.optc.sup, 'm+-');
h5b = semilogx(svd_ratios, bounds.optc.inf, 'm+-');

xlabel("singular value ratio");
grid on
legend([h1,h2,h3,h4,h5],'true','Method 1','Method 1 CON','Method 2','Method 2 CON');




















































































%.