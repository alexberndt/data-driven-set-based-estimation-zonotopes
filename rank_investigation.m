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

ctrb(A_true,B_true)

useRandPointExtreme = false;

% Init matrices
x       = zeros(2,T+1);
z       = zeros(2,T);
t       = zeros(1,T+1);
gam     = zeros(2,T);
w       = zeros(2,T);
u       = zeros(1,T);
u_const = zeros(1,T);

% Define initial condition
Z_X_0   = zonotope([0;0], blkdiag(15,15));
x(:,1)  = randPoint(Z_X_0);
x(:,1)  = [10.0;10.0];

% Define input zonotope
Z_u         = zonotope( 0, 10);

% Define noise zonotopes
c_gam       = 0.1;
c_w         = 0.1;
Z_gamma     = zonotope( [0.0;0.0], blkdiag(c_gam, c_gam) );
Z_w         = zonotope( [0.0;0.0], [c_w 0 0.2*c_w; 0 0.9*c_w -0.05*c_w] ); 

% generate random sequences a-priori
for k = 1:T
    u_const(k) = randPoint(Z_u);
    if useRandPointExtreme
        gam(:,k)    = randPointExtreme(Z_gamma);
        w(:,k)      = randPointExtreme(Z_w);
    else
        gam(:,k)    = randPoint(Z_gamma);
        w(:,k)      = randPoint(Z_w);
    end
end

%%

u_bound_vals = [0.01, 0.05, 0.1, 0.15, 0.2, 0.5, 1.0, 2.0];
% u_bound_vals = [0.1, 0.2, 1.0];

results = struct();
results.svd = cell(1,size(u_bound_vals,2));
u_idx = 1;

for u_bnd = u_bound_vals
    disp("u_bnd: " + num2str(u_bnd))

    % Init matrices
    x       = zeros(2,T+1);
    z       = zeros(2,T);
    x(:,1)  = [10.0;10.0];
    u       = zeros(1,T);

    % Generate training data
    for k = 1:T
        % random input sequence
        u(k)        = u_bnd*u_const(k);
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
    rank(Z_U_minus)

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

end

%%
for u_bnd = u_bound_vals
    
    %% Run set-based observer

    N   = 4;
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
    c_v_meas{1}  = 8.5;
    c_v_meas{2}  = 8.5;
    c_v_meas{3}  = 8.5;

    Z_v_meas     = cell(1,q);
    Z_v_meas{1}  = zonotope(0,c_v_meas{1});
    Z_v_meas{2}  = zonotope(0,c_v_meas{2});
    Z_v_meas{3}  = zonotope([0;0],blkdiag(c_v_meas{3},c_v_meas{3}));
    % Z_v_meas{4}  = zonotope(2,c_v_meas{4});

    x_est       = zeros(nx,N+1);
    u_est       = zeros(nu,N+1);
    x_est(:,1)  = [-9.9;-7]; %randPointExtreme(Z_X_0);   %           % true initial condition
    w_est       = zeros(nx,N);
    y_est       = cell(1,q);
    v_est       = cell(1,q);

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
        y_est{j}    = zeros(C_sz,N);
        v_est{j}    = zeros(C_sz,N);
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
       R = blkdiag(R,c_v_meas{j}*15*eye(p_i));
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

    for k = 1:N
        %% Simulate true system
        disp("Timestep: " + t(k) )

        u_est(k)        = 2.5; %randPoint(Z_u);
        u_zon       = zonotope(u_est(k),0.01);

        w_est(:,k)      = randPoint(Z_w);
        x_est(:,k+1)    = A_true*x_est(:,k) + B_true*u_est(k) + w_est(:,k);
        for i = 1:q
            % v{i}(k) = randPoint(Z_v_meas{i});
            y_est{i}(:,k) = C{i}*x_est(:,k) + randPoint(Z_v_meas{i});
        end

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

        x_est_svd_kp1   = M_Sigma*cartProd(x_est_svd,u_zon) + Z_w;
        x_est_svd       = reduce(x_est_svd_kp1,'girard',5);

        %% ZONOTOPE SVD METHOD - NO AV ASSUMPION
    %     disp(' - Zonotope SVD (no AV assumption) method');
    %     x_est_svd_AV_prev_z{k} = x_est_svd_AV;
    %     x_est_svd_AV_prev      = x_est_svd_AV;
    %     
    %     for i = 1:q
    %         Z_x_y_i{i,k}   = measurement_zonotope(y{i}(:,k), C{i}, Z_v_meas{i});
    %         meas_zonotope  = Z_x_y_i{i,k};
    %         x_est_svd_AV   = and(Z_x_y_i{i,k},x_est_svd_AV);
    %     end
    %     x_est_svd_AV_z{k}  = x_est_svd_AV;
    % 
    %     % propogate current state estimate
    %     x_est_svd_AV_kp1   = M_dash * (x_est_svd_AV + Z_gamma) + Z_AV + Z_w;
    %     x_est_svd_AV       = reduce(x_est_svd_AV_kp1,'girard',5);

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

        yl = {y_est{1}(:,k), y_est{2}(:,k), y_est{3}(:,k) };   
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

        yl = {y_est{1}(:,k), y_est{2}(:,k), y_est{3}(:,k) }; 
        x_est_opt_con = intersectConZonoZono2(x_est_opt_con, C, Z_v_meas,yl,'frobenius');
        x_est_opt_con_z{k}  = x_est_opt_con;

        % propogate current state estimate 
        % x_est_opt_con_kp1   = M_dash * (x_est_opt_con + Z_gamma) + Z_AV + Z_w;
        x_est_opt_con_kp1   = M_Sigma*cartProd(x_est_opt_con,u_zon)  + Z_w;

        % x_est_opt_con       = x_est_opt_con_kp1;
        x_est_opt_con       = reduce(x_est_opt_con_kp1,'girard',5);

        %% SYSID APPROACH

        disp(' - KF');
        %Static KF
        e_k_vec = [];
        for i = 1:q
            e_k_vec = [e_k_vec; y_est{i}(:,k) - C{i}*x_hat(:,k)];
        end
        x_hat(:,k+1) = m.A*x_hat(:,k) + m.B*u_est(k) + K*e_k_vec;

        %% Advance time
        t(k+1)      = t(k) + Ts;

    end


    results.svd{u_idx} = x_est_svd;
    results.svdc{u_idx} = x_est_svd_con;

    results.opt{u_idx} = x_est_opt;
    results.optc{u_idx} = x_est_opt_con;

    u_idx = u_idx + 1;

end

%%
u_idx = 1;
coord = 1;

figure(1);
hold on
grid on

for u_bnd = u_bound_vals
    
    int_svd     = interval(results.svd{u_idx});
    bounds.svd.sup(u_idx) = int_svd.sup(coord);
    bounds.svd.inf(u_idx) = int_svd.inf(coord);
       
    int_opt     = interval(results.opt{u_idx});
    bounds.opt.sup(u_idx) = int_opt.sup(coord);
    bounds.opt.inf(u_idx) = int_opt.inf(coord);
    
    u_idx = u_idx + 1;
end

plot(u_bound_vals, bounds.svd.sup, 'm-');
plot(u_bound_vals, bounds.svd.inf, 'm-');

plot(u_bound_vals, bounds.opt.sup, 'g-');
plot(u_bound_vals, bounds.opt.inf, 'g-');


























































%.
