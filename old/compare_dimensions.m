% Larger state space


%% Training phase to identify the system
clc
clear
rng("default");
addpath('./functions/');

% Run the 'training' phase
% A_true  = [0.9455   -0.2426;
%            0.2486    0.9455];       
% A_true = [-1 -4  0  0  0; 
%           4  -1  0  0  0; 
%           0   0 -3  1  0; 
%           0   0 -1 -3  0; 
%           0   0  0  0 -2];
      
A_true = [0.9323 -0.189 0 0 0;
          0.1890 0.9323 0 0 0;
          0        0    0.8596 0.0430 0;
          0        0    -0.0430  0.8596 0; 
          0        0     0    0 0.9048];
% B_true = ones(5,1);
nx = 5;

T       = 500;
Ts      = 1.0;

useRandPointExtreme = true;

% init matrices
x       = zeros(nx,T+1);
z       = zeros(nx,T);
t       = zeros(1,T+1);
gam     = zeros(nx,T);
w       = zeros(nx,T);
u       = zeros(1,T);


% Z_X_0  = zonotope([0;0], blkdiag(10,10));
% x(:,1)  = randPoint(Z_X_0);
x(:,1)  = -1*ones(5,1);

% define noise zonotopes
c_gam       = 0.0001;
c_w         = 0.0001;
Z_gamma     = zonotope( zeros(nx,1), c_gam*eye(nx) );
Z_w         = zonotope( zeros(nx,1), c_w*eye(nx)  ); 

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
m = n4sid(sysid_data,5,'ssp','can'); % canonical form for C = eye(2)
[pvec,pvec_sd] = getpvec(m);
unc_stddev = reshape(pvec_sd,[nx,3*nx]);
M_n4sid = intervalMatrix(m.A,3*unc_stddev(:,1:nx));

% IS the true A contained within the N4SID interval matrix?
A_true > M_n4sid.Inf;
A_true < M_n4sid.Sup;

%% Zonotope identification
Z_minus         = z(:,1:T-1);
Z_plus          = z(:,2:T);

% construct M_v - matrix zonotope of measurement noise
C_gam         = repmat(Z_gamma.center,1,T-1); %zeros(2,T-1);
G_gam         = cell(1,nx*(T-1));
Gen_gam     = Z_gamma.generators;
for i = 1:T-1
    for j = 1:nx
        G_gam{i+(T-1)*(j-1)}            = zeros(nx,T-1);
        G_gam{i+(T-1)*(j-1)}(:,i)       = Gen_gam(:,j); 
%         G_gam{i+T-1}        = zeros(2,T-1);
%         G_gam{i+T-1}(:,i)   = Gen_gam(:,2);
    end
end
M_gam         = matZonotope(C_gam,G_gam);

% construct M_w - matrix zonotope of process noise
C_w         = repmat(Z_w.center,1,T-1); %zeros(2,T-1);
G_w         = cell(1,nx*(T-1));
Gen_w       = Z_w.generators;
for i = 1:T-1
    for j = 1:nx
        G_w{i+(T-1)*(j-1)}            = zeros(nx,T-1);
        G_w{i+(T-1)*(j-1)}(:,i)       = Gen_w(:,j); 
%         G_gam{i+T-1}        = zeros(2,T-1);
%         G_gam{i+T-1}(:,i)   = Gen_gam(:,2);
    end
end
M_w         = matZonotope(C_w,G_w);

% construct M_w - matrix zonotope of process noise
% C_w         = repmat(Z_w.center,1,T-1);
% G_w         = cell(1,nx*(T-1));
% Gen_w       = Z_w.generators;
% for i = 1:T-1
%     G_w{i}            = zeros(2,T-1);
%     G_w{i}(:,i)       = Gen_w(:,1); 
%     G_w{i+T-1}        = zeros(2,T-1);
%     G_w{i+T-1}(:,i)   = Gen_w(:,2);
% end
% M_w         = matZonotope(C_w,G_w);

% determine propogations matrices
M_dash      = (Z_plus - C_gam - C_w)*pinv(Z_minus);
M_AV        = Z_plus - M_dash*Z_minus + (-1)*M_w + (-1)*M_gam;
Int_Mat_AV  = intervalMatrix(M_AV);
M_v_sup     = Int_Mat_AV.Sup;
M_v_inf     = Int_Mat_AV.Inf;
Z_max       = max(M_v_sup,[],2);
Z_min       = min(M_v_inf,[],2);
Z_AV        = zonotope(interval(Z_min, Z_max));

%% REACHABILITY
x_start = [-1; -1; 0; 1; 0];
z_start = zonotope(x_start,[0.1*eye(nx) [0.1;0.05;0.05;0.05;0.01]]);
% % z_start = zonotope([-10; 10],0.1*eye(2));
z_n4sid = M_n4sid*z_start + Z_w; 
z_true  = A_true*z_start + Z_w;
z_zon = M_dash*(z_start + Z_gamma) + Z_AV + Z_w;

figure(1);
clf;
hold on
grid on
plot(z_start,[1 2],'k-');
plot(z_true,[1 2],'b-');
plot(z_n4sid,[1 2],'r*-');
plot(z_zon,[1 2],'g*-');
for i = 1:1
    % True 
    plot(z_true,[1 2],'b-');
    z_true = A_true*z_true + Z_w;
%     z_true = reduce(z_true,'girard',5);
    % with some points
%     pont_1 = A_true*pont_1 + Z_w;
%     plot(pont_1,[1 2],'k-');
%     pont_2 = A_true*pont_2 + Z_w;
%     plot(pont_2,[1 2],'m-');
%     pont_3 = A_true*pont_3 + Z_w;
%     plot(pont_3,[1 2],'g-');
%     pont_4 = A_true*pont_4 + Z_w;
%     plot(pont_4,[1 2],'y-');
    % N4SID method
    plot(z_n4sid,[1 2],'r*-');
    z_n4sid = M_n4sid*z_n4sid + Z_w;
%     z_n4sid = reduce(z_n4sid,'girard',5);
    % Zonotope
    plot(z_zon,[1 2],'g*-');
    z_zon = M_dash*(z_zon + Z_gamma) + Z_AV + Z_w;
%     z_zon = reduce(z_zon,'girard',5);
end

legend('start','true','n4sid','zon');

z_true_rad   = radius(z_true)
z_zon_rad    = radius(z_zon)
z_n4sid_rad  = radius(z_n4sid)
Z_gamma_rad  = radius(Z_gamma)
Z_w_rad      = radius(Z_w)







































