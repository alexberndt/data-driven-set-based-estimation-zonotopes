%% Reachability Analysis comparison of N4SID and Matrix Zonotopes Methods
%
% We are considering the system 
%  x(k+1) = Ax(k) + Bu(k) + w(k)
%  z(k)   = x(k) + gamma(k)
%
% A = [0.9455   -0.2426;
%      0.2486    0.9455];
%
% B = [0.1; 0];
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

% Define initial condition
Z_X_0   = zonotope([0;0], blkdiag(15,15));
x(:,1)  = randPoint(Z_X_0);
x(:,1)  = [10.0;10.0];

% Define input zonotope
Z_u         = zonotope( 0, 10);

% Define noise zonotopes
c_gam       = 0.1;
c_w         = 0.1;
cw_offset   = 0;
cgam_offset = 0;
Z_gamma     = zonotope( [cgam_offset;cgam_offset], blkdiag(c_gam, c_gam) );
Z_w         = zonotope( [cw_offset;cw_offset], blkdiag(c_w, c_w) ); %[c_w 0 0.2*c_w; 0 0.9*c_w -0.05*c_w] ); 

% Generate training data
for k = 1:T
    % random input sequence
    u(k)        = 0.1*randPoint(Z_u);
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
               
[~,S,~] = svd(Z_U_minus,'econ')
svd_ratio = S(3,3)/S(2,2);

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
gcf = figure();
clf;
hold on
grid on

% z_zon       = M_dash*(z_start + Z_gamma) + Z_AV + Z_w;
% z_zon_AV    = M_Sigma*(z_start) + Z_w;
% 
% z_zon       = reduce(z_zon,'girard',5);
% z_zon_AV    = reduce(z_zon_AV,'girard',5);

plot(z_start,[1 2],'k-');

u_zon = cell(1,5);
z_true = z_start;
z_n4sid = z_start;
z_zon = z_start;
z_zon_AV = z_start;

for i = 1:3
    
    % Generate random input
    u_zon{i} = zonotope(1,0.01);
    
    % True
    z_true = [A_true B_true]*cartProd(z_true, u_zon{i}) + Z_w;
    z_true = reduce(z_true,'girard',3);
    plot(z_true,[1 2],'b-');
    
    % N4SID method
    z_n4sid = M_3sigma*cartProd(z_n4sid, u_zon{i}) + Z_w;
    z_n4sid = reduce(z_n4sid,'girard',3);
    plot(z_n4sid,[1 2],'r*-');
    
    % Zonotope
%     z_zon = M_dash*cartProd((z_zon + Z_gamma), u_zon{i}) + Z_AV + Z_w;
%     z_zon = reduce(z_zon,'girard',3);
%     plot(z_zon,[1 2],'m*-');
    
    % AV Zonotope
    z_zon_AV = M_Sigma*cartProd(z_zon_AV,u_zon{i}) + Z_w;
    z_zon_AV = reduce(z_zon_AV,'girard',3);
    plot(z_zon_AV,[1 2],'g+-');
end

if useRandPointExtreme
    samplemethod = "RandPointExtreme";
else
    samplemethod = "RandPoint";
end

xlabel("$x_1$",'Interpreter','latex');
ylabel("$x_2$",'Interpreter','latex');
xlim([-16 -5]);
ylim([-2 12]);
title("$c_w$ = " + c_w + ", $c_{\gamma}$ = " + c_gam + ", $T$ = " + T + ", bias $c_w$ = " + cw_offset + ", bias $c_\gamma$ = " + cgam_offset,'Interpreter','latex'); %+ ", ExtremeSampling = " + num2str(useRandPointExtreme)
title("$c_w$ = " + c_w + ", $c_{\gamma}$ = " + c_gam + ", $T$ = " + T + ", svd ratio = " + svd_ratio,'Interpreter','latex'); %+ ", ExtremeSampling = " + num2str(useRandPointExtreme)
legend('Start','True','N4SID','Mat Zon'); %,'M dash','M Sigma');


save_loc = '/home/alberndt/Documents/research/data_driven/berndt2020zonotope_analysis/figures/';
fig_name = strrep(strcat('res_cw_',num2str(c_w),'_cgam_',num2str(c_gam),'_T_',num2str(T),'_ExtremeSample_',num2str(useRandPointExtreme)),'.','_');
fig_name = strrep(strcat('res_cw_',num2str(c_w),'_cgam_',num2str(c_gam),'_T_',num2str(T),'_sample_comp'),'.','_');
fig_name = strrep(strcat('res_cw_',num2str(c_w),'_cgam_',num2str(c_gam),'_T_',num2str(T),'_offset_cw_',num2str(cw_offset),'_offset_cgam_',num2str(cgam_offset)),'.','_');
fig_name = strrep(strcat('res_cw_',num2str(c_w),'_cgam_',num2str(c_gam),'_T_',num2str(T),'_svd_',num2str(svd_ratio)),'.','_');
fig_name = strcat(fig_name, '.eps');

% saveas(gcf,strcat(save_loc,fig_name), 'epsc')

