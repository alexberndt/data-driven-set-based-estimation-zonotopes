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
cd '/home/alberndt/Documents/research/data_driven/code/data_driven_set_based_estimation_zonotopes'
rng("default");
addpath('./functions/');

% Training phase

% Run the 'training' phase
A_true  = [0.9455   -0.2426;
           0.2486    0.9455];
B_true  = [0.1; 0];
T       = 5000;
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

input_ratios = logspace(-1,2,20);

bounds = struct();
bounds.true.sup = zeros(1,length(input_ratios));
bounds.true.inf = zeros(1,length(input_ratios));
bounds.n4sid.sup = zeros(1,length(input_ratios));
bounds.n4sid.inf = zeros(1,length(input_ratios));
bounds.zono.sup = zeros(1,length(input_ratios));
bounds.zono.inf = zeros(1,length(input_ratios));

svd_ratios = zeros(1,length(input_ratios));

idx = 1;
for u_multi = input_ratios

    % Generate training data
    for k = 1:T
        % random input sequence
        u(k)        = u_multi*randPoint(Z_u);
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

    x_start = [0; 0];
    z_start = zonotope(x_start,[0.1 0 0.04; 0 0.1 -0.12]);

    % Matrix zonotope identification
    U_minus         = u(1:T-1);
    Z_minus         = z(:,1:T-1);
    Z_plus          = z(:,2:T);

    Z_U_minus       = [Z_minus;
                       U_minus];

    [~,S,~] = svd(Z_U_minus,'econ');
    svd_ratio = S(3,3)/S(2,2);
    svd_ratios(idx) = svd_ratio;

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

    

    % z_zon       = M_dash*(z_start + Z_gamma) + Z_AV + Z_w;
    % z_zon_AV    = M_Sigma*(z_start) + Z_w;
    % 
    % z_zon       = reduce(z_zon,'girard',5);
    % z_zon_AV    = reduce(z_zon_AV,'girard',5);

%     plot(z_start,[1 2],'k-');

    u_zon = cell(1,5);
    z_true = z_start;
    z_n4sid = z_start;
    z_zon = z_start;
    z_zon_AV = z_start;

    for i = 1:1

        % Generate random input
        u_zon{i} = zonotope(1,0.01);

        % True
        z_true = [A_true B_true]*cartProd(z_true, u_zon{i}) + Z_w;
%         z_true = reduce(z_true,'girard',3);
%         plot(z_true,[1 2],'b-');

        % N4SID method
        z_n4sid = M_3sigma*cartProd(z_n4sid, u_zon{i}) + Z_w;
%         z_n4sid = reduce(z_n4sid,'girard',3);
%         plot(z_n4sid,[1 2],'r*-');

        % Zonotope
    %     z_zon = M_dash*cartProd((z_zon + Z_gamma), u_zon{i}) + Z_AV + Z_w;
    %     z_zon = reduce(z_zon,'girard',3);
    %     plot(z_zon,[1 2],'m*-');

        % AV Zonotope
        z_zon_AV = M_Sigma*cartProd(z_zon_AV,u_zon{i}) + Z_w;
%         z_zon_AV = reduce(z_zon_AV,'girard',3);
%         plot(z_zon_AV,[1 2],'g+-');
    end
    
    i_true      = interval(z_true);
    i_n4sid     = interval(z_n4sid);
    i_zon_AV    = interval(z_zon_AV);
    
    bounds.true.sup(idx)    = i_true.sup(1);
    bounds.true.inf(idx)    = i_true.inf(1);
    bounds.n4sid.sup(idx)   = i_n4sid.sup(1);
    bounds.n4sid.inf(idx)   = i_n4sid.inf(1);
    bounds.zono.sup(idx)    = i_zon_AV.sup(1);
    bounds.zono.inf(idx)    = i_zon_AV.inf(1);
    idx = idx + 1
end

% if useRandPointExtreme
%     samplemethod = "RandPointExtreme";
% else
%     samplemethod = "RandPoint";
% end

% plot the results
gcf = figure();
clf;

h1 = semilogx(svd_ratios,bounds.true.sup,'k-');
hold on
h2 = semilogx(svd_ratios,bounds.true.inf,'k-');

h3 = semilogx(svd_ratios,bounds.n4sid.sup,'b-');
h4 = semilogx(svd_ratios,bounds.n4sid.inf,'b-');

h5 = semilogx(svd_ratios,bounds.zono.sup,'m-');
h6 = semilogx(svd_ratios,bounds.zono.inf,'m-');

grid on
% set(gcf, 'XScale', 'log');

xlabel("singular value ratio = (last singular value)/(second-to-last singular value)",'Interpreter','latex');
ylabel("min-max bounds of $\mathcal{X}_{k+1}$",'Interpreter','latex');
% xlim([-16 -5]);
% ylim([-2 12]);
% title("$c_w$ = " + c_w + ", $c_{\gamma}$ = " + c_gam + ", $T$ = " + T + ", bias $c_w$ = " + cw_offset + ", bias $c_\gamma$ = " + cgam_offset,'Interpreter','latex'); %+ ", ExtremeSampling = " + num2str(useRandPointExtreme)
% title("$c_w$ = " + c_w + ", $c_{\gamma}$ = " + c_gam + ", $T$ = " + T + ", svd ratio = " + svd_ratio,'Interpreter','latex'); %+ ", ExtremeSampling = " + num2str(useRandPointExtreme)

title("$T$ = " + T,'Interpreter','latex');

hleglines = [h1(1) h3(1) h5(1)];
legend(hleglines,'True','N4SID','Mat Zon'); %,'M dash','M Sigma');

x0 = 10;
y0 = 10;
width = 450;
height = 200;

set(gcf,'units','points','position',[x0,y0,width,height])

save_loc = '/home/alberndt/Documents/research/data_driven/berndt2020zonotope_analysis/figures/';
fig_name = strrep(strcat('res_svd_T_',num2str(T),'_svdratiorange_',num2str(min(svd_ratios)),'_to_',num2str(max(svd_ratios))),'.','_');
fig_name = strcat(fig_name, '.eps');

saveas(gcf,strcat(save_loc,fig_name), 'epsc')
