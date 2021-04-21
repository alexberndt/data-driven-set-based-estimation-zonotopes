%% Illustration of SVD proof
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
%% SVD matrix 
clc
clear
cd '/home/alberndt/Documents/research/data_driven/code/data_driven_set_based_estimation_zonotopes'
rng(1);
addpath('./functions/');

y = 80;
C = [0.4 -3.8];
Z_v = zonotope(100.0, 10.0);

[p,n] = size(C);
[U,~,V]     = svd(C);
V1              = V(:,1:p);
V2              = V(:,p+1:n);
U1              = U;

c_k = [20;10];
G_k = [0    2.1    0.1;
       1   0   2  ];
R_k = zonotope(c_k, G_k);

K = 300;

r_R_k = radius(R_k);

% M = radius(V2'*R_k) + abs(V2'*c_k);
% M = radius(R_k) + norm(c_k);

M = r_R_k + norm(V2'*c_k);
% M = K;

Z_x_y = measurement_zonotope_M(y,C,Z_v,M);

figure(1);
clf;
plot(Z_x_y, [1 2], 'k-');
hold on
plot(R_k, [1 2], 'b-');
h1 = ellipse(K,K,2*pi,0,0,'r');
h2 = ellipse(r_R_k,r_R_k,2*pi,c_k(1),c_k(2),'g-');
legend("$Z_{x|y}$","$R_k$","$K$ bound","radius($R_k$)","Interpreter","latex");

%% Long skinny C
clc
clear
cd '/home/alberndt/Documents/research/data_driven/code/data_driven_set_based_estimation_zonotopes'
rng(1);
addpath('./functions/');

% % wide short C
% y =  4 ;
% C = [ 1    0.3 ];
% Z_v = zonotope(1, 10);
% % 
% % skinny tall C
% y = [0;1;3];
% C = [ 1    0.7 ;
%       0     0.6  ;
%       0.3  0.1];
% Z_v = zonotope([0;0;0], 10*eye(3));

% rank-deficient C
y = [0;0;0];
C = [ 6    0 ;
      6    0
      3    0 ];
Z_v = zonotope([0;0;0], 10*eye(3));
% % 
% % square invertible C
% y = [0;1;3];
% C = [ 1    0.7  0.1;
%       0    0.6  0.2;
%       0.3  0    0.9];
% Z_v = zonotope([0;0;0], 10*eye(3));

% [p,n] = size(C);
% [U,Sigma,V]     = svd(C);
% V1              = V(:,1:p);
% V2              = V(:,p+1:n);
% U1              = U;

c_k = [20;10];
G_k = [0    2.1    0.1;
       1   0   2  ];
R_k = zonotope(c_k, G_k);
K = 300;
r_R_k = radius(R_k);

% M = radius(V2'*R_k) + abs(V2'*c_k);
% M = radius(R_k) + norm(c_k);

% M = r_R_k + norm(V2'*c_k);
M = K;

Z_x_y = measurement_zonotope_M(y,C,Z_v,M);

figure(1);
clf;
plot(Z_x_y, [1 2], 'k-');
hold on
plot(R_k, [1 2], 'b-');
h1 = ellipse(K,K,2*pi,0,0,'r');
h2 = ellipse(r_R_k,r_R_k,2*pi,c_k(1),c_k(2),'g-');
legend("$Z_{x|y}$","$R_k$","$K$ bound","radius($R_k$)","Interpreter","latex");

%% 3D case

C = [1 0.4 0.7 ;
     0 0 0.2 ];
Z_v = zonotope([-5.0; 10.0], [1 0 0.1; 0 1 -0.3]); 
y = [70;18];

[p,n]           = size(C);
[U,~,V]         = svd(C);
V1              = V(:,1:p);
V2              = V(:,p+1:n);
U1              = U;

c_k = [26;18;40];
G_k = [1    0    0.1;
       0   2.1   2  ;
       0    0     1];
R_k = zonotope(c_k, G_k);

r_R_k = radius(R_k);

K = 50;
M = r_R_k + norm(V2'*c_k);
% M = K;

Z_x_y = measurement_zonotope_M(y,C,Z_v,M);

figure(2);
clf;
plot(Z_x_y, [2 3], 'k-');
hold on
plot(R_k, [2 3], 'b-');
h1 = ellipse(K,K,2*pi,0,0,'r');
h2 = ellipse(r_R_k,r_R_k,2*pi,c_k(2),c_k(3),'g-');
legend("$\mathcal{Z}_{x|y}$","$\mathcal{R}_k$","Interpreter","latex");
