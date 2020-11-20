%% Compare strips, zonotopes
%
%
%
%% init
clc
clear

% define observability matrices
q = 3;
C       = cell(1,q);
C{1}    = [1 0];
C{2}    = [0 1];
C{3}    = [1 1];
Z_v_meas     = cell(1,q);
Z_v_meas{1}  = zonotope(-5,2);
Z_v_meas{2}  = zonotope(0,2);
Z_v_meas{3}  = zonotope(0,2);

% define initial zonotope representing R_k
Z = zonotope([1 2 2 2 6 2 8;3 2 2 0 5 0 6 ]);
Z_svd = Z;

% measurement values
yl{1} = -3;
yl{2} = 2;
yl{3} = 2;

%% UNCOMMENT FOR 2x2 C matrix

% C{3}    = [0.7071 -0.7071;
%            0 1];
% Z_v_meas{3}  = zonotope([0;0],blkdiag(3,3));
% yl{3} = [2; 2];

%%

% calculate Z_svd using SVD approach 
Z_x_y_i = cell(1,q);
for i = 1:q
    Z_x_y_i{i}   = measurement_zonotope(yl{i}, C{i}, Z_v_meas{i});
    Z_svd = and(Z_svd,Z_x_y_i{i});
end

Z_zonoZono = intersectZonoZono(Z,C,Z_v_meas,yl,'radius');

% andAveraging
Z_x_y_i{4} = Z;
Z_svd_andavg = andAveraging1(Z_x_y_i,'radius'); %,true,1.0);

%% Plot the results
figure(10); 
clf;
hold on 
plot(Z,[1 2],'g-+');
% plot(res_zono,[1 2],'b-*');
plot(Z_svd,[1 2],'r-*');
% plot(poly,[1 2],'k-*');
plot(Z_svd_andavg,[1 2],'m-*');

% plot Zonotope-strips
for i = 1:q
    plot(Z_x_y_i{i},[1 2],'g-');
end
plot(Z_zonoZono,[1 2],'k--*');
xlim([-10 10]);
ylim([-10 10]);
legend('R_k zonotope, zonotope-strips','SVD <and> approach','intersectZonoZono','SVD and Averaging');









































































%.