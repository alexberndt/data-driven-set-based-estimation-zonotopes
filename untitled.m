% test

clc
clear

c = [0;0];
G = [1 0 1; 0 1 1];
A = [1 0 0];
b = [0];

% A = [];
% b = [];

conZon = conZonotope(c,G,A,b);

c = [0;0];
G = [0 1; 1 1];
A = [];
b = [];
conZon2 = conZonotope(c,G,A,b);

figure(1);
plot(conZon, [1 2], 'b-');
plot(conZon2, [1 2], 'r-');

% inter_zon = conZonotope(Z_x_y_i{3,idx});
% 
% for j = 1:q
%     plot(Z_x_y_i{j,idx},[1 2],'g-');
% end
% 
% for j = 1:2
%     inter_zon = and(inter_zon,conZonotope(Z_x_y_i{j,idx}));
% end
% 
% plot(inter_zon,[1 2],'r-');
