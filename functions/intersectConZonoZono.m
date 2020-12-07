function conZres = intersectConZonoZono(z1,Cl,Z_v,yl,varargin)

if nargin==4
    %The optimization function is based on norm of the generators
    method='normGen';
elseif nargin==5
    method =varargin{1};
end
%     H = generators(z1);
    H = z1.Z(:,2:end);
if strcmp(method,'frobenius') 
    
    C_sizes = cell(1,length(Cl));
    sum_of_p = 0;
    for idx = 1:length(Cl)
        C_z_v = Cl{idx};
        [p_idx,~] = size(C_z_v);
        C_sizes{idx} = p_idx;
        sum_of_p = sum_of_p + p_idx;
    end
    
    lambda0 = rand(length(z1.center),sum_of_p);
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
    %find the weights
    lambda = lambda0; %fminunc(@fun,lambda0, options);
else
    disp('Method is not supported');
    return;
end

%prepare center Eq 21.
c_new=z1.center;
for i=1:length(Cl)
    p = C_sizes{i};
    c_new = c_new + lambda(:,i:(i+p-1))*( yl{i} - Cl{i}*z1.center - Z_v{i}.center);
end

%prepare generators Eq 22.
part1 = eye(length(z1.center));
p_so_far = 1;
for ii=1:length(Z_v)
    p = C_sizes{ii};
    part1 = part1 - lambda(:,p_so_far:(p_so_far+p-1))*Cl{ii};
    G_v_i = Z_v{ii}.generators;
    part2(:,p_so_far:(p_so_far+p-1)) = G_v_i*lambda(:,p_so_far:(p_so_far+p-1));
    p_so_far = p_so_far + p;
end

part1 = part1 * H;
H_new = [part1 part2];
Zres = zonotope([c_new H_new]);

% determine A,b
% Abar_k = [];

% for i=1:length(Z_v)
%     G_v_i = Z_v{i}.generators;
%     C_k_i = Cl{i};
%     
%     C_k_i*G_v_i G_v_i 
% end

b_height = length(z1.b) + sum_of_p;
bbar_k = [z1.b];
current_b_size = length(z1.b);
bbar_k = zeros(b_height,1);
bbar_k(1:current_b_size) = z1.b;
p_so_far = 0;

for i=1:length(Z_v)
%     bbar_k = [bbar_k;
    b_comp = yl{i} - Cl{i}*z1.center + Z_v{i}.center;
        
    p = C_sizes{i};
    
    j_s = current_b_size + 1 + p_so_far;
    j_f = current_b_size + 1 + p_so_far + p - 1;
    
    bbar_k(j_s:j_f) = b_comp;
end

[~,cols] = size(H_new);

current_A_size = size(z1.A,1);
[A_row, A_col] = size(z1.A);

A_height = current_A_size + sum_of_p;

Abar_k = zeros(A_height,cols);

if A_row > 0 && A_col > 0
    Abar_k(1:A_row,1:A_col) = z1.A;
end

p_so_far = 0;
width_so_far = 0;

width_G = size(H,2);

for i=1:length(Cl)
    p = C_sizes{i};
    
    j_s = current_A_size + 1 + p_so_far;
    j_f = current_A_size + 1 + p_so_far + p - 1;
    Abar_k( j_s:j_f, 1:width_G ) = Cl{i}*H;
    Abar_k( j_s:j_f, width_G + width_so_far + 1 : width_G + width_so_far + size(Z_v{i}.generators,2) ) = -Z_v{i}.generators;
    
    width_so_far = width_so_far + size(Z_v{i}.generators,2);
    
    p_so_far = p_so_far + p;
end

% disp('hello');
conZres = conZonotope(Zres.center, Zres.generators, Abar_k, bbar_k);
disp("conZon done");

    function nfro = fun(lambda)
        part1 = eye(length(z1.center));
        
        p_so_far = 1;
        for ii=1:length(Z_v)
            p = C_sizes{ii};
            part1 = part1 - lambda(:,ii:(ii+p-1))*Cl{ii};
            G_v_i = Z_v{ii}.generators;
            
            
            part2(:,p_so_far:(p_so_far+p-1)) = G_v_i*lambda(:,p_so_far:(p_so_far+p-1));
            
            p_so_far = p_so_far + p;
        end
        part1 = part1 * H;
        H_new = [part1 part2];
        if strcmp(method,'frobenius')
            nfro = trace(H_new*H_new');
        end
        
    end


end
