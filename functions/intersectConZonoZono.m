function conZres = intersectConZonoZono(z1,Cl,Z_v,yl,varargin)

if nargin==4
    %The optimization function is based on norm of the generators
    method='normGen';
elseif nargin==5
    method =varargin{1};
end
%     H = generators(z1);
    H = z1.Z(:,2:end);
if strcmp(method,'svd') || strcmp(method,'radius') || strcmp(method,'frobenius') 
%     lambda0 = zeros(length(z1.center),length(Z_v));
    
%     lambda_len = length(z1.center);
    
    C_sizes = cell(1,length(Cl));
    
    sum_of_p = 0;
    for idx = 1:length(Cl)
        C_z_v = Cl{idx};
        [p_idx,~] = size(C_z_v);
        C_sizes{idx} = p_idx;
        sum_of_p = sum_of_p + p_idx;
    end
    
    lambda0 = zeros(length(z1.center),sum_of_p);
        
    options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
    %find the weights
    lambda = fminunc(@fun,lambda0, options);
elseif strcmp(method,'normGen')
    % Find the analytical solution
    h_combined=[];
    for i=1:length(Cl)
        h_combined = [ h_combined ; Cl{i}];
    end    
    gamma=eye(length(Cl));
    num= H*H'*h_combined';
    den = h_combined * H*H' * h_combined' ;
    for i=1:length(Cl)
        G_v_i = Z_v{i}.generators;
        den = den + gamma(:,i) *G_v_i^2* gamma(:,i)';
    end
    
    lambda = num * den^-1;
else
    disp('Method is not supported');
    return;
end

% get dimensions

%prepare center
c_new=z1.center;
for i=1:length(Cl)
    p = C_sizes{i};
    c_new = c_new + lambda(:,i:(i+p-1))*( yl{i} - Cl{i}*z1.center - Z_v{i}.center);
end

%prepare generators
part1 = eye(length(z1.center));
for ii=1:length(Z_v)
    p = C_sizes{ii};
    part1 = part1 - lambda(:,ii:(ii+p-1))*Cl{ii};
    G_v_i = Z_v{ii}.generators;
    part2(:,ii:(ii+p-1)) = G_v_i*lambda(:,ii:(ii+p-1));
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

bbar_k = [z1.b];
for i=1:length(Z_v)
    bbar_k = [bbar_k;
    yl{i} - Cl{i}*z1.center + Z_v{i}.center];
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
conZres = conZonotope(Zres.center, Zres.generators ,Abar_k, bbar_k);
disp("conZon done");

    function nfro = fun(lambda)
        part1 = eye(length(z1.center));
        for ii=1:length(Z_v)
            p = C_sizes{ii};
            part1 = part1 - lambda(:,ii:(ii+p-1))*Cl{ii};
            G_v_i = Z_v{ii}.generators;
            part2(:,ii:(ii+p-1)) = G_v_i*lambda(:,ii:(ii+p-1));
        end
        part1 = part1 * H;
        H_new = [part1 part2];
        if strcmp(method,'svd')
            nfro = sum(svd(H_new));
        elseif strcmp(method,'radius')
            nfro = radius(zonotope([zeros(length(z1.center),1) H_new]));
        elseif strcmp(method,'frobenius')
            nfro = trace(H_new*H_new');
        end
        
    end


end
