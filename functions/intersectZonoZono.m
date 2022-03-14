function Zres = intersectZonoZono(z1,Cl,Z_v,yl,varargin)
% intersectZonoStrip - computes the intersection between one zonotope and
%    list of strips according to [1]
%    the strip is defined as | hx-y | <= d
%
% Syntax:  
%    Zres = intersectZonoStrip(z1,hl,Rl,yl,varargin)
%
% Inputs:
%    z1 - zonotope object
%    h1 - 
%    R1 - 
%    y1 -
%    varargin - methods to calculate the weights
%               'normGen' default and has analytical solution
%               'svd'
%               'radius'
%
% Outputs:
%    res - boolean whether obj is contained in Z, or not
%
% Example: (three strips and one zonotope)
%    hl{1} = [1 0];
%    Rl{1} = 5;
%    yl{1} = -2;
% 
%    hl{2} = [0 1];
%    Rl{2} = 3;
%    yl{2} = 2;
% 
%    hl{3} = [1 1];
%    Rl{3} = 3;
%    yl{3} = 2;
% 
%    Z = zonotope([1 2 2 2 6 2 8;1 2 2 0 5 0 6 ]);
%    res_zono = intersectZonoStrip(Z,hl,Rl,yl);
% 
%    % just for comparison
%    poly = mptPolytope([1 0;-1 0; 0 1;0 -1; 1 1;-1 -1],[3;7;5;1;5;1]);
%    Zpoly = Z & poly;
% 
%    figure; hold on 
%    plot(Z,[1 2],'r-+');
%    plot(poly,[1 2],'r-*');
%    plot(Zpoly,[1 2],'b-+');
%    plot(res_zono,[1 2],'b-*');
% 
%    legend('zonotope','strips','zono&poly','zonoStrips');
%
%
% References:
%    [1] Amr Alanwar, Jagat Jyoti Rath, Hazem Said, Matthias Althoff
%       Distributed Set-Based Observers Using Diffusion Strategy
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Amr Alanwar
% Written:       09-Mar-2020
% Last update:   ---              
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin==4
    %The optimization function is based on norm of the generators
    method='normGen';
elseif nargin==5
    method =varargin{1};
end
    H = generators(z1);
if strcmp(method,'svd') || strcmp(method,'radius') || strcmp(method,'frobenius') 
    lambda0 = zeros(length(z1.center),length(Z_v));
    
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
    part2(:,ii:(ii+p-1)) = - G_v_i*lambda(:,ii:(ii+p-1));
end

part1 = part1 * H;
H_new = [part1 part2];
Zres = zonotope([c_new H_new]);



    function nfro = fun(lambda)
        part1 = eye(length(z1.center));
        for ii=1:length(Z_v)
            p = C_sizes{ii};
            part1 = part1 - lambda(:,ii:(ii+p-1))*Cl{ii};
            G_v_i = Z_v{ii}.generators;
            part2(:,ii:(ii+p-1)) = - G_v_i*lambda(:,ii:(ii+p-1));
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
