function conZres = intersectConZonoZono2(z1,Cl,Z_v,yl,method)

    G_k = z1.Z(:,2:end);
    c_k = z1.center;
    
    n = size(Cl{1},2); % number of states
    
    if strcmp(method,'frobenius') 
        
        C_sizes = cell(1,length(Cl));
        sum_of_p_is = 0;
        for idx = 1:length(Cl)
            C_z_v = Cl{idx};
            [p_i,~] = size(C_z_v);
            C_sizes{idx} = p_i;
            sum_of_p_is = sum_of_p_is + p_i;
        end
        
        % concatenation of lambdas
        lambda0 = rand(length(z1.center),sum_of_p_is);
        options = optimoptions(@fminunc,'Algorithm', 'quasi-newton','Display','off');
        % optimize to obtain the weights
        lambda = fminunc(@fun,lambda0, options);
    else
        disp('Method is not supported');
        return;
    end

    % Prepare center Eq 21.
    cbar_k = c_k;
    
    p_so_far = 1;
    for i=1:length(Cl)
        p = C_sizes{i};
        cbar_k = cbar_k + lambda(:,p_so_far:(p_so_far+p-1))*( yl{i} - Cl{i}*z1.center - Z_v{i}.center);
        p_so_far = p_so_far + p;
    end

    % Prepare generators Eq 22.
    LHS = eye(length(z1.center));
    RHS = zeros(n,sum_of_p_is);
    
    p_so_far = 1;
    for idx = 1:length(Z_v)
        p = C_sizes{idx};
        LHS = LHS - lambda(:,p_so_far:(p_so_far+p-1))*Cl{idx};
        
        G_v_i = Z_v{idx}.generators;
        
        RHS(:, p_so_far : (p_so_far+p-1) ) =- lambda( :, p_so_far: (p_so_far+p-1) ) * G_v_i;
        p_so_far = p_so_far + p;
    end

    LHS = LHS * G_k;
    Gbar_k = [LHS RHS];
    Zres = zonotope([cbar_k Gbar_k]);

    % determine A,b
    % Abar_k = [];

    % for i=1:length(Z_v)
    %     G_v_i = Z_v{i}.generators;
    %     C_k_i = Cl{i};
    %     
    %     C_k_i*G_v_i G_v_i 
    % end

    b_height = length(z1.b) + sum_of_p_is;
    bbar_k = [z1.b];
    current_b_size = length(z1.b);
    bbar_k = zeros(b_height,1);
    bbar_k(1:current_b_size) = z1.b;
    
    p_so_far = 0;
    for i = 1:length(Z_v)
        
        b_comp = yl{i} - Cl{i}*z1.center + Z_v{i}.center;
        p = C_sizes{i};
        j_s = current_b_size + 1 + p_so_far;
        j_f = current_b_size + 1 + p_so_far + p - 1;
        
        p_so_far = p_so_far + p;
        
        bbar_k(j_s:j_f) = b_comp;
    end

    [~,cols] = size(Gbar_k)

    current_A_size = size(z1.A,1);
    [A_row, A_col] = size(z1.A);

    A_height = current_A_size + sum_of_p_is;

    Abar_k = zeros(A_height,cols);

    if A_row > 0 && A_col > 0
        Abar_k(1:A_row,1:A_col) = z1.A;
    end

    p_so_far = 0;
    width_so_far = 0;

    width_G = size(G_k,2);

    for i=1:length(Cl)
        p = C_sizes{i};

        j_s = current_A_size + 1 + p_so_far;
        j_f = current_A_size + 1 + p_so_far + p - 1;
        Abar_k( j_s:j_f, 1:width_G ) = Cl{i}*G_k;
        Abar_k( j_s:j_f, width_G + width_so_far + 1 : width_G + width_so_far + size(Z_v{i}.generators,2) ) = Z_v{i}.generators;

        width_so_far = width_so_far + size(Z_v{i}.generators,2);

        p_so_far = p_so_far + p;
    end

    % disp('hello');
    conZres = conZonotope(Zres.center, Zres.generators, Abar_k, bbar_k);
%     disp("conZon done");

    function nfro = fun(lambda)
        part1 = eye(length(z1.center));
        
        p_so_far_f = 1;
        for ii=1:length(Z_v)
            p_f = C_sizes{ii};
            part1 = part1 - lambda(:,ii:(ii+p_f-1))*Cl{ii};
            G_v_i = Z_v{ii}.generators;
            
            
            part2(:,p_so_far_f:(p_so_far_f+p_f-1)) = G_v_i*lambda(:,p_so_far_f:(p_so_far_f+p_f-1));
            
            p_so_far_f = p_so_far_f + p_f;
        end
        part1 = part1 * G_k;
        Gbar_k = [part1 part2];
        if strcmp(method,'frobenius')
            nfro = trace(Gbar_k*Gbar_k');
        end
        
    end


end
