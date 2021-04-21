function Z_x_y_i = measurement_zonotope_M(y,C,Z_v,M)
    
    [p,n] = size(C);
    r = rank(C);
    
    % extract zonotope info 
    c_v = Z_v.center;
    G_v = Z_v.generators;
    
    [U,Sigma,V]     = svd(C);
    Sigma_sq_inv    = inv(Sigma(1:r,1:r));
    
    if (r == n) && (r < p)
        % skinny tall C
        disp("skinny tall C");
        U1 = U(:,1:r);
        U2 = U(:,r+1:p);
        V1 = V;
        V2 = [];
    elseif (r < n) && (r == p)
        % wide short C
        disp("wide short C");
        U1 = U;
        U2 = [];
        V1 = V(:,1:r);
        V2 = V(:,r+1:n);
    elseif (r < n) && (r < p)
        % rank-deficient C
        disp("rank-deficient C");
        U1 = U(:,1:r);
        U2 = U(:,r+1:p);
        V1 = V(:,1:r);
        V2 = V(:,r+1:n);
    elseif (r == n) && (r == p)
        % square invertible C
        disp("invertible C");
        U1 = U;
        U2 = [];
        V1 = V;
        V2 = [];
    else
        disp("Error");
    end
    
    c_x             = V1*Sigma_sq_inv*U1'*(y-c_v);
    G_x_range       = V1*Sigma_sq_inv*U1'*G_v;

    G_x_null        = V2*M;
    
    G_x             = [G_x_range, G_x_null];

    Z_x_y_i         = zonotope(c_x, G_x);
end