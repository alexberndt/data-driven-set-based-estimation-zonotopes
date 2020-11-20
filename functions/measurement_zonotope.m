function Z_x_y_i = measurement_zonotope(y,C,Z_v)
    
    [p,n] = size(C);
    
    % extract zonotope info 
    c_v = Z_v.center;
    G_v = Z_v.generators;
    
    [U,Sigma,V]     = svd(C);
    V1              = V(:,1:p);
    V2              = V(:,p+1:n);
    U1              = U;
    
    Sigma_sq_inv    = inv(Sigma(1:p,1:p));
    
    c_x             = V1*Sigma_sq_inv*U1'*(y-c_v);
    G_x_range       = V1*Sigma_sq_inv*U1'*G_v;

    M               = 100;
    G_x_null        = V2*M;
    
    G_x             = [G_x_range, G_x_null];

    Z_x_y_i         = zonotope(c_x, G_x);
end