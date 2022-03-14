function M_x_y_i = measurement_zonotope_matzono(Y,C,M_v)
    
    [p,n] = size(C);
    
    % extract zonotope info 
    C_v = M_v.center;
    G_v = M_v.generator;
    
    [U,Sigma,V]     = svd(C);
    V1              = V(:,1:p);
    V2              = V(:,p+1:n);
    U1              = U;
    
    Sigma_sq_inv    = inv(Sigma(1:p,1:p));
    
    c_x             = V1*Sigma_sq_inv*U1'*(Y-C_v);
    M               = 100;
    G_x_null        = V2*M;
    
    
    for i=1:length(G_v)
        G_x{i}       = V1*Sigma_sq_inv*U1'*G_v{i};
        % G_x{i}             = [G_x_range{i}];
    end
    G_x{i} = repmat(G_x_null,1,size(G_x{1},2));
    
    M_x_y_i         = matZonotope(c_x, G_x);
end