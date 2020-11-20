function res = check_in_interval(A_Matrix, Interval_Matrix)
    % Function that checks if a matrix A_Matrix 
    %  is inside an intervale matrix Interval_Matrix
    %
    % Input:
    %  A_Matrix
    %  Interval_Matrix
    %
    % Output:
    %  res : boolean - true if A_Matrix is in Interval_Matrix 


    [n,m] = size(A_Matrix);

    res = true;
    
    Inf_Mat = Interval_Matrix.Inf;
    Sup_Mat = Interval_Matrix.Sup;
    
    for i = 1:n
        for j = 1:m
            
            if A_Matrix(i,j) <= Inf_Mat(i,j) ||  A_Matrix(i,j) >= Sup_Mat(i,j)
                res = false;
                break;
            end
        end
    end
end