function [L, D] = myldlt(A)
    n = size(A, 1)
    L = zeros(n, n)
    D = zeros(n, n)
    
    D(1, 1) = A(1, 1)
    L(1, 1) = 1
    
    for i = 2 : n
        for j = 1 : i - 1
            somme = 0
            for k = 1 : j - 1
                somme = somme + L(i, k) * L(j, k) * D(k, k)
            end
            L(i, j) = (1 / D(j , j)) * (A(i, j) - somme )
        end

        somme = 0
        for k = 1 : i - 1
            somme = somme + (L(i, k) * L(i, k)) * D(k, k)
        end
        D(i, i) = A(i, i) - somme

        L(i, i) = 1
    end
end


function [L, U] = mylu(A)
    [n, m] = size(A)
    for k = 1 : n - 1
        for i = k + 1 : n
            A(i, k) = A(i, k) / A(k, k)
        end
        for i = k + 1 : n
            for j = k + 1 : n
                A(i, j) = A(i, j) - A(i, k) * A(k, j)
            end
        end
    end
    
    U = triu(A)
    L = tril(A)
    
    for i = 1 : n
        L(i, i) = 1
    end
end



function [A] = symetrise(A)
    n = size(A, 1)
    for i = 1 : n
        for j = i : n
            A(i, j) = A(j, i)
        end
    end
end
